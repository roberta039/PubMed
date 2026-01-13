import streamlit as st
from Bio import Entrez
from groq import Groq

# Configurare pagină
st.set_page_config(page_title="Medical Assistant", page_icon="🩺", layout="wide")

# --- 1. SECRETS MANAGEMENT ---
try:
    if "GROQ_API_KEY" in st.secrets:
        api_key = st.secrets["GROQ_API_KEY"]
    else:
        st.error("Lipsește GROQ_API_KEY din secrets.")
        st.stop()
        
    if "EMAIL_ADRESS" in st.secrets:
        email_address = st.secrets["EMAIL_ADRESS"]
    else:
        email_address = st.secrets.get("EMAIL_ADDRESS", "email@test.com")
        
except FileNotFoundError:
    st.error("Configurează secrets.toml!")
    st.stop()

# Inițializare client Groq
client = Groq(api_key=api_key)

# --- 2. FUNCȚII ---

def search_pubmed(query, email, max_results=5):
    """Caută pe PubMed direct cu termenul introdus"""
    Entrez.email = email
    try:
        # Căutare
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        st.write(f"ℹ️ S-au găsit {len(id_list)} studii pentru: *{query}*")

        if not id_list:
            return None

        # Descărcare detalii
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        articles_text = handle.read()
        handle.close()
        return articles_text
    except Exception as e:
        st.error(f"Eroare PubMed: {e}")
        return None

def generate_answer(query, context):
    """Generează răspunsul final folosind Llama 3.3"""
    prompt = f"""
    Ești un asistent medical expert. Sarcina ta este să sintetizezi informația din studiile de mai jos.
    
    Întrebarea utilizatorului: {query}
    
    Context (Studii PubMed):
    {context}
    
    Instrucțiuni:
    1. Răspunde în LIMBA ROMÂNĂ.
    2. Folosește doar informațiile din context.
    3. Dacă studiile nu sunt relevante, spune asta.
    """
    
    try:
        # FOLOSIM NOUA VERSIUNE DE MODEL: llama-3.3-70b-versatile
        completion = client.chat.completions.create(
            model="llama-3.3-70b-versatile", 
            messages=[
                {"role": "system", "content": "Ești un medic cercetător care răspunde în limba română."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3,
        )
        return completion.choices[0].message.content
    except Exception as e:
        return f"Eroare AI: {e}"

# --- 3. INTERFAȚA ---
st.title("🩺 PubMed AI Assistant (Groq Free)")
st.markdown("""
Acest asistent caută pe PubMed și sintetizează rezultatele.
**Sfat:** Pentru cele mai bune rezultate, introduceți termenii de căutare în **Engleză** (ex: *aspirin side effects*), dar AI-ul va răspunde în Română.
""")

query = st.text_input("Termen de căutare (preferabil în Engleză):", placeholder="ex: immunotherapy lung cancer")

if st.button("Caută"):
    if not query:
        st.warning("Scrie o întrebare.")
    else:
        # 1. Căutare directă
        with st.spinner("Căutăm studii pe PubMed..."):
            pubmed_data = search_pubmed(query, email_address)
            
        # 2. Analiză
        if pubmed_data:
            with st.expander("Vezi rezumatele studiilor (Engleză)"):
                st.text(pubmed_data)
                
            with st.spinner("Llama 3.3 analizează datele..."):
                ans = generate_answer(query, pubmed_data)
                st.markdown("### Răspuns Sintetizat (Română):")
                st.write(ans)
        else:
            st.error("Nu am găsit studii. Încearcă să folosești termeni în engleză.")
