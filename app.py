import streamlit as st
from Bio import Entrez
import openai

# Configurare pagină
st.set_page_config(page_title="Asistent Medical PubMed", page_icon="🩺", layout="wide")

# --- 1. SECRETS MANAGEMENT ---
try:
    # Verificăm dacă cheile există
    if "OPENAI_API_KEY" in st.secrets:
        api_key = st.secrets["OPENAI_API_KEY"]
    else:
        st.error("Lipsește OPENAI_API_KEY din secrets.")
        st.stop()
        
    if "EMAIL_ADRESS" in st.secrets:
        email_address = st.secrets["EMAIL_ADRESS"]
    else:
        # Fallback dacă ai scris greșit ADDRESS sau ADRESS
        email_address = st.secrets.get("EMAIL_ADDRESS", "email_generic@test.com")
        
except FileNotFoundError:
    st.error("Fișierul secrets.toml nu a fost găsit! (Local)")
    st.stop()

# --- 2. FUNCȚII ---

def translate_to_english(text, api_key):
    """Traduce întrebarea în engleză pentru PubMed folosind GPT"""
    client = openai.OpenAI(api_key=api_key)
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "You are a translator. Translate the following medical query to English keywords suitable for PubMed search. Return ONLY the English keywords."},
                {"role": "user", "content": text}
            ]
        )
        return response.choices[0].message.content
    except Exception as e:
        st.error(f"Eroare la traducere: {e}")
        return text # Returnăm textul original dacă eșuează

def search_pubmed(query, email, max_results=5):
    """Caută pe PubMed"""
    Entrez.email = email
    
    try:
        # Pasul A: Căutare ID-uri
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        
        # DEBUG: Vedem câți am găsit
        st.write(f"ℹ️ PubMed a găsit {len(id_list)} articole pentru termenul: *{query}*")

        if not id_list:
            return None

        # Pasul B: Descărcare conținut
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        articles_text = handle.read()
        handle.close()
        return articles_text
        
    except Exception as e:
        st.error(f"⚠️ Eroare critică PubMed: {e}")
        return None

def generate_answer(query, context, api_key):
    """Generează răspunsul final"""
    client = openai.OpenAI(api_key=api_key)
    
    prompt = f"""
    Ești un asistent medical expert. Răspunde la întrebarea utilizatorului în LIMBA ROMÂNĂ.
    Folosește informațiile din rezumatele de mai jos. Citează sursele (Autor, An).
    
    Întrebare originală: {query}
    
    Context (Studii PubMed):
    {context}
    """

    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "Ești un medic specialist care răspunde în limba română."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Eroare AI: {e}"

# --- 3. INTERFAȚA ---
st.title("🩺 Asistent Medical AI")
st.markdown("Scrie întrebarea în **Română**. AI-ul o va traduce, va căuta studii internaționale și îți va răspunde în Română.")

query = st.text_input("Întrebare:", placeholder="ex: Care sunt riscurile aspirinei la copii?")

if st.button("Caută Răspuns"):
    if not query:
        st.warning("Te rog scrie o întrebare.")
    else:
        # 1. Traducem
        with st.spinner("Traducem întrebarea pentru PubMed..."):
            english_query = translate_to_english(query, api_key)
            st.caption(f"Termeni căutare (Engleză): {english_query}")
        
        # 2. Căutăm
        with st.spinner("Căutăm studii științifice..."):
            pubmed_data = search_pubmed(english_query, email_address)
        
        # 3. Analizăm
        if pubmed_data:
            with st.expander("📄 Vezi datele brute (Abstracte în Engleză)"):
                st.text(pubmed_data)
            
            with st.spinner("🧠 AI-ul sintetizează răspunsul..."):
                answer = generate_answer(query, pubmed_data, api_key)
                st.markdown("### 📝 Răspuns:")
                st.write(answer)
        else:
            st.error("Nu s-au găsit studii relevante. Încearcă să reformulezi întrebarea.")
