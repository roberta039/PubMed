import streamlit as st
from Bio import Entrez
from groq import Groq # Folosim Groq în loc de OpenAI

# Configurare pagină
st.set_page_config(page_title="Asistent Medical Llama3", page_icon="🩺", layout="wide")

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

def translate_to_english(text):
    """Traduce întrebarea în engleză folosind Llama 3"""
    try:
        completion = client.chat.completions.create(
            model="llama3-70b-8192", # Model mare și gratuit
            messages=[
                {"role": "system", "content": "You are a translator. Translate the medical query to English keywords for PubMed. Return ONLY the keywords, nothing else."},
                {"role": "user", "content": text}
            ],
            temperature=0,
        )
        return completion.choices[0].message.content
    except Exception as e:
        st.error(f"Eroare traducere: {e}")
        return text

def search_pubmed(query, email, max_results=5):
    """Caută pe PubMed"""
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

        # Descărcare
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        articles_text = handle.read()
        handle.close()
        return articles_text
    except Exception as e:
        st.error(f"Eroare PubMed: {e}")
        return None

def generate_answer(query, context):
    """Generează răspunsul final"""
    prompt = f"""
    Ești un asistent medical expert. Răspunde la întrebare în LIMBA ROMÂNĂ.
    Folosește DOAR contextul de mai jos. Citează sursele (Autor, An).
    
    Întrebare: {query}
    
    Context (Abstracte PubMed):
    {context}
    """
    
    try:
        completion = client.chat.completions.create(
            model="llama3-70b-8192", # Folosim modelul 70b pentru acuratețe mai mare
            messages=[
                {"role": "system", "content": "Răspunzi mereu în limba română, profesionist."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3,
        )
        return completion.choices[0].message.content
    except Exception as e:
        return f"Eroare AI: {e}"

# --- 3. INTERFAȚA ---
st.title("🩺 Asistent Medical (Llama 3 + PubMed)")
st.markdown("Acest instrument este **gratuit** și folosește modelul Llama 3 pentru a analiza studii medicale.")

query = st.text_input("Întrebare (în Română):", placeholder="ex: Eficiența metforminei în prediabet")

if st.button("Caută"):
    if not query:
        st.warning("Scrie o întrebare.")
    else:
        # 1. Traducere
        with st.spinner("Traducem..."):
            eng_query = translate_to_english(query)
            st.caption(f"Căutăm: {eng_query}")
            
        # 2. Căutare
        with st.spinner("Căutăm studii..."):
            pubmed_data = search_pubmed(eng_query, email_address)
            
        # 3. Analiză
        if pubmed_data:
            with st.expander("Vezi rezumatele studiilor"):
                st.text(pubmed_data)
                
            with st.spinner("Llama 3 analizează datele..."):
                ans = generate_answer(query, pubmed_data)
                st.markdown("### Răspuns:")
                st.write(ans)
        else:
            st.error("Nu am găsit studii.")
