import streamlit as st
from Bio import Entrez
import openai
import pandas as pd

# Configurare pagină
st.set_page_config(page_title="Asistent Medical PubMed", page_icon="🩺", layout="wide")

# --- SECRETS MANAGEMENT ---
# Încercăm să luăm cheile din secrets
try:
    api_key = st.secrets["OPENAI_API_KEY"]
    email_address = st.secrets["EMAIL_ADRESS"]
except FileNotFoundError:
    st.error("Cheile API nu sunt configurate! Te rog setează secrets.")
    st.stop()
# ---------------------------

st.title("🩺 Asistent Medical AI - PubMed Search")
st.markdown("Căutare automată de studii și sinteză cu AI.")

# Sidebar simplificat (nu mai cerem cheia)
with st.sidebar:
    st.header("Opțiuni Căutare")
    max_results = st.slider("Număr de studii de analizat", 1, 10, 5)
    st.info("Aplicația folosește o cheie API pre-configurată.")

# Funcția de căutare pe PubMed
def search_pubmed(query, email, max_results=5):
    Entrez.email = email
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            return None

        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        articles_text = handle.read()
        handle.close()
        return articles_text
    except Exception as e:
        st.error(f"Eroare PubMed: {e}")
        return None

# Funcția AI
def generate_answer(query, context, api_key):
    client = openai.OpenAI(api_key=api_key)
    
    prompt = f"""
    Ești un asistent medical expert. Răspunde la întrebare folosind DOAR contextul de mai jos.
    Citează autorii și anii studiilor.
    
    Întrebare: {query}
    
    Context (Studii):
    {context}
    """

    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "Ești un asistent util și precis."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Eroare AI: {e}"

# Interfața
query = st.text_input("Întrebare medicală:", placeholder="ex: Managementul diabetului tip 2 la pacienți vârstnici")

if st.button("Caută"):
    if not query:
        st.warning("Introduceți o întrebare.")
    else:
        with st.spinner("Căutăm pe PubMed..."):
            pubmed_data = search_pubmed(query, email_address, max_results)
        
        if pubmed_data:
            with st.expander("Vezi rezumatele studiilor (Raw Data)"):
                st.text(pubmed_data)
            
            with st.spinner("Generăm răspunsul..."):
                answer = generate_answer(query, pubmed_data, api_key)
                st.markdown("### Răspuns Sintetizat:")
                st.write(answer)
        else:
            st.error("Nu s-au găsit studii.")
