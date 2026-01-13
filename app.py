import streamlit as st
from Bio import Entrez
import openai
import pandas as pd

# Configurare pagină
st.set_page_config(page_title="Asistent Medical PubMed", page_icon="🩺", layout="wide")

# Titlu și descriere
st.title("🩺 Asistent Medical AI - PubMed Search")
st.markdown("""
Acest asistent caută cele mai recente studii pe **PubMed** și folosește AI pentru a sintetiza informația.
**Atenție:** Acest instrument este doar pentru informare și nu înlocuiește judecata clinică profesională.
""")

# Sidebar pentru setări
with st.sidebar:
    st.header("Setări")
    api_key = st.text_input("Introdu cheia OpenAI API", type="password")
    email = st.text_input("Email (cerut de PubMed)", placeholder="doctor@exemplu.com")
    max_results = st.slider("Număr de studii de analizat", 1, 10, 5)
    st.markdown("---")
    st.markdown("Obține o cheie API de la [OpenAI Platform](https://platform.openai.com/).")

# Funcția de căutare pe PubMed
def search_pubmed(query, email, max_results=5):
    Entrez.email = email
    try:
        # 1. Căutare ID-uri
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            return None

        # 2. Descărcare detalii studii
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        articles_text = handle.read()
        handle.close()
        
        return articles_text
    except Exception as e:
        st.error(f"Eroare la conectarea cu PubMed: {e}")
        return None

# Funcția AI (GPT)
def generate_answer(query, context, api_key):
    client = openai.OpenAI(api_key=api_key)
    
    prompt = f"""
    Ești un asistent medical expert. Folosește DOAR următoarele rezumate din studii științifice pentru a răspunde la întrebarea medicului.
    Dacă informația nu există în studii, spune asta. Citează studiile (Autor, An) când este posibil.

    Întrebare: {query}

    Studii PubMed (Context):
    {context}
    """

    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo", # Poți schimba cu gpt-4 dacă ai acces
            messages=[
                {"role": "system", "content": "Ești un asistent de cercetare medicală."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Eroare AI: {e}"

# Interfața principală
query = st.text_input("Ce informație medicală căutați?", placeholder="ex: Tratamentul actual pentru hipertensiune rezistentă")

if st.button("Caută și Analizează"):
    if not api_key:
        st.warning("Te rog introdu cheia OpenAI API în meniul din stânga.")
    elif not email:
        st.warning("Te rog introdu o adresă de email pentru PubMed în meniul din stânga.")
    elif not query:
        st.warning("Te rog introdu o întrebare.")
    else:
        with st.spinner("Căutăm pe PubMed..."):
            pubmed_data = search_pubmed(query, email, max_results)
        
        if pubmed_data:
            with st.expander("Vezi datele brute (Abstracte PubMed)"):
                st.text(pubmed_data)
            
            with st.spinner("AI-ul analizează studiile..."):
                answer = generate_answer(query, pubmed_data, api_key)
                st.success("Analiză Finalizată")
                st.markdown("### Răspuns Sintetizat:")
                st.write(answer)
        else:
            st.error("Nu s-au găsit articole pe PubMed pentru această căutare.")
