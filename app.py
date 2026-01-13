import streamlit as st
from Bio import Entrez
from groq import Groq

# Configurare pagină
st.set_page_config(page_title="Asistent Medical Complet", page_icon="🩺", layout="wide")

# --- 1. SECRETS ---
try:
    if "GROQ_API_KEY" in st.secrets:
        api_key = st.secrets["GROQ_API_KEY"]
    else:
        st.error("Lipsește GROQ_API_KEY.")
        st.stop()
        
    if "EMAIL_ADRESS" in st.secrets:
        email_address = st.secrets["EMAIL_ADRESS"]
    else:
        email_address = st.secrets.get("EMAIL_ADDRESS", "email@test.com")
except FileNotFoundError:
    st.error("Configurează secrets.toml!")
    st.stop()

client = Groq(api_key=api_key)

# --- 2. FUNCȚII ---

def translate_to_english(text):
    """Traduce întrebarea din Română în Engleză pentru PubMed"""
    try:
        completion = client.chat.completions.create(
            model="llama-3.3-70b-versatile", # Modelul NOU și funcțional
            messages=[
                {"role": "system", "content": "You are a medical translator. Translate the user's query from Romanian to English keywords suitable for a PubMed search. Return ONLY the English keywords, no explanation."},
                {"role": "user", "content": text}
            ],
            temperature=0,
        )
        return completion.choices[0].message.content.strip()
    except Exception as e:
        st.error(f"Eroare traducere: {e}")
        return text

def get_pubmed_data(query, email, max_results=5):
    """Caută, extrage textul și creează linkurile"""
    Entrez.email = email
    try:
        # Căutare ID-uri
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        st.write(f"ℹ️ PubMed a găsit {len(id_list)} studii.")

        if not id_list:
            return None, None

        # Descărcare XML
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        papers = Entrez.read(handle)
        handle.close()

        ai_context = ""
        ui_references = []

        for i, paper in enumerate(papers['PubmedArticle']):
            try:
                medline = paper['MedlineCitation']
                article = medline['Article']
                
                pmid = str(medline['PMID'])
                title = article['ArticleTitle']
                
                try:
                    abstract_list = article['Abstract']['AbstractText']
                    abstract = " ".join(abstract_list)
                except KeyError:
                    abstract = "Abstract indisponibil."

                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

                # Context pentru AI
                ai_context += f"Studiu [{i+1}]:\nTitlu: {title}\nID: {pmid}\nAbstract: {abstract}\n\n"
                
                # Date pentru UI (Linkuri)
                ui_references.append({"index": i+1, "title": title, "url": link})

            except Exception as e:
                continue

        return ai_context, ui_references

    except Exception as e:
        st.error(f"Eroare PubMed: {e}")
        return None, None

def generate_answer(query, context):
    prompt = f"""
    Ești un asistent medical expert. 
    Întrebarea utilizatorului (Română): {query}
    
    Mai jos ai o listă de studii (marcate cu Studiu [1], [2] etc.).
    
    Instrucțiuni:
    1. Răspunde în LIMBA ROMÂNĂ.
    2. Sintetizează informația medicală clar și concis.
    3. Când folosești o informație, pune referința în text folosind paranteze pătrate, ex: [1], [2].
    
    Context:
    {context}
    """
    
    try:
        completion = client.chat.completions.create(
            model="llama-3.3-70b-versatile", 
            messages=[
                {"role": "system", "content": "Ești un medic cercetător care răspunde în română și citează sursele."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3,
        )
        return completion.choices[0].message.content
    except Exception as e:
        return f"Eroare AI: {e}"

# --- 3. INTERFAȚA ---
st.title("🩺 PubMed AI Assistant (RO)")
st.markdown("Introduceți o întrebare medicală în **Română**. Sistemul o traduce, caută studii și oferă răspunsul cu linkuri.")

query_ro = st.text_input("Întrebare:", placeholder="ex: care sunt efectele secundare ale imunoterapiei?")

if st.button("Caută Răspuns"):
    if not query_ro:
        st.warning("Te rog scrie o întrebare.")
    else:
        # 1. Traducere
        with st.spinner("Traducem întrebarea..."):
            query_en = translate_to_english(query_ro)
            st.info(f"Termeni căutare generați (Engleză): **{query_en}**")
        
        # 2. Căutare PubMed
        with st.spinner("Căutăm studii și generăm linkuri..."):
            context_text, references = get_pubmed_data(query_en, email_address)
            
        if context_text:
            # 3. Generare Răspuns
            with st.spinner("Llama 3 sintetizează informația..."):
                ans = generate_answer(query_ro, context_text)
                
                st.markdown("### 📝 Răspuns Sintetizat:")
                st.write(ans)
            
            st.divider()
            
            # 4. Afișare Linkuri
            st.markdown("### 🔗 Bibliografie (Click pentru articol):")
            for ref in references:
                st.markdown(f"**[{ref['index']}]** [{ref['title']}]({ref['url']})")
                
            with st.expander("Vezi textul abstractelor (pentru verificare)"):
                st.text(context_text)
        else:
            st.error("Nu am găsit studii relevante pe PubMed.")
