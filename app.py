import streamlit as st
from Bio import Entrez
from groq import Groq

# Configurare pagină
st.set_page_config(page_title="Medical Assistant & Links", page_icon="🩺", layout="wide")

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

def get_pubmed_data(query, email, max_results=5):
    """
    Caută pe PubMed și returnează două lucruri:
    1. Textul pentru AI (String)
    2. Lista de linkuri și titluri pentru UI (Listă)
    """
    Entrez.email = email
    try:
        # Pas 1: Căutare ID-uri
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        st.write(f"ℹ️ Găsit {len(id_list)} studii.")

        if not id_list:
            return None, None

        # Pas 2: Descărcare detalii în format XML (mai ușor de procesat)
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        papers = Entrez.read(handle)
        handle.close()

        # Procesare date
        ai_context = ""
        ui_references = []

        for i, paper in enumerate(papers['PubmedArticle']):
            try:
                # Extragem datele
                medline = paper['MedlineCitation']
                article = medline['Article']
                
                pmid = str(medline['PMID'])
                title = article['ArticleTitle']
                
                # Încercăm să luăm abstractul (unele studii nu au abstract)
                try:
                    abstract_list = article['Abstract']['AbstractText']
                    abstract = " ".join(abstract_list)
                except KeyError:
                    abstract = "Abstract indisponibil."

                # Construim linkul
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

                # Pregătim textul pentru AI
                ai_context += f"Studiu [{i+1}]:\nTitlu: {title}\nID: {pmid}\nAbstract: {abstract}\n\n"
                
                # Salvăm pentru afișare
                ui_references.append({"index": i+1, "title": title, "url": link})

            except Exception as e:
                continue # Dacă un articol are format ciudat, îl sărim

        return ai_context, ui_references

    except Exception as e:
        st.error(f"Eroare PubMed: {e}")
        return None, None

def generate_answer(query, context):
    prompt = f"""
    Ești un asistent medical expert. 
    Întrebarea utilizatorului: {query}
    
    Mai jos ai o listă de studii (marcate cu Studiu [1], [2] etc.).
    
    Instrucțiuni:
    1. Răspunde în LIMBA ROMÂNĂ.
    2. Sintetizează informația medicală.
    3. Când folosești o informație, pune referința în text folosind paranteze pătrate, ex: [1], [2].
    
    Context:
    {context}
    """
    
    try:
        completion = client.chat.completions.create(
            model="llama-3.3-70b-versatile", 
            messages=[
                {"role": "system", "content": "Ești un medic cercetător care răspunde în română și citează sursele cu numere [1]."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3,
        )
        return completion.choices[0].message.content
    except Exception as e:
        return f"Eroare AI: {e}"

# --- 3. INTERFAȚA ---
st.title("🩺 PubMed AI Assistant + Surse")
st.markdown("Caută studii, generează o sinteză și oferă linkuri directe către articole.")

query = st.text_input("Termen de căutare (Engleză):", placeholder="ex: vitamin d deficiency symptoms")

if st.button("Caută"):
    if not query:
        st.warning("Scrie o întrebare.")
    else:
        with st.spinner("Căutăm și procesăm linkurile..."):
            context_text, references = get_pubmed_data(query, email_address)
            
        if context_text:
            # 1. Generare Răspuns
            with st.spinner("Generăm sinteza..."):
                ans = generate_answer(query, context_text)
                st.markdown("### 📝 Răspuns Sintetizat:")
                st.write(ans)
            
            st.divider()
            
            # 2. Afișare Linkuri
            st.markdown("### 🔗 Bibliografie și Linkuri:")
            for ref in references:
                st.markdown(f"**[{ref['index']}]** [{ref['title']}]({ref['url']})")
                
            with st.expander("Vezi textul brut trimis la AI"):
                st.text(context_text)
        else:
            st.error("Nu am găsit studii.")
