from Bio import Entrez

def search_pubmed(email, query, max_results=5):
    Entrez.email = email

    handle = Entrez.esearch(db='pubmed', term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record['IdList']

    # Fetch XML data
    fetch_handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="xml", retmode="xml")
    fetch_records = Entrez.read(fetch_handle)

    # collect the title and abstract
    entries = []
    for article in fetch_records["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]

        abstract_parts = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [])
        abstract = " ".join(str(part) for part in abstract_parts)

        entries.append({
            "title": title,
            "abstract": abstract
        })

    # format to return one data string
    data = "\n\n".join([
        f"Title: {entry['title']}\nAbstract: {entry['abstract']}"
        for entry in entries
    ])

    return data