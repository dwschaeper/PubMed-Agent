from Bio import Entrez
from langchain.docstore.document import Document


def search_pubmed(email, query, max_results=5):
    Entrez.email = email

    handle = Entrez.esearch(db='pubmed', term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record['IdList']

    # Fetch XML data
    fetch_handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="xml", retmode="xml")
    fetch_records = Entrez.read(fetch_handle)

    # collect the title and abstract
    documents = []
    for article in fetch_records["PubmedArticle"]:
        pmid = article["MedlineCitation"]["PMID"]

        title = article["MedlineCitation"]["Article"]["ArticleTitle"]

        abstract_parts = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [])
        abstract = " ".join(str(part) for part in abstract_parts)

        date = 'unknown'
        article_date = article["MedlineCitation"]["Article"].get("ArticleDate", [])
        if article_date:
            year = article_date[0].get("Year", "")
            month = article_date[0].get("Month", "")
            day = article_date[0].get("Day", "")
            date = f"{year}-{month.zfill(2)}-{day.zfill(2)}"
        else:
            journal_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
            date = journal_date.get("Year", "unknown")

        doc = Document(page_content=f"Title: {title}\nAbstract: {abstract}\nPMID: {pmid}",
                       metadata={"title": title, "date": date, "PMID": str(pmid)}
                       )
        documents.append(doc)

    return documents