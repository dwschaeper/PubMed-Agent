from .input_processing import retry_entrez_query

from Bio import Entrez
from urllib.parse import quote_plus
from langchain_core.documents import Document

def get_ids(email: str, entrez_query: str, max_results: int = 50):
    """
    Send entrez_query to get the IDs, WebEnv, and QueryKey.

    Args:
        email (str):        e-mail addressed used with the Entrez query
        entrez_query (str): The query sent to Entrez
        max_results (int):  The maxiumum number of abstracts to return from Entrez

    Return:
        list: The list of abstract IDs
        str:  The WebEnv
        str:  The QueryKey
    """
    Entrez.email = email

    handle = Entrez.esearch(db='pubmed', 
        term=entrez_query, 
        retmax=max_results, 
        retmode="xml", 
        usehistory='y'
        )

    record = Entrez.read(handle)
    handle.close()

    ids = record["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]

    return ids, webenv, query_key

def search_pubmed(entrez_query: str, email: str, max_results: int = 50):
    """
    Generate the list of search Entrez to get the list of abstrac IDs, fetch data, and
    generate documents.

    Args:
        email (str): e-mail addressed used with the Entrez query
        entrez_query (str): The query sent to Entrez
        max_results (int):  The maxiumum number of abstracts to return from Entrez

    Return:
        list: The collection of documents
    """
    ids, webenv, query_key = get_ids(email, entrez_query, max_results)

    if len(ids) == 0:
        entrez_query = retry_entrez_query(entrez_query)
        ids, webenv, query_key = get_ids(email, entrez_query, max_results)
    if len(ids) == 0:
        raise ValueError("Entrez query returned no results after retry.")

    # Fetch XML data
    fetch_handle = Entrez.efetch(db="pubmed", 
        id=",".join(ids), 
        rettype="xml",
        retmode="xml",
        webenv=webenv,
        query_key=query_key)
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

        doc = Document(page_content=f"{abstract}",
                       metadata={"title": title, "date": date, "PMID": str(pmid)}
                       )
        documents.append(doc)
        
    return documents
