from .input_processing import make_entrez_query, remake_entrez_query
from .pubmed_search import search_pubmed
from .summarizer import summarize_abstracts

from langchain.tools import tool

@tool("make_entrez_query")
def make_entrez_query_tool(query: str) -> str:
    """
    Convert a natrual language question into a high-quality Entrez 
    efetch search string.
    """
    return make_entrez_query(query)

@tool("remake_entrez_query')
def remake_entrez_query_tool(bad_query: str) -> str:
    """
    Loosen, simplify, or reowkr a PubMed/Entrez query that returned 0 results. Remove overly strict filters, rare fields,
    and excessive AND clauses.
    """
    return remake_entrez_query(bad_query)

@tool("search_pubmed")
def search_pubmed_tool(email: str, query: str) -> list:
    """
    Search PubMed using BioPython Entrez and return a list of documents.
    Each entry in the returned list called document_details is a dictionary that has keys
        - abstract: The value is the abstract of the article
        - PMID: The value is the PMID of the article
        - publish_date: The value is the date the article was published, if it is known
    """
    documents = search_pubmed(query, email=email, max_results=5)
    document_details = [{'abstract': d.page_content,
            'PMID': d.metadata.get("PMID"),
            'title': d.metadata.get("title"),
            'publish_date': d.metadata.get("date")} for d in documents]
    return document_details
    

@tool("summarize_abstracts")
def summarize_abstracts_tool(document_details) -> str:
    """
    Summarize a list of PubMed article abstracts into a single combined scientific summary.
    This takes the document_details output from search_pubmed_tool
    """
    summary = summarize_abstracts(document_details)
    return summary
