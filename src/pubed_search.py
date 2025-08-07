from .input_processing import retry_query

import requests
import xml.etree.ElementTree as ET
from langchain.docstore.document import Document

def get_ids(email, query, max_results=50):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "email": email
    }

    response = requests.get(base_url, params=params)
    response.raise_for_status()

    data = response.json()
    id_list = data["esearchresult"]["idlist"]

    return id_list


def get_data(ids, email):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    id_string = ",".join(ids)

    params = {
        "db": "pubmed",
        "id": id_string,
        "retmode": "xml",
        "rettype": "abstract",
        "email": email
    }

    response = requests.get(base_url, params=params)
    response.raise_for_status()

    return response.text

def parse_data(data):
    root = ET.fromstring(data)
    articles = []

    for article in root.findall(".//PubmedArticle"):
        pmid_elem = article.find(".//PMID")
        title_elem = article.find(".//ArticleTitle")
        abstract_elems = article.findall(".//Abstract/AbstractText")
        date_elem = article.find(".//Article/ArticleDate")
        journal_date = article.find(".//JournalIssue/PubDate")

        pmid = pmid_elem.text if pmid_elem is not None else "unknown"
        title = title_elem.text if title_elem is not None else "No Title"
        abstract = " ".join([a.text or "" for a in abstract_elems]) if abstract_elems else "No Abstract"

        if date_elem is not None:
            year = date_elem.findtext("Year", default="0000")
            month = date_elem.findtext("Month", default="01")
            day = date_elem.findtext("Day", default="01")
        elif journal_date is not None:
            year = journal_date.findtext("Year", default="0000")
            month = journal_date.findtext("Month", default="01")
            day = journal_date.findtext("Day", default="01")
        else:
            year, month, day = "0000", "01", "01"

        date = f"{year}-{month.zfill(2)}-{day.zfill(2)}"

        articles.append({
            "PMID": pmid,
            "title": title,
            "abstract": abstract,
            "date": date
        })

    return articles

def search_pubmed(email, query, max_results=50):
    ids = get_ids(email, query, max_results)

    if len(ids) == 0:
        query = retry_query(query)
        ids, webenv, query_key = get_ids(email, query, max_results)
    if len(ids) == 0:
        print("Query was not successful after retry.")
        exit(1)

    xml_data = get_data(ids, email)
    articles = parse_data(xml_data)

    documents = []
    for article in articles:
        doc = Document(page_content=f"Title: {article['title']}\nAbstract: {article['abstract']}\nPMID: {article['PMID']}",
                       metadata={"title": article['title'], "date": article['date'], "PMID": article['PMID']})
        documents.append(doc)

    return documents
