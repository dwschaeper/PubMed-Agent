from src.input_processing import make_query
from src.pubed_search import search_pubmed
from src.summarizer import summarize_abstracts

import argparse
from dotenv import load_dotenv


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--email', required=True, type=str)
    parser.add_argument('--query', required=True, type=str)

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse()
    load_dotenv()

    query = make_query(args.query)

    documents = search_pubmed(args.email, query=query, max_results=50)

    summary = summarize_abstracts(documents)

    print(summary)