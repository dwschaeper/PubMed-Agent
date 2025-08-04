from src.pubed_search import search_pubmed

import argparse


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--email', required=True, type=str)
    parser.add_argument('--query', required=True, type=str)

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse()

    data = search_pubmed(args.email, query=args.query)

    print(data)
