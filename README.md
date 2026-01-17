# PubMed Agent

PubMed Agent is an agentic biomedical research assistant that autonomously turns a natural
language query into a PubMed search, retrieves abstracts, and produces a structured
scientific summary using LLM reasoning and tool orchestration.

The LLMs served via [Groq](https://console.groq.com/home) so a Groq API key is needed.

The agent operates in a loop where an LLM planner inspects its current state and decides which
tool to call next. The initial version provides the agent with four tools:

1. make_entrez_query - generate a high-quality query for Entrez
2. search_pubmed - Use the query to search PubMed
3. remake_entrez_query - re-create the query in the case it was not successful
4. summarize_abstracts - Generate a summary of the abstracts

The agent continues until it successfully produces a summary or determines
it cannot progress further.

## Setup

1. Clone
```
git clone https://github.com/dwschaeper/PubMed-Agent
cd PubMed-Agent
```
2. Install requirements
```
pip install -r requirements.txt
```
3. Provide Groq API Key
```commandline
echo 'GROQ_API_KEY=[your_key_here]' > .env
```
4. Run
```
python agent.py --email [email] --query [query]
```

## Output
The program prints:

 - Step-by-step agent reasoning and current state

 - Tools selected at each step

 - Final summarized report of PubMed abstracts

## Future Improvements
 - Add relevance-scoring tool

 - Use tool feedback to tune queries

 - Optional human-in-the-loop