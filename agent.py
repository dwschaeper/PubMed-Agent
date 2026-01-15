import argparse
from dotenv import load_dotenv
from langchain.chat_models import init_chat_model
from langchain.agents import create_agent
from pydantic import BaseModel, Field

from src.tools import make_pubmed_query_tool, search_pubmed_tool, summarize_abstracts_tool


class SummaryResponse(BaseModel):
    """Summary of abstracts"""
    summary: str = Field(description="The final summary")



def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--email', required=True, type=str)
    parser.add_argument('--query', required=True, type=str)

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    SYSTEM_PROMPT = """You are an expect scientific researcher talented at summarizing research articles. 
                       You have access to two tools:
                           - make_pubmed_query_tool: use this to convert the query privided into a high quality PubMed query
                           - search_pubmed_tool: use this to take the created high quality PubMed query to search Entrez with Biopython. Only use the supplied email address
                           - summarize_abstracts_tool: use this to summarize the abstracts that are output from search_pubmed_tool
                       The tools should be called in this order to complete the task: make_pubmed_query_tool, search_pubmed_tool, summarize_abstracts_tool. Do not invent parameters, only use those supplied.
                    """
    args = parse()
    load_dotenv()

    llm = init_chat_model("llama-3.3-70b-versatile", model_provider="groq")
    tools = [make_pubmed_query_tool, search_pubmed_tool, summarize_abstracts_tool]

    agent = create_agent(model=llm, tools=tools, response_format=SummaryResponse)
    response = agent.invoke({"messages":[{'role': 'user', 'content': f'Make a high quality pubmed query from this query: {args.query}. Use this email "{args.email}" to query PubMed when it is time to call the function.'}]})


    print('\n--- FINAL SUMMARY ---\n')
    print(response['structured_response'].summary)

