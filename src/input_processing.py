from langchain.chat_models import init_chat_model
from langchain_core.prompts import PromptTemplate
from langchain_core.output_parsers import StrOutputParser


def make_entrez_query(user_query: str):
    """
    Convert the input query into a high-quality Entrez/PubMed query.

    Args:
        user_query (str): The user provided query

    Return:
        str: Optimized query
    """
    # define model
    llm = init_chat_model("llama-3.1-8b-instant", model_provider="groq")

    prompt = PromptTemplate(input_variables=["user_query"], template=(
        "Take the query and convert it into a high quality search for PubMed, no code. Return ONLY the query and nothing else.\n{user_query}"))

    chain = prompt | llm | StrOutputParser()
    result = chain.invoke({"user_query": user_query})

    return result.strip()


def remake_entrez_query(bad_query: str):
    """
    Revise an Entrez/PubMed query that returned no results.

    Args:
        bad_query (str): The previously used query that did not find abstracts
    Return:
        str: Revised query
    """
    llm = init_chat_model("llama-3.1-8b-instant", model_provider="groq")

    prompt = PromptTemplate(input_variables=["bad_query"], template=("""This query '{bad_query}' for PubMed did not successfully return any abstracts. Follow this strategy to make it work
        - remove overly strict fields
        - reduce excessive AND constraints
        - prefer OR when reasonable
        - keep original intent"""))

    chain = prompt | llm | StrOutputParser()
    result = chain.invoke({"bad_query": bad_query})

    return result.strip()
