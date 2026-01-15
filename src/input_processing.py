from langchain.chat_models import init_chat_model
from langchain_core.prompts import PromptTemplate
from langchain_core.output_parsers import StrOutputParser

def make_entrez_query(user_query: str):
    """
    Convert the input query into a high-quality Entrez/PubMed query.

    Args:
        query (str): The user provided query
    Return:
        str: Optimized query
    """
    # define model
    llm = init_chat_model("llama-3.3-70b-versatile", model_provider="groq")

    prompt = PromptTemplate(input_variables=["user_query"], template=("Take the query and convert it into a high quality Entrez efetch search from biopython for PubMed. Return ONLY the query and nothing else. Return ONLY the query and none of your extra text as the beginning or end of the output.\n{user_query}"))

    chain = prompt | llm | StrOutputParser()
    result = chain.invoke({"user_query": user_query})

    return result.strip()


def retry_entrez_query(bad_query: str):
    """
    Revise an Entrez/PubMed query that failed to return any abstracts.

    Args:
        bad_query (str): The preiously used query that did not find abstracts
    Return:
        str: Revised query
    """
    llm = init_chat_model("llama-3.1-8b-instant", model_provider="groq")

    prompt = PromptTemplate(input_variables=["bad_query"], template=("This query for Entrez efetch search from biopython did not successfully return any abstracts. Adjust it to have it work.\n{bad_query}"))

    chain = prompt | llm | StrOutputParser()
    result = chain.invoke({"bad_query": bad_query})

    return result.strip()

