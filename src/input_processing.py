from langchain.chat_models import init_chat_model
from langchain.prompts import PromptTemplate
from langchain_core.output_parsers import StrOutputParser


def make_query(query):
    # define model
    llm = init_chat_model("llama3-8b-8192", model_provider="groq")

    prompt = PromptTemplate(input_variables=["text"],
                            template="Take the query and convert it into a high quality Entrez efetch search"
                                     " for PubMed. Return ONLY the query and nothing else. Return ONLY the "
                                     "query and none of your extra text at the beginning or end of the "
                                     "output.\n{query}")

    chain = prompt | llm | StrOutputParser()
    result = chain.invoke({"query": query})

    return result.strip()


def retry_query(query):
    llm = init_chat_model("llama3-8b-8192", model_provider="groq")

    prompt = PromptTemplate(input_variables=["query"], template=(
        "This query for Entrez efetch search did not successfully return any abstracts. Adjust it to"
        " have it work. Return ONLY the query and none of your extra text at the beginning or end of the output."
        "\n{query}"))

    chain = prompt | llm | StrOutputParser()
    result = chain.invoke({"query": query})

    return result.strip()
