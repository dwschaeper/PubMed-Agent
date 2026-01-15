from langchain.chat_models import init_chat_model
from langchain_core.prompts import PromptTemplate
from langchain_core.output_parsers import StrOutputParser


def single_summary(chain, document):
    """
    Generate a summary of a single abstract.

    Args:
        chain (Chain):  The chain used to summarize an abstract.
        document (str): The document data passed as a str from previous tool

    Return
        dict: A dictionary that has the summary as the key, and metadata as the value
    """
    result = chain.invoke({'text': str(document)})
    return result.strip()


def summarize_abstracts(documents):
    """
    Summarize of LangChain Document objects containing PubMed abstracts.
    
    Args:
        documents (list): List of documents data from the previous tool containing abstracts
    
    Return:
        str: Complete summary of all contained abstracts 
    """
    # define model
    llm = init_chat_model("llama-3.3-70b-versatile", model_provider="groq")

    # define prompts
    summarize_prompt = PromptTemplate(input_variables=["text"],
                                      template="Summarize this abstract with specific details. Include PMID in summary.\n{text}")

    # make chains
    summarize_chain = summarize_prompt | llm | StrOutputParser()

    # make summary
    summaries = [single_summary(summarize_chain, document) for document in documents]
    result = '\n'.join(summaries)

    return result
