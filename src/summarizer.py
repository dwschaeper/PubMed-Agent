from langchain.chat_models import init_chat_model
from langchain.prompts import PromptTemplate
from langchain_core.output_parsers import StrOutputParser


def single_summary(chain, document):
    result = chain.invoke({'text': document.page_content})
    return {'text': result.strip(), 'metadata': document.metadata}


def combine_summaries(chain, summaries):
    while len(summaries) > 1:
        new_summaries = []
        for i in range(0, len(summaries), 2):
            if i + 1 < len(summaries):
                combined_text = summaries[i]['text'] + "\n\n" + summaries[i + 1]['text']
                result = chain.invoke({'text': combined_text})
                new_summaries.append({'text': result.strip()})
            else:
                new_summaries.append(summaries[i])

        summaries = new_summaries
    return summaries[0]['text']


def summarize_abstracts(documents):
    # define model
    llm = init_chat_model("llama3-8b-8192", model_provider="groq")

    # define prompts
    summarize_prompt = PromptTemplate(input_variables=["text"], template=("You're a biomedical research assistant. "
                                                                          "Summarize this abstract with specific "
                                                                          "details. You MUST include the PMID in summary."
                                                                          " ONLY print summary, don't say here it is.\n{text}"))

    combine_prompt = PromptTemplate(input_variables=["text"],
                                    template=("You are a scientific editor who has received the "
                                              "following summaries that have multiple articles "
                                              "in them:\n\n{text}\n\nCombine all articles into"
                                              " a single summary. Group the related papers "
                                              "together. You MUST include PMID for EACH paper "
                                              "in summary. If combining similar summaries, "
                                              "use all PMIDs associated with them. ONLY print"
                                              " summary, don't say here it is."))

    # make chains
    summarize_chain = summarize_prompt | llm | StrOutputParser()
    combine_chain = combine_prompt | llm | StrOutputParser()

    # make summary
    summaries = [single_summary(summarize_chain, document) for document in documents]
    result = combine_summaries(combine_chain, summaries)

    return result
