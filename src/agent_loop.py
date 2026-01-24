from copy import deepcopy

from .agent_state import AgentState
from .tools import make_entrez_query_tool, remake_entrez_query_tool, search_pubmed_tool, summarize_abstracts_tool
from dotenv import load_dotenv
from langchain.chat_models import init_chat_model
from langchain_core.prompts import PromptTemplate
from langchain_core.output_parsers import JsonOutputParser


load_dotenv()
llm = init_chat_model("llama-3.3-70b-versatile", model_provider="groq")

PLANNER_PROMPT = """
    You are an autonomous biomedical research agent with the goal of making detailed summaries of research articles.

    Current state:
    {state_json}

    Available tools:
    - make_entrez_query(query)
    - remake_entrez_query(bad_query)
    - search_pubmed(email, query)
    - summarize_abstracts(documents)
    - request_human_review_query(query)
    - request_human_review_documents(documents)

    Use this logic to make a decision on what to do next:
    - Utilize the action history, looping is an unacceptable behavior. This will be especially helpful if a query fails to know what step was previously performed.
    - If entrez_query in current state is None, generate a query with make_entrez_query.
    - A human must review the query with request_human_review.
    - Take entrez_query and use it to search PubMed for abstracts with search_pubmed. pubmed_search_failed state will reflect the result.
    - If PubMed search failed, retry making the query with remake_entrez_query or make_entrez_query.
    - A human must review the query with request_human_review_documents.
    - After the documents are approved, summarize. Check action history, if documents have not been sent for approval, do not summarize.
    - If summary exists, mark done=True.
    - Provide a brief justification referencing the state fields used.


    Return JSON only in this format:
    {{
       "action": "name of the tool to call next",
       "args": {{}},
       "reasoning": "Explain why this tool was chosen next."
    }}
        
    ONLY OUTPUT THE JSON AS SPECIFIED ABOVE AS RAW TEXT, NOT CONTAINED IN MARKDOWN OR CODE
    """

def planner(state: AgentState):
    """
    Call LLM to decide next action.
	
    Args:
        state (AgentState): The current state of the agent
    """
    prompt = PromptTemplate(
        input_variables=["state_json"],
        template=PLANNER_PROMPT
    )
    chain = prompt | llm | JsonOutputParser()
    decision = chain.invoke({"state_json": state.model_dump_json()})
    return decision


def state_view(state: AgentState):
    """
    Print the current state of the agent.

    Args:
        state (AgentState): The current state of the agent
    """
    print('[CURRENT STATE]')
    print('\tUser query:', state.user_query)
    print('\tEmail:', state.email)
    print('\tEntrez query:', state.entrez_query)
    print('\tQuery ready:', state.query_ready)
    print('\tQuery approved by user:', state.query_approved)
    print('\tNumber of documents:', state.num_documents)
    print('\tPubMed searched:', state.pubmed_searched)
    print('\tPubMed search failed:', state.pubmed_search_failed)
    print('\tDocuments approved by user:', state.documents_approved)
    print('\tErrors:', state.errors)
    print('\tDone:', state.done)


def run_agent(user_query: str, email: str):
    """
    Begin the agent.

    Args:
        user_query (str): The query the agent is asked to search
        email (str):      The email to use when searching PubMed

    Return:
        AgentState: The current state of the agent
    """
    state = AgentState(user_query=user_query, email=email)
    previous_state = deepcopy(state)

    while not state.done:
        decision = planner(state)
        reasoning = decision.get("reasoning", "")

        action = decision["action"].strip()
        args = decision.get("args", {})

        # --- status update ---
        print("\n[STATUS] Agent reasoning:")
        print("-", reasoning)
        print('Taking action: ', action)
        state_view(state)
        

        try:
            state.ordered_action_history.append(action)
            if action == "make_entrez_query":
                state.entrez_query = make_entrez_query_tool.func(**args)
                state.query_ready = True
            elif action == "remake_entrez_query":
                state.entrez_query = remake_entrez_query_tool.func(**args)
            elif action == "search_pubmed":
                state.documents = search_pubmed_tool.func(**args)
                state.pubmed_searched = True
                if len(state.documents) == 0:
                    state.pubmed_search_failed = True
                else:
                    state.pubmed_search_failed = False
                    state.num_documents = len(state.documents)
            elif action == "summarize_abstracts":
                state.summary = summarize_abstracts_tool.func(state.documents)
                state.done = True
            elif action == 'request_human_review_query':
                print("\n[HUMAN REVIEW REQUIRED]")
                print("Proposed query:\n", args["query"])

                decision = input("Approve this query? (y / edit / n): ")

                if decision.lower() == 'y':
                    state.query_approved = True
                elif decision.lower() == 'edit':
                    edited = input('Enter revised query:\n')
                    state.entrez_query = edited
                    state.query_approved = True
                else:
                    state.entrez_query = None
                    state.human_demand = 'make_entrez_query'
                    state.query_ready = False
            elif action == "request_human_review_documents":
                print("\n[HUMAN REVIEW REQUIRED]")
                print("Proposed documents:\n")
                for document in state.documents:
                    print('\t', document['title'])

                decision = input("Approve these documents? (y / edit / n): ")

                if decision.lower() == 'y':
                    state.query_ready = True
                    state.documents_approved = True
                else:
                    state.entrez_query = None
                    state.human_demand = 'make_entrez_query'

        except Exception as e:
            state.errors.append(str(e))

        # state updating
        if state.human_demand == action:
            state.human_demand = None
        if not previous_state.query_approved and (action == 'make_entrez_query' or action == 'remake_entrez_query'):
            state.query_approved = None


    return state