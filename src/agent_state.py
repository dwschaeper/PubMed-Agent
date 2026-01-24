from pydantic import BaseModel
from typing import List, Optional

class AgentState(BaseModel):
    user_query: str
    email: str
    entrez_query: str = None
    query_ready: bool = False
    query_approved: bool = False
    documents: List[dict] = None
    num_documents: int = None
    pubmed_searched: bool = False
    documents_approved: bool = False
    pubmed_search_failed: bool = None
    summary: str = None
    ordered_action_history: List[str] = []
    human_demand: str = None
    errors: List[str] = []
    done: bool = False
