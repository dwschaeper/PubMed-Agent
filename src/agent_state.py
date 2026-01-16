from pydantic import BaseModel
from typing import List, Optional

class AgentState(BaseModel):
    user_query: str
    email: str
    entrez_query: Optional[str] = None
    query_ready: bool = False
    documents: Optional[List[dict]] = None
    pubmed_searched: bool = False
    pubmed_search_failed: bool = None
    summary: Optional[str] = None
    ordered_action_history: List[list] = []
    errors: List[str] = []
    done: bool = False
