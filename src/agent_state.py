from pydantic import BaseModel
from typing import List, Optional

class AgentState(BaseModel):
    user_query: str
    email: str
    entrez_query: Optional[str] = None
    query_ready: bool = False
    documents: Optional[List[dict]] = None
    pubmed_searched: bool = False
    summary: Optional[str] = None
    errors: List[str] = []
    done: bool = False
    thoughts: List[str] = []
