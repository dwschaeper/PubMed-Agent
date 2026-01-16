from pydantic import BaseModel
from typing import List, Optional

class AgentState(BaseModel):
    user_query: str
    email: str
    entrez_query: Optional[str] = None
    query_ready: bool = False
    documents: Optional[List[dict]] = None
    num_documents: int = 0
    summary: Optional[str] = None
    errors: List[str] = []
    done: bool = False
    thoughts: List[str] = []
