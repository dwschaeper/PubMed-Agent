import argparse
from dotenv import load_dotenv
from pydantic import BaseModel, Field

from src.agent_loop import run_agent


class SummaryResponse(BaseModel):
    """Summary of abstracts"""
    summary: str = Field(description="The final summary")


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--email', required=True, type=str)
    parser.add_argument('--query', required=True, type=str)

    arguments = parser.parse_args()

    return arguments


if __name__ == '__main__':
    args = parse()
    load_dotenv()

    final_state = run_agent(user_query=args.query, email=args.email)

    print("\n--- FINAL SUMMARY ---\n")
    print(final_state.summary)
