import click

# NOTE: To use environment variables as click.option default values, .env file should be loaded before importing modules
from dotenv import load_dotenv

load_dotenv()

from rxndata2.presentation.draw_reactions import draw_reactions


@click.group()
def visualize():
    pass


visualize.add_command(draw_reactions)  # type: ignore

if __name__ == "__main__":
    visualize()
