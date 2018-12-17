import click
from Translator import Translator


@click.command()
@click.argument('sequence')
def main(sequence):
    print(Translator(sequence).to_protein)


if __name__ == '__main__':
    main()
