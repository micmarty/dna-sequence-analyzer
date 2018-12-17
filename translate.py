import click
from Translator import Translator


@click.command()
@click.argument('sequence')
@click.option('-r', '--rna', is_flag=True)
def main(sequence, rna):
    if rna:
        print(Translator.rna_sequence_to_protein(sequence))


if __name__ == '__main__':
    main()
