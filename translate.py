import click
from Translator import Translator


@click.command()
@click.argument('sequence', required=False)
@click.option('-i', '--input-file', type=click.Path(exists=True, dir_okay=False, readable=True), 
    help='Path text file containing long nucleotide sequence')
def main(sequence, input_file=None):
    if input_file:
        with open(input_file) as f:
            lines = f.readlines()
            print_protein = lambda seq: print(Translator(seq).to_protein)
            [print_protein(line.strip()) for line in lines]
    elif sequence:
        print(Translator(sequence).to_protein)
    else:
        click.echo('Missing argument. Please run --help.')


if __name__ == '__main__':
    main()
