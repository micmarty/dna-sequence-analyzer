import click
from SequenceAnalyzer import SequencesAnalyzer


@click.command()
@click.argument('sequence_a')
@click.argument('sequence_b')
@click.option('-S', '--summary', is_flag=True)
@click.option('-s', '--similarity', is_flag=True)
@click.option('-e', '--edit-distance', is_flag=True)
@click.option('-a', '--alignment', type=click.Choice(['global', 'local']))
def main(summary, similarity, edit_distance, sequence_a, sequence_b, alignment):
    analyzer = SequencesAnalyzer(sequence_a, sequence_b)

    if summary:
        analyzer.edit_distance()
        analyzer.similarity()
        analyzer.local_alignment()
        analyzer.global_alignment()
    if similarity:
        analyzer.similarity()
    if edit_distance:
        analyzer.edit_distance()

    if alignment == 'local':
        analyzer.local_alignment()
    elif alignment == 'global':
        analyzer.global_alignment()


if __name__ == '__main__':
    main()
