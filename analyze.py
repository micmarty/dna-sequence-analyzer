import click
from SequenceAnalyzer import SequencesAnalyzer
from HirschbergAlgorithm import HirschbergAlgorithm
from ScoringSystem import ScoringSystem

@click.command()
@click.argument('sequence_a')
@click.argument('sequence_b')
@click.option('-S', '--summary', is_flag=True)
@click.option('-s', '--similarity', is_flag=True)
@click.option('-e', '--edit-distance', is_flag=True)
@click.option('-a', '--alignment', type=click.Choice(['global', 'local']))

@click.option('--gap-start', default=-1, help='When gap starts, that score will be applied (besides `gap` penalty).')
@click.option('--gap', default=-1, help='Penalty for a single gap.')
@click.option('--match', default=1, help='Score for identical letters in both seqences.')
@click.option('--mismatch', default=-1, help='Score for mismatched letters in both seqences.')
@click.option('--load-csv', is_flag=True, help='Load scores.csv and edit_cost.csv')
@click.option('--verbose', is_flag=True, help='Verbose mode.')
def main(load_csv, summary, similarity, edit_distance, sequence_a, sequence_b, alignment, gap_start, gap, match, mismatch, verbose):
    if verbose:
        print(
            f'Gap start penalty: {gap_start}\n'
            f'Gap penalty: {gap}\n'
            f'Match: {match}\n'
            f'Mismatch: {mismatch}\n'
        )
    analyzer = SequencesAnalyzer(sequence_a, sequence_b, load_csv=load_csv)

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
        print('--------------------------')
        #analyzer.hirschberg_algorithm(X=analyzer.seq_a, Y=analyzer.seq_b)
        scoring_sys = ScoringSystem(match=2, mismatch=-1, gap=-2)
        HirschbergAlgorithm(scoring_sys).align(sequence_a, sequence_b)


if __name__ == '__main__':
    main()
