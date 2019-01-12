import click
from SequenceAnalyzer import SequencesAnalyzer
from HirschbergAlgorithm import HirschbergAlgorithm
from NeedlemanWunschAlgorithm import NeedlemanWunschAlgorithm
from ScoringSystem import ScoringSystem

@click.command()
@click.argument('sequence_a')
@click.argument('sequence_b')
@click.option('-S', '--summary', is_flag=True)
@click.option('-s', '--similarity', is_flag=True)
@click.option('-e', '--edit-distance', is_flag=True)
@click.option('-a', '--alignment', type=click.Choice(['global', 'local']))

@click.option('--gap-start', default=0, help='When gap starts, that score will be applied (besides `gap` penalty).')
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
        print('Needleman-Wunsch (with start gap penalty)')
        if load_csv:
            scoring_sys = ScoringSystem(gap_start=gap_start)
            scoring_sys.load_csv('scores.csv')
            # Show what's inside the files
            print('[Custom scoring system]\n', scoring_sys)
        else:
            scoring_sys = ScoringSystem(match=match, mismatch=mismatch, gap_start=gap_start, gap=gap)
        NeedlemanWunschAlgorithm(scoring_sys).align(sequence_a, sequence_b, start_gap_penalty_mode=True)

        print('Hirschberg (without start gap penalty)')
        HirschbergAlgorithm(scoring_sys).align(sequence_a, sequence_b)

if __name__ == '__main__':
    main()
