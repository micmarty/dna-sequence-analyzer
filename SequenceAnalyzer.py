import numpy as np
from typing import Tuple, Dict, Any
'''
Author: Michal Martyniak (github: @micmarty)

Helpful resources: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/local.pdf
Note: Smith Waterman and Needleman Wunsch algorithms are very, very similar but it's meant to be separated
'''


class SequencesAnalyzer:
    # similarity_matrix = {
    #     'A': {'A': 1, 'G': 2, 'C': 2, 'T': 2, 'U': 2, '-': 0},
    #     'G': {'A': -1, 'G': 1, 'C': 2, 'T': 2, 'U': 2, '-': 0},
    #     'C': {'A': -1, 'G': 2, 'C': 1, 'T': 2, 'U': 2, '-': 0},
    #     'T': {'A': -1, 'G': 2, 'C': 2, 'T': 1, 'U': 2, '-': 0},
    #     'U': {'A': -1, 'G': -1, 'C': 2, 'T': 2, 'U': 1, '-': 0},
    #     '-': {'A': 0, 'G': 0, 'C': 0, 'T': 0, 'U': 0, '-': 1}
    # }
    # replacement_cost = {
    #     'A': {'A': 0, 'G': 2, 'C': 2, 'T': 2, 'U': 2, '-': 1},
    #     'G': {'A': 2, 'G': 0, 'C': 2, 'T': 2, 'U': 2, '-': 1},
    #     'C': {'A': 2, 'G': 2, 'C': 0, 'T': 2, 'U': 2, '-': 1},
    #     'T': {'A': 2, 'G': 2, 'C': 2, 'T': 0, 'U': 2, '-': 1},
    #     'U': {'A': 2, 'G': 2, 'C': 2, 'T': 2, 'U': 0, '-': 1},
    #     '-': {'A': 1, 'G': 1, 'C': 1, 'T': 1, 'U': 1, '-': 0}
    # }

    traceback_symbols = {
        0: '↖',
        1: '↑',
        2: '←',
        3: '•'
    }

    def _traceback(self, result_matrix, traceback_matrix, start_pos: Tuple[int, int]) -> Tuple[str, str]:
        '''Use both matrices to replay the optimal route'''
        seq_a_aligned = ''
        seq_b_aligned = ''

        # 1. Select starting point
        row, col = start_pos

        # 2. Terminate when 0 is reached (end of path)
        while result_matrix[row, col] != 0:
            symbol = traceback_matrix[row, col]

            # 3. Use arrows to navigate and collect letters (in reversed order)
            # Shift indexes by one (matrix has additional row and column)
            if symbol == '↖':
                seq_a_aligned += self.seq_a[row - 1]
                seq_b_aligned += self.seq_b[col - 1]
                row -= 1
                col -= 1
            elif symbol == '↑':
                seq_a_aligned += self.seq_a[row - 1]
                seq_b_aligned += '-'
                row -= 1
            elif symbol == '←':
                seq_a_aligned += '-'
                seq_b_aligned += self.seq_b[col - 1]
                col -= 1
        # Reverse strings (traceback goes from bottom-right to top-left)
        return seq_a_aligned[::-1], seq_b_aligned[::-1]

    # TODO Refactor score and edit_cost into separate class, allow for loading from file
    def score(self, a, b):
        assert isinstance(a, str) and isinstance(b, str)
        assert len(a) == 1 and len(b) == 1

        match, mismatch, gap = 1, -1, -1
        if a == b:
            return match
        else:
            if '-' in [a, b]:
                return gap
            return mismatch

    def edit_cost(self, a, b):
        assert isinstance(a, str) and isinstance(b, str)
        assert len(a) == 1 and len(b) == 1

        match, mismatch, gap = 0, 1, 1
        if a == b:
            return match
        else:
            if '-' in [a, b]:
                return gap
            return mismatch

    def __init__(self, seq_a: str, seq_b: str) -> None:
        self.seq_a = seq_a
        self.seq_b = seq_b

    def needleman_wunsch_algorithm(self, minimize: bool) -> Dict[str, Any]:
        '''
        Dynamic programming technique
        Reference: [l3a.pdf, slide #5] + [https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm]
        Algorithm: Needleman-Wunch
        Time complexity: O(nm)
        Space complexity: O(nm)

        `minimize` is a flag which needs to be enabled when calculating edit distance
        '''
        # 1. Prepare dimensions (required additional 1 column and 1 row)
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1

        # 2. Initialize matrices
        # Use grid/matrix as graph-like acyclic digraph (array cells are vertices)
        H = np.zeros(shape=(rows, cols), dtype=int)
        traceback = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))

        # Put sequences' letters into first row and first column (better visualization)
        traceback[0, 1:] = np.array(list(self.seq_b), dtype=str)
        traceback[1:, 0] = np.array(list(self.seq_a), dtype=str)

        # 3. Top row and leftmost column, like: 0, 1, 2, 3, etc.
        H[0, :] = np.arange(start=0, stop=cols)
        H[:, 0] = np.arange(start=0, stop=rows)

        for row in range(1, rows):
            for col in range(1, cols):
                # Current pair of letters from sequence A and B
                a = self.seq_a[row - 1]
                b = self.seq_b[col - 1]

                if minimize:
                    # Required for edit cost calculation
                    score_func = self.edit_cost
                else:
                    # Required for similarity calculation
                    score_func = self.score

                leave_or_replace_letter = H[row-1, col-1] + score_func(a, b)
                delete_indel = H[row-1, col] + score_func('-', b)
                insert_indel = H[row, col-1] + score_func(a, '-')

                scores = [leave_or_replace_letter, delete_indel, insert_indel]

                if minimize:
                    best_action = np.argmin(scores)
                else:
                    best_action = np.argmax(scores)

                H[row, col] = scores[best_action]
                traceback[row, col] = self.traceback_symbols[best_action]

        return {
            'result_matrix': H,
            'traceback_matrix': traceback,
            'score': H[-1, -1],
            'score_pos': (rows - 1, cols - 1)
        }

    def smith_waterman_algorithm(self) -> Dict[str, Any]:
        # TODO Add description similar to needleman-wunsch
        # 1. Prepare dimensions (required additional 1 column and 1 row)
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1

        # 2. Initialize matrices
        # Use grid/matrix as graph-like acyclic digraph (array cells are vertices)
        H = np.zeros(shape=(rows, cols), dtype=int)
        traceback = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))

        # Put sequences' letters into first row and first column (better visualization)
        traceback[0, 1:] = np.array(list(self.seq_b), dtype=str)
        traceback[1:, 0] = np.array(list(self.seq_a), dtype=str)

        # 3. Top row and leftmost colum are already 0
        for row in range(1, rows):
            for col in range(1, cols):
                # Alias: current pair of letters
                a = self.seq_a[row - 1]
                b = self.seq_b[col - 1]

                leave_or_replace_letter = H[row-1, col-1] + self.score(a, b)
                delete_indel = H[row-1, col] + self.score('-', b)
                insert_indel = H[row, col-1] + self.score(a, '-')

                # Zero is required - ignore negative numbers
                scores = [leave_or_replace_letter,
                          delete_indel, insert_indel, 0]
                best_action = np.argmax(scores)

                H[row, col] = scores[best_action]
                traceback[row, col] = self.traceback_symbols[best_action]

        return {
            'result_matrix': H,
            'traceback_matrix': traceback,
            'score': H.max(),
            'score_pos': np.unravel_index(np.argmax(H, axis=None), H.shape)
        }

    def local_alignment(self) -> Tuple[str, str]:
        result = self.smith_waterman_algorithm()
        alignment_a, alignment_b = self._traceback(
            result_matrix=result['result_matrix'],
            traceback_matrix=result['traceback_matrix'],
            start_pos=result['score_pos'])

        print(result['result_matrix'])
        print(result['traceback_matrix'])
        print('[Local Alignment] Score={}'.format(result['score']))
        print(alignment_a)
        print(alignment_b)
        return alignment_a, alignment_b

    def similarity(self) -> int:
        result = self.needleman_wunsch_algorithm(minimize=False)

        print(result['result_matrix'])
        print(result['traceback_matrix'])
        print('[Similarity] Score={}'.format(result['score']))
        return result['score']

    def edit_distance(self) -> int:
        result = self.needleman_wunsch_algorithm(minimize=True)

        print(result['result_matrix'])
        print(result['traceback_matrix'])
        print('[Edit distance] Cost={}'.format(result['score']))
        return result['score']
