import numpy as np
from typing import Tuple
'''
Author: Michal Martyniak (github: @micmarty)

Helpful resources:
https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/local.pdf
'''


class SequencesAnalyzer:
    similarity_matrix = {
        'A': {'A': 1, 'G': 2, 'C': 2, 'T': 2, 'U': 2, '-': 0},
        'G': {'A': -1, 'G': 1, 'C': 2, 'T': 2, 'U': 2, '-': 0},
        'C': {'A': -1, 'G': 2, 'C': 1, 'T': 2, 'U': 2, '-': 0},
        'T': {'A': -1, 'G': 2, 'C': 2, 'T': 1, 'U': 2, '-': 0},
        'U': {'A': -1, 'G': -1, 'C': 2, 'T': 2, 'U': 1, '-': 0},
        '-': {'A': 0, 'G': 0, 'C': 0, 'T': 0, 'U': 0, '-': 1}
    }
    replacement_cost = {
        'A': {'A': 0, 'G': 2, 'C': 2, 'T': 2, 'U': 2, '-': 1},
        'G': {'A': 2, 'G': 0, 'C': 2, 'T': 2, 'U': 2, '-': 1},
        'C': {'A': 2, 'G': 2, 'C': 0, 'T': 2, 'U': 2, '-': 1},
        'T': {'A': 2, 'G': 2, 'C': 2, 'T': 0, 'U': 2, '-': 1},
        'U': {'A': 2, 'G': 2, 'C': 2, 'T': 2, 'U': 0, '-': 1},
        '-': {'A': 1, 'G': 1, 'C': 1, 'T': 1, 'U': 1, '-': 0}
    }

    traceback_symbols = {
        0: '↖',
        1: '↑',
        2: '←',
        3: '•'
    }

    def traceback(self, result_matrix, traceback_matrix, start_pos: Tuple[int, int]) -> Tuple[str, str]:
        # 1.
        seq_a_aligned = ''
        seq_b_aligned = ''
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

    def score(self, a, b):
        #match, mismatch, gap = 3, -2, -1
        match, mismatch, gap = 2, -2, -1
        if a == b:
            return match
        elif '-' in [a, b]:
            return gap
        else:
            return mismatch

    def editing_cost(self, a, b):
        match, mismatch, gap = 0, 1, 1
        if a == b:
            return match
        elif '-' in [a, b]:
            return gap
        else:
            return mismatch

    def __init__(self, seq_a: str, seq_b: str) -> None:
        self.seq_a = seq_a
        self.seq_b = seq_b

    def needleman_wunsch_algorithm(self):
        '''
        Reference: https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm
        '''
        pass

    def local_alignment(self):
        result = self.smith_waterman_algorithm()
        alignment_a, alignment_b = self.traceback(
            result_matrix=result['result_matrix'],
            traceback_matrix=result['traceback_matrix'],
            start_pos=result['score_pos'])
        print('[Local Alignment] Score={}'.format(result['score']))
        print(alignment_a)
        print(alignment_b)
        
        

    def smith_waterman_algorithm(self) -> tuple:
        # 1. Prepare dimensions (required additional 1 column and 1 row)
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1

        # 2. Initialize matrices
        H = np.zeros(shape=(rows, cols), dtype=int)
        traceback = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))

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

    @property
    def similarity(self) -> int:
        '''
        Dynamic programming technique
        Reference: l3a.pdf, slide #5 
        Algorithm: Needleman-Wunch
        Time complexity: O(nm)
        Space complexity: O(nm)
        '''

        # Use grid/matrix as graph-like acyclic digraph (array cells are vertices)
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1
        D = np.zeros(shape=(rows, cols), dtype=int)
        backtrace = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))

        # 1. Initialize matrix
        # Top row and leftmost column, like: 0, 1, 2, 3, etc.
        D[0, :] = np.arange(start=0, stop=cols)
        D[:, 0] = np.arange(start=0, stop=rows)

        for row in range(1, rows):
            for col in range(1, cols):
                # Current pair of letters from sequence A and B
                a = self.seq_a[row - 1]
                b = self.seq_b[col - 1]

                # Determine what's the editing cost
                cost = self.replacement_cost[a][b]
                print(f'Current pair: "{a}" -> "{b}", cost={cost}')

                letter_replacement = D[row-1, col-1] + 1  # cost(a[i], b[i])
                delete_indel = D[row-1, col] + 1  # cost    # cost('-', b[i])
                insert_indel = D[row, col-1] + 1  # cost # cost(a[i], '-')

                operations_costs = [letter_replacement,
                                    delete_indel, insert_indel]
                min_operation = np.argmax(operations_costs)

                D[row, col] = operations_costs[min_operation]

                if min_operation == 0:
                    op = f'↖ ({cost})'
                elif min_operation == 1:
                    op = f'↑ ({cost})'
                elif min_operation == 2:
                    op = f'← ({cost})'
                backtrace[row, col] = op
                # print(
                #     f'[0] Letter replacement={letter_replacement} | ' \
                #     f'[1] Replace indel={delete_indel} | ' \
                #     f'[2] Insert indel={insert_indel} | ' \
                #     f'Operation: [{min_operation}]')
        # Bottom-right cell contains optimal distance
        min_distance = D[-1, -1]

        print('[Similarity matrix]\n', D)
        print('[Backtrace matrix] Follow the arrows (cost in brackets)\n', backtrace)
        print('Similarity =', min_distance)
        return min_distance

    @property
    def edit_distance(self) -> int:
        '''
        Dynamic programming technique
        Reference: l3a.pdf, slide #5 
        Algorithm: Needleman-Wunch
        Time complexity: O(nm)
        Space complexity: O(nm)
        '''

        # Use grid/matrix as graph-like acyclic digraph (array cells are vertices)
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1
        D = np.zeros(shape=(rows, cols), dtype=int)
        backtrace = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))

        # 1. Initialize matrix
        # Top row and leftmost column, like: 0, 1, 2, 3, etc.
        D[0, :] = np.arange(start=0, stop=cols)
        D[:, 0] = np.arange(start=0, stop=rows)

        for row in range(1, rows):
            for col in range(1, cols):
                # Current pair of letters from sequence A and B
                a = self.seq_a[row - 1]
                b = self.seq_b[col - 1]

                # Determine what's the editing cost
                #cost = self.replacement_cost[a][b]
                cost = self.editing_cost(a, b)
                print(f'Current pair: "{a}" -> "{b}", cost={cost}')

                letter_replacement = D[row-1, col-1] + cost  # cost(a[i], b[i])
                delete_indel = D[row-1, col] + cost    # cost('-', b[i])
                insert_indel = D[row, col-1] + cost  # cost(a[i], '-')

                operations_costs = [letter_replacement,
                                    delete_indel, insert_indel]
                min_operation = np.argmin(operations_costs)

                D[row, col] = operations_costs[min_operation]

                if min_operation == 0:
                    op = f'↖ ({cost})'
                elif min_operation == 1:
                    op = f'↑ ({cost})'
                elif min_operation == 2:
                    op = f'← ({cost})'
                backtrace[row, col] = op
                # print(
                #     f'[0] Letter replacement={letter_replacement} | ' \
                #     f'[1] Replace indel={delete_indel} | ' \
                #     f'[2] Insert indel={insert_indel} | ' \
                #     f'Operation: [{min_operation}]')
        # Bottom-right cell contains optimal distance
        min_distance = D[-1, -1]

        print('[Distance matrix]\n', D)
        print('[Backtrace matrix] Follow the arrows (cost in brackets)\n', backtrace)
        print('Min edit distance =', min_distance)
        return min_distance


#SequencesAnalyzer.edit_distance('AGTTT', 'AT-')
#SequencesAnalyzer('AGCACACA', 'ACACACTA').similarity




#SequencesAnalyzer('AACTTAC', 'TGAATTT').local_alignment()
SequencesAnalyzer('xyaxbacsll', 'pqraxabcstvq' ).local_alignment()
