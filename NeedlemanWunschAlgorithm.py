import numpy as np
from utils import traceback, traceback_symbols

class NeedlemanWunschAlgorithm:
    def __init__(self, scoring_sys) -> None:
        self.aligned_seq_a: str = ''
        self.aligned_seq_b: str = ''
        self.scoring_sys = scoring_sys
    
    def align(self, seq_a: str, seq_b: str, start_gap_penalty_mode=False):
        if start_gap_penalty_mode:
            result = self.execute_with_start_gap_penalty(seq_a, seq_b)
            alignment_a, alignment_b = traceback(
                seq_a, 
                seq_b,
                result_matrix=result['result_matrix'],
                traceback_matrix=result['traceback_matrix'],
                start_pos=result['score_pos'],
                global_alignment=True
            )

            print(
                f"Score={result['score']}\n"
                f"Result:\n {result['result_matrix']}\n"
                f"Traceback:\n {result['traceback_matrix']}\n"
                f"Alignment:\n {alignment_a}\n {alignment_b}\n"
            )
        else:
            self.execute(seq_a, seq_b)
            print(self.aligned_seq_a)
            print(self.aligned_seq_b)

    def execute(self, seq_a: str, seq_b: str):
        # 1. Prepare dimensions (required additional 1 column and 1 row)
        rows, cols = len(seq_a) + 1, len(seq_b) + 1

        # 2. Initialize matrices
        # Use grid/matrix as graph-like acyclic digraph (array cells are vertices)
        H = np.zeros(shape=(rows, cols), dtype=int)
        score_func = self.scoring_sys.score

        # 3. 1st row and column need to have negative values
        # Top row and leftmost column, are like: 0, -1*d, -2*d, -3*d, etc.
        # d is gap penalty
        d = self.scoring_sys.gap
        H[0, :] = np.arange(start=0, stop=d * cols, step=d)
        H[:, 0] = np.arange(start=0, stop=d * rows, step=d)

        for row in range(1, rows):
            for col in range(1, cols):
                # Current pair of letters from sequence A and B
                a = seq_a[row - 1]
                b = seq_b[col - 1]

                leave_or_replace_letter = H[row - 1, col - 1] + score_func(a, b)
                delete_indel = H[row - 1, col] +  score_func('-', b)
                insert_indel = H[row, col - 1] + score_func(a, '-')

                scores = [leave_or_replace_letter, delete_indel, insert_indel]
                best_action = np.argmax(scores)

                H[row, col] = scores[best_action]
        return H
    
    def execute_with_start_gap_penalty(self, seq_a, seq_b):
        '''It's huge mess...'''
        rows, cols = len(seq_a) + 1, len(seq_b) + 1
        H = np.zeros(shape=(rows, cols), dtype=int)
        T = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))
        G = np.zeros(shape=(rows, cols), dtype=int)
        score_func = self.scoring_sys.score

        T[0, 1:] = np.array(list(seq_b), dtype=str)
        T[1:, 0] = np.array(list(seq_a), dtype=str)
    
        d = self.scoring_sys.gap # d is gap penalty
        H[0, :] = np.arange(start=0, stop=d * cols, step=d)
        H[:, 0] = np.arange(start=0, stop=d * rows, step=d)

        for row in range(1, rows):
            for col in range(1, cols):
                # Current pair of letters from sequence A and B
                a = seq_a[row - 1]
                b = seq_b[col - 1]

                leave_or_replace_letter = H[row - 1, col - 1] + score_func(a, b)

                if G[row - 1, col] == 0:
                    score = score_func('-', b)
                else:
                    score = 0
                delete_indel = H[row - 1, col] + score

                if G[row, col - 1] == 0:
                    score = score_func(a, '-')
                else:
                    score = 0
                insert_indel = H[row, col - 1] + score

                scores = [leave_or_replace_letter, delete_indel, insert_indel]

                best_action = np.argmax(scores)
                if best_action in [1, 2]:
                    G[row, col] = 1

                H[row, col] = scores[best_action]
                T[row, col] = traceback_symbols[best_action]

        print(G)
        return {
            'result_matrix': H,
            'traceback_matrix': T,
            'score': H[-1, -1],                 # Always right-bottom corner
            'score_pos': (rows - 1, cols - 1)   # as above...
        }