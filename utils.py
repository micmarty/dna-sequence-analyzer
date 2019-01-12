from typing import Dict, Any, List, Tuple

# Useful for visualization
traceback_symbols = {
    0: '↖',
    1: '↑',
    2: '←',
    3: '•'
}

def traceback(seq_a, seq_b, result_matrix, traceback_matrix, start_pos: Tuple[int, int], global_alignment: bool) -> Tuple[str, str]:
    seq_a_aligned = ''
    seq_b_aligned = ''

    # 1. Select starting point
    row, col = start_pos

    if global_alignment:
        # Terminate when top left corner (0,0) is reached (end of path)
        end_condition_reached = lambda row, col: row == 0 and col == 0
    else:
        # Terminate when 0 is reached
        end_condition_reached = lambda row, col: result_matrix[row, col] == 0

    while not end_condition_reached(row, col):
        symbol = traceback_matrix[row, col]
        if row == 0:
            symbol = '←'
        if col == 0:
            symbol = '↑'
        # Use arrows to navigate and collect letters (in reversed order)
        # Shift/reverse indexes by one beforehand (we want to get the letter that arrow points to)
        if symbol == '↖':
            row -= 1
            col -= 1
            letter_a, letter_b = seq_a[row], seq_b[col]
        elif symbol == '↑':
            row -= 1
            letter_a, letter_b = seq_a[row], '-'
        elif symbol == '←':
            col -= 1
            letter_a, letter_b = '-', seq_b[col]

        # Acumulate letter (in reversed order)
        seq_a_aligned += letter_a
        seq_b_aligned += letter_b

    # Reverse strings (traceback goes from bottom-right to top-left)
    return seq_a_aligned[::-1], seq_b_aligned[::-1]