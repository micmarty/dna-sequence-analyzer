import pandas as pd

class ScoringSystem:
    def __init__(self, match: int=1, mismatch: int=-1, gap: int=-1, load_from_csv: bool=False) -> None:
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.custom_scoring = None
        if load_from_csv:
            # TODO File name cannot be hardcoded
            self.custom_scoring = pd.read_csv('scores.csv', header=0, index_col=0, sep=' ')

    def _default_scoring(self, a: str, b: str) -> int:
        if a == b:
            return self.match
        elif a == '-' or b == '-':
            return self.gap
        return self.mismatch

    def score(self, a: str, b: str) -> int:
        assert isinstance(a, str) and isinstance(b, str)
        assert len(a) == 1 and len(b) == 1

        if self.custom_scoring is not None:
            try:
                return self.custom_scoring[a][b]
            except KeyError:
                return self._default_scoring(a, b)
        
    def __str__(self):
        return f'Match: {self.match}, Mismatch: {self.mismatch}, Gap: {self.gap}'

ScoringSystem(0,0,0).load_from_file()