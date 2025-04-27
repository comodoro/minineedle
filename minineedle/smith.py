from typing import Any, Callable, Optional, Sequence

from minineedle.core import OptimalAlignment
from minineedle.typesvars import ItemToAlign, UP, LEFT, DIAG, NONE


class SmithWaterman(OptimalAlignment[ItemToAlign]):
    """
    Smith-Waterman algorithm
    """

    def __init__(self, seq1: Sequence[ItemToAlign], seq2: Sequence[ItemToAlign], similarity_func: Optional[Callable[[Any, Any], float]]=None) -> None:
        super().__init__(seq1, seq2, similarity_func)

    def _add_gap_penalties(self) -> None:
        """
        Fills number matrix first row and first column with the gap penalties.
        """
        for i in range(1, len(self.seq1) + 1):
            self._nmatrix[0][i] = 0

        for j in range(1, len(self.seq2) + 1):
            self._nmatrix[j][0] = 0

    def _get_last_cell_position(self) -> tuple[int, int]:
        """
        Returns the cell row and column of the last cell in the matrix in which
        the alignment ends. For Needleman-Wunsch this will be the last cell of the matrix,
        for Smith-Waterman will be the cell with the highest score.
        """
        imax, jmax = 0, 0
        max_score = 0
        for irow in range(0, len(self._nmatrix)):
            for jcol in range(0, len(self._nmatrix[0])):
                score = self._nmatrix[irow][jcol]
                if score > max_score:
                    imax = irow
                    jmax = jcol
                    max_score = score
        return imax, jmax

    def _check_best_score(self, diagscore: int, topscore: int, leftscore: int, irow: int, jcol: int) -> None:
        best_pointer: Optional[str] = NONE
        best_score = 0

        if diagscore >= topscore:
            if diagscore >= leftscore:
                best_pointer, best_score = (DIAG, diagscore)
            else:
                best_pointer, best_score = (LEFT, leftscore)
        else:
            if topscore > leftscore:
                best_pointer, best_score = (UP, topscore)
            else:
                best_pointer, best_score = (LEFT, leftscore)

        if best_score < 0:
            best_pointer = None
            best_score = 0

        self._pmatrix[irow + 1][jcol + 1] = best_pointer
        self._nmatrix[irow + 1][jcol + 1] = best_score
