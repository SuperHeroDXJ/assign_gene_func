def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    n = len(seq1)
    m = len(seq2)
    gap_penalty = -1  # Gap penalty inferred from the doctest example

    # Initialize score matrix and traceback matrix
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[''] * (m + 1) for _ in range(n + 1)]

    # Initialize the first column and first row (filling with gaps)
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap_penalty
        traceback[i][0] = 'U'  # Up
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap_penalty
        traceback[0][j] = 'L'  # Left

    # Fill the score matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Calculate scores for three operations (match/mismatch, deletion, insertion)
            match = score_matrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty

            max_score = max(match, delete, insert)
            score_matrix[i][j] = max_score

            # Record traceback direction (prioritize diagonal, then up, then left in case of a tie)
            if max_score == match:
                traceback[i][j] = 'D'  # Diagonal
            elif max_score == delete:
                traceback[i][j] = 'U'  # Up
            else:
                traceback[i][j] = 'L'  # Left

    # Traceback to find aligned sequences
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n, m

    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 'D':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and traceback[i][j] == 'U':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        elif j > 0 and traceback[i][j] == 'L':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1

    # Traceback goes from the end to the beginning, so reverse the strings at the end
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, float(score_matrix[n][m])


def local_alignment(seq1, seq2, scoring_function, gap_penalty=-1):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    n = len(seq1)
    m = len(seq2)

    # Initialize score matrix and traceback matrix with zeros
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[''] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    # Fill the score matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Calculate scores for three operations
            match = score_matrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty

            # Local alignment resets to 0 if the score is negative
            score = max(0, match, delete, insert)
            score_matrix[i][j] = score

            # Record traceback direction
            if score == 0:
                traceback[i][j] = '0'
            elif score == match:
                traceback[i][j] = 'D'  # Diagonal
            elif score == delete:
                traceback[i][j] = 'U'  # Up
            else:
                traceback[i][j] = 'L'  # Left

            # Keep track of the maximum score and its position
            if score > max_score:
                max_score = score
                max_i, max_j = i, j

    # Traceback from the highest score until we hit a zero score
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_i, max_j

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if traceback[i][j] == 'D':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 'U':
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback[i][j] == 'L':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
        else:
            break

    # Reverse the strings
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, float(max_score)


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
