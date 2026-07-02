from typing import Any, Callable


def nw(s1: str, s2: str, gap: float, score_fn: Callable[[str, str], float]):
    m = len(s1)
    n = len(s2)

    grid = [[0.0 for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(1, m + 1):
        grid[i][0] = i * -gap
    for j in range(1, n + 1):
        grid[0][j] = j * -gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = grid[i - 1][j - 1] + score_fn(s1[i - 1], s2[j - 1])
            vert = grid[i - 1][j] - gap
            hori = grid[i][j - 1] - gap
            grid[i][j] = max(diag, vert, hori)

    a1 = []
    a2 = []
    i, j = m, n
    while i > 0 and j > 0:
        score = grid[i][j]
        diag = grid[i - 1][j - 1]
        vert = grid[i - 1][j]
        hori = grid[i][j - 1]

        if score == diag + score_fn(s1[i - 1], s2[j - 1]):
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif score == vert - gap:
            a1.append(s1[i - 1])
            a2.append("-")
            i -= 1
        elif score == hori - gap:
            a1.append("-")
            a2.append(s2[j - 1])
            j -= 1

    while i > 0:
        a1.append(s1[i - 1])
        a2.append("-")
        i -= 1
    while j > 0:
        a1.append("-")
        a2.append(s2[j - 1])
        j -= 1

    return grid, "".join(reversed(a1)), "".join(reversed(a2))


def sw(s1: str, s2: str, gap: float, score_fn: Callable[[str, str], float]):
    m = len(s1)
    n = len(s2)

    grid = [[0.0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = grid[i - 1][j - 1] + score_fn(s1[i - 1], s2[j - 1])
            vert = grid[i - 1][j] - gap
            hori = grid[i][j - 1] - gap
            best = max(diag, vert, hori)
            grid[i][j] = best if best > 0 else 0

    max_i, max_j = max_2d(grid)
    seg1, seg2 = [], []
    i, j = max_i, max_j
    while i > 0 and j > 0 and grid[i][j] > 0:
        score = grid[i][j]
        diag = grid[i - 1][j - 1]
        vert = grid[i - 1][j]
        hori = grid[i][j - 1]

        if score == diag + score_fn(s1[i - 1], s2[j - 1]):
            seg1.append(s1[i - 1])
            seg2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif score == vert - gap:
            seg1.append(s1[i - 1])
            seg2.append("-")
            i -= 1
        elif score == hori - gap:
            seg1.append("-")
            seg2.append(s2[j - 1])
            j -= 1

    seg1.reverse()
    seg2.reverse()
    pref1, pref2 = s1[:i], s2[:j]
    suff1, suff2 = s1[max_i:], s2[max_j:]
    full1, full2 = [], []
    lp1, lp2 = len(pref1), len(pref2)
    if lp1 > lp2:
        full1.extend(pref1[: lp1 - lp2])
        full2.extend("-" * (lp1 - lp2))
        full1.extend(pref1[lp1 - lp2 :])
        full2.extend(pref2)
    else:
        full1.extend("-" * (lp2 - lp1))
        full2.extend(pref2[: lp2 - lp1])
        full1.extend(pref1)
        full2.extend(pref2[lp2 - lp1 :])
    full1.extend(seg1)
    full2.extend(seg2)
    ls1, ls2 = len(suff1), len(suff2)
    if ls1 > ls2:
        full1.extend(suff1[:ls2])
        full2.extend(suff2)
        full1.extend(suff1[ls2:])
        full2.extend("-" * (ls1 - ls2))
    else:
        full1.extend(suff1)
        full2.extend(suff2[:ls1])
        full1.extend("-" * (ls2 - ls1))
        full2.extend(suff2[ls1:])
    return grid, "".join(full1), "".join(full2)


def dna_score_fn(match: float, mismatch: float) -> Callable[[str, str], float]:
    return lambda a, b: match if a == b else -mismatch


def protein_score_fn(name: str) -> Callable[[str, str], float]:
    from src.matrices import get_matrix

    matrix = get_matrix(name)
    return lambda a, b: matrix.get((a, b), matrix.get((b, a), -4))


def max_2d(arr: list[list[Any]]) -> tuple[int, int]:
    row_idx = arr.index(max(arr, key=max))
    col_idx = arr[row_idx].index(max(arr[row_idx]))
    return row_idx, col_idx
