from typing import Any


def nw(s1: str, s2: str, match: int, mismatch: int, gap: int):
    m = len(s1)
    n = len(s2)

    grid = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(1, m + 1):
        grid[i][0] = i * -gap
    for j in range(1, n + 1):
        grid[0][j] = j * -gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = grid[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else -mismatch)
            vert = grid[i - 1][j] - gap
            hori = grid[i][j - 1] - gap
            grid[i][j] = max(diag, vert, hori)

    a1, a2 = "", ""
    i, j = m, n
    while i > 0 and j > 0:
        score = grid[i][j]
        diag = grid[i - 1][j - 1]
        vert = grid[i - 1][j]
        hori = grid[i][j - 1]

        if score == diag + (match if s1[i - 1] == s2[j - 1] else -mismatch):
            a1 = s1[i - 1] + a1
            a2 = s2[j - 1] + a2
            i -= 1
            j -= 1
        elif score == vert - gap:
            a1 = s1[i - 1] + a1
            a2 = "-" + a2
            i -= 1
        elif score == hori - gap:
            a1 = "-" + a1
            a2 = s2[j - 1] + a2
            j -= 1

    while i > 0:
        a1 = s1[i - 1] + a1
        a2 = "-" + a2
        i -= 1
    while j > 0:
        a1 = "-" + a1
        a2 = s2[j - 1] + a2
        j -= 1

    return grid, a1, a2


def sw(s1: str, s2: str, match: int, mismatch: int, gap: int):
    m = len(s1)
    n = len(s2)

    grid = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = grid[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else -mismatch)
            vert = grid[i - 1][j] - gap
            hori = grid[i][j - 1] - gap
            grid[i][j] = max(diag, vert, hori) if max(diag, vert, hori) > 0 else 0

    a1, a2 = "", ""
    i, j = max_2d(grid)
    while i > 0 and j > 0 and grid[i][j] > 0:
        score = grid[i][j]
        diag = grid[i - 1][j - 1]
        vert = grid[i - 1][j]
        hori = grid[i][j - 1]

        if score == diag + (match if s1[i - 1] == s2[j - 1] else -mismatch):
            a1 = s1[i - 1] + a1
            a2 = s2[j - 1] + a2
            i -= 1
            j -= 1
        elif score == vert - gap:
            a1 = s1[i - 1] + a1
            a2 = "-" + a2
            i -= 1
        elif score == hori - gap:
            a1 = "-" + a1
            a2 = s2[j - 1] + a2
            j -= 1

    return grid, a1, a2


def max_2d(arr: list[list[Any]]) -> tuple[int, int]:
    row_idx = arr.index(max(arr, key=max))
    col_idx = arr[row_idx].index(max(arr[row_idx]))
    return row_idx, col_idx
