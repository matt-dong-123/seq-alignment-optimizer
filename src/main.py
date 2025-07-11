import argparse
import numpy as np


def nw(s1: str, s2: str, match: int, mismatch: int, gap: int):
    m = len(s1)
    n = len(s2)

    grid = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(1, m + 1):
        grid[i, 0] = i * gap
    for j in range(1, n + 1):
        grid[0, j] = j * gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = grid[i - 1, j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            vert = grid[i - 1, j] + gap
            hori = grid[i, j - 1] + gap
            grid[i, j] = max(diag, vert, hori)

    a1, a2 = "", ""
    i, j = m, n
    while i > 0 and j > 0:
        score = grid[i, j]
        diag = grid[i - 1, j - 1]
        vert = grid[i - 1, j]
        hori = grid[i, j - 1]

        if score == diag + (match if s1[i - 1] == s2[j - 1] else mismatch):
            a1 = s1[i - 1] + a1
            a2 = s2[j - 1] + a2
            i -= 1
            j -= 1
        elif score == vert + gap:
            a1 = s1[i - 1] + a1
            a2 = "-" + a2
            i -= 1
        elif score == hori + gap:
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


def main():
    parser = argparse.ArgumentParser(
        description="Needleman-Wunsch algorithm implementation"
    )
    parser.add_argument("file1", help="File containing sequence 1")
    parser.add_argument("file2", help="File containing sequence 2")
    parser.add_argument("--match", type=int, default=1, help="Match score")
    parser.add_argument("--mismatch", type=int, default=-1, help="Mismatch score")
    parser.add_argument("--gap", type=int, default=-1, help="Gap score")
    args = parser.parse_args()

    with open(args.file1, "r") as f1, open(args.file2, "r") as f2:
        s1 = f1.read().replace("\n", "").strip()
        s2 = f2.read().replace("\n", "").strip()

    grid, a1, a2 = nw(s1, s2, args.match, args.mismatch, args.gap)
    print(grid[-1][-1])
    print(a1)
    print(a2)


if __name__ == "__main__":
    main()
