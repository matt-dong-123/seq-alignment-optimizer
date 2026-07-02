import math

from src.alignment import nw, sw, dna_score_fn, protein_score_fn

# match is fixed at 1.0 (scaling is irrelevant — DP max() only cares about ratios).
# Only the mismatch:match and gap:match ratios matter.
DNA_MATCH_FIXED = 1.0
DNA_MISMATCH_BOUNDS = (0.5, 5.0)
DNA_GAP_BOUNDS = (0.5, 10.0)

PROTEIN_GAP_BOUNDS = (4.0, 24.0)


def alignment_stats(a1: str, a2: str) -> dict:
    matches = sum(1 for c1, c2 in zip(a1, a2) if c1 == c2 and c1 != "-")
    mismatches = sum(
        1 for c1, c2 in zip(a1, a2) if c1 != c2 and c1 != "-" and c2 != "-"
    )
    gaps = sum(1 for c1, c2 in zip(a1, a2) if c1 == "-" or c2 == "-")
    length = len(a1)
    return {
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "length": length,
        "identity": round(matches / length * 100, 2) if length > 0 else 0,
    }


def match_indicator(a1: str, a2: str) -> str:
    chars = []
    for c1, c2 in zip(a1, a2):
        if c1 == "-" or c2 == "-":
            chars.append(" ")
        elif c1 == c2:
            chars.append("|")
        else:
            chars.append(".")
    return "".join(chars)


def format_alignment(a1: str, a2: str, width: int = 60) -> list[str]:
    indicator = match_indicator(a1, a2)
    rows = []
    for i in range(0, len(a1), width):
        rows.append(a1[i : i + width])
        rows.append(indicator[i : i + width])
        rows.append(a2[i : i + width])
        rows.append("")
    return rows


def run_with_weights(
    seq1: str,
    seq2: str,
    gap: float,
    score_fn,
    seq_type: str = "dna",
    matrix_name: str | None = None,
    match: float | None = None,
    mismatch: float | None = None,
) -> dict:
    def build_result(alg_func, alg_type):
        grid, a1, a2 = alg_func(seq1, seq2, gap, score_fn)
        stats = alignment_stats(a1, a2)
        if alg_type == "nw":
            score = grid[-1][-1]
        else:
            score = max(max(row) for row in grid)
        if seq_type == "dna":
            weights = {"match": match, "mismatch": mismatch, "gap": gap}
        else:
            weights = {"gap": gap}
        return {
            "weights": weights,
            "matrix_name": matrix_name,
            "aligned_s1": a1,
            "aligned_s2": a2,
            "score": score,
            "stats": stats,
            "fitness": stats["matches"] - stats["gaps"],
            "alignment_rows": format_alignment(a1, a2),
        }

    return {
        "needleman_wunsch": build_result(nw, "nw"),
        "smith_waterman": build_result(sw, "sw"),
    }


def build_dna_score_fn(match, mismatch):
    return dna_score_fn(match, mismatch)


def build_protein_score_fn(matrix_name):
    return protein_score_fn(matrix_name)


def _golden_section_max(fn, left, right, n_iter):
    phi = (math.sqrt(5) - 1) / 2
    a, b = left, right
    c = b - phi * (b - a)
    d = a + phi * (b - a)
    fc = fn(c)
    fd = fn(d)
    for _ in range(n_iter - 1):
        if fc > fd:
            b = d
            d = c
            fd = fc
            c = b - phi * (b - a)
            fc = fn(c)
        else:
            a = c
            c = d
            fc = fd
            d = a + phi * (b - a)
            fd = fn(d)
    return (a + b) / 2 if fc > fd else (c + d) / 2


def _nw_score_dna(seq1, seq2, match, mismatch, gap):
    m, n = len(seq1), len(seq2)
    grid = [[0.0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        grid[i][0] = i * -gap
    for j in range(1, n + 1):
        grid[0][j] = j * -gap
    for i in range(1, m + 1):
        gi, gi_1, si = grid[i], grid[i - 1], seq1[i - 1]
        for j in range(1, n + 1):
            diag = gi_1[j - 1] + (match if si == seq2[j - 1] else -mismatch)
            vert = gi_1[j] - gap
            hori = gi[j - 1] - gap
            gi[j] = diag if diag > vert else vert
            if hori > gi[j]:
                gi[j] = hori
    return grid[-1][-1]


def coordinate_descent_dna(seq1, seq2):
    match = DNA_MATCH_FIXED
    mismatch = 2.0
    gap = 3.0
    n_iter = 3

    mismatch = _golden_section_max(
        lambda m: _nw_score_dna(seq1, seq2, match, m, gap),
        *DNA_MISMATCH_BOUNDS,
        n_iter,
    )
    gap = _golden_section_max(
        lambda g: _nw_score_dna(seq1, seq2, match, mismatch, g),
        *DNA_GAP_BOUNDS,
        n_iter,
    )

    def nw_full():
        m, n = len(seq1), len(seq2)
        grid = [[0.0] * (n + 1) for _ in range(m + 1)]
        for i in range(1, m + 1):
            grid[i][0] = i * -gap
        for j in range(1, n + 1):
            grid[0][j] = j * -gap
        for i in range(1, m + 1):
            gi, gi_1, si = grid[i], grid[i - 1], seq1[i - 1]
            for j in range(1, n + 1):
                diag = gi_1[j - 1] + (match if si == seq2[j - 1] else -mismatch)
                vert = gi_1[j] - gap
                hori = gi[j - 1] - gap
                gi[j] = diag if diag > vert else vert
                if hori > gi[j]:
                    gi[j] = hori
        a1, a2 = [], []
        i, j = m, n
        while i > 0 and j > 0:
            diag = grid[i - 1][j - 1]
            vert = grid[i - 1][j]
            hori = grid[i][j - 1]
            if grid[i][j] == diag + (
                match if seq1[i - 1] == seq2[j - 1] else -mismatch
            ):
                a1.append(seq1[i - 1])
                a2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif grid[i][j] == vert - gap:
                a1.append(seq1[i - 1])
                a2.append("-")
                i -= 1
            else:
                a1.append("-")
                a2.append(seq2[j - 1])
                j -= 1
        while i > 0:
            a1.append(seq1[i - 1])
            a2.append("-")
            i -= 1
        while j > 0:
            a1.append("-")
            a2.append(seq2[j - 1])
            j -= 1
        return grid[-1][-1], "".join(reversed(a1)), "".join(reversed(a2))

    def sw_full():
        m, n = len(seq1), len(seq2)
        grid = [[0.0] * (n + 1) for _ in range(m + 1)]
        max_val, max_i, max_j = 0.0, 0, 0
        for i in range(1, m + 1):
            gi, gi_1, si = grid[i], grid[i - 1], seq1[i - 1]
            for j in range(1, n + 1):
                diag = gi_1[j - 1] + (match if si == seq2[j - 1] else -mismatch)
                vert = gi_1[j] - gap
                hori = gi[j - 1] - gap
                best = diag if diag > vert else vert
                if hori > best:
                    best = hori
                best = best if best > 0 else 0
                gi[j] = best
                if best > max_val:
                    max_val, max_i, max_j = best, i, j
        a1, a2 = [], []
        i, j = max_i, max_j
        while i > 0 and j > 0 and grid[i][j] > 0:
            diag = grid[i - 1][j - 1]
            vert = grid[i - 1][j]
            hori = grid[i][j - 1]
            if grid[i][j] == diag + (
                match if seq1[i - 1] == seq2[j - 1] else -mismatch
            ):
                a1.append(seq1[i - 1])
                a2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif grid[i][j] == vert - gap:
                a1.append(seq1[i - 1])
                a2.append("-")
                i -= 1
            else:
                a1.append("-")
                a2.append(seq2[j - 1])
                j -= 1
        return max_val, "".join(reversed(a1)), "".join(reversed(a2))

    nw_score, nw_a1, nw_a2 = nw_full()
    sw_score, sw_a1, sw_a2 = sw_full()

    def p(alg, score, a1, a2):
        stats = alignment_stats(a1, a2)
        return {
            "weights": {"match": match, "mismatch": mismatch, "gap": gap},
            "matrix_name": None,
            "aligned_s1": a1,
            "aligned_s2": a2,
            "score": score,
            "stats": stats,
            "fitness": stats["matches"] - stats["gaps"],
            "alignment_rows": format_alignment(a1, a2),
        }

    return {
        "needleman_wunsch": p("nw", nw_score, nw_a1, nw_a2),
        "smith_waterman": p("sw", sw_score, sw_a1, sw_a2),
    }


def coordinate_descent_protein(seq1, seq2, alg_func, alg_type, matrix_name):
    score_fn = build_protein_score_fn(matrix_name)
    gap = 12.0

    def score_for(g):
        grid, a1, a2 = alg_func(seq1, seq2, g, score_fn)
        if alg_type == "nw":
            return grid[-1][-1]
        return max(max(row) for row in grid)

    gap = _golden_section_max(score_for, *PROTEIN_GAP_BOUNDS, n_iter=10)
    grid, a1, a2 = alg_func(seq1, seq2, gap, score_fn)
    stats = alignment_stats(a1, a2)
    if alg_type == "nw":
        score = grid[-1][-1]
    else:
        score = max(max(row) for row in grid)
    return {
        "weights": {"gap": gap},
        "matrix_name": matrix_name,
        "aligned_s1": a1,
        "aligned_s2": a2,
        "score": score,
        "stats": stats,
        "fitness": stats["matches"] - stats["gaps"],
        "alignment_rows": format_alignment(a1, a2),
    }


def optimize_both(
    seq1: str,
    seq2: str,
    match: float | None = None,
    mismatch: float | None = None,
    gap: float | None = None,
    seq_type: str = "dna",
    matrix_name: str | None = None,
) -> dict:
    if seq_type == "dna":
        if match is not None and mismatch is not None and gap is not None:
            score_fn = build_dna_score_fn(match, mismatch)
            return run_with_weights(
                seq1,
                seq2,
                gap,
                score_fn,
                seq_type="dna",
                match=match,
                mismatch=mismatch,
            )
        return coordinate_descent_dna(seq1, seq2)
    else:
        if gap is not None:
            score_fn = build_protein_score_fn(matrix_name or "BLOSUM62")
            return run_with_weights(
                seq1,
                seq2,
                gap,
                score_fn,
                seq_type="protein",
                matrix_name=matrix_name or "BLOSUM62",
            )
        return {
            "needleman_wunsch": coordinate_descent_protein(
                seq1, seq2, nw, "nw", matrix_name or "BLOSUM62"
            ),
            "smith_waterman": coordinate_descent_protein(
                seq1, seq2, sw, "sw", matrix_name or "BLOSUM62"
            ),
        }
