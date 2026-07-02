from src.alignment import nw, sw, dna_score_fn, protein_score_fn

DNA_SEARCH_RANGES = {
    "match": [1.0,2.0],
    "mismatch": [1.0,2.0],
    "gap": [2.0,3.0],
}

PROTEIN_SEARCH_RANGES = {
    "gap": [4, 6, 8, 10, 12, 14, 16, 18, 20],
}


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


def grid_search_dna(seq1: str, seq2: str, alg_func, alg_type: str) -> dict:
    best_fitness = -float("inf")
    best = {
        "weights": {"match": 0, "mismatch": 0, "gap": 0},
        "matrix_name": None,
        "aligned_s1": "",
        "aligned_s2": "",
        "score": 0,
        "stats": {"matches": 0, "mismatches": 0, "gaps": 0, "length": 0, "identity": 0},
        "fitness": 0,
        "alignment_rows": [],
    }

    for match in DNA_SEARCH_RANGES["match"]:
        for mismatch in DNA_SEARCH_RANGES["mismatch"]:
            for gap in DNA_SEARCH_RANGES["gap"]:
                score_fn = build_dna_score_fn(match, mismatch)
                grid, a1, a2 = alg_func(seq1, seq2, gap, score_fn)
                stats = alignment_stats(a1, a2)

                if alg_type == "nw":
                    score = grid[-1][-1]
                else:
                    score = max(max(row) for row in grid)

                fitness = stats["matches"] - stats["gaps"]

                if fitness > best_fitness:
                    best_fitness = fitness
                    best = {
                        "weights": {
                            "match": match,
                            "mismatch": mismatch,
                            "gap": gap,
                        },
                        "matrix_name": None,
                        "aligned_s1": a1,
                        "aligned_s2": a2,
                        "score": score,
                        "stats": stats,
                        "fitness": fitness,
                        "alignment_rows": format_alignment(a1, a2),
                    }

    return best


def grid_search_protein(
    seq1: str, seq2: str, alg_func, alg_type: str, matrix_name: str
) -> dict:
    score_fn = build_protein_score_fn(matrix_name)
    best_fitness = -float("inf")
    best = {
        "weights": {"gap": 0},
        "matrix_name": matrix_name,
        "aligned_s1": "",
        "aligned_s2": "",
        "score": 0,
        "stats": {"matches": 0, "mismatches": 0, "gaps": 0, "length": 0, "identity": 0},
        "fitness": 0,
        "alignment_rows": [],
    }

    for gap in PROTEIN_SEARCH_RANGES["gap"]:
        grid, a1, a2 = alg_func(seq1, seq2, gap, score_fn)
        stats = alignment_stats(a1, a2)

        if alg_type == "nw":
            score = grid[-1][-1]
        else:
            score = max(max(row) for row in grid)

        fitness = stats["matches"] - stats["gaps"]

        if fitness > best_fitness:
            best_fitness = fitness
            best = {
                "weights": {"gap": gap},
                "matrix_name": matrix_name,
                "aligned_s1": a1,
                "aligned_s2": a2,
                "score": score,
                "stats": stats,
                "fitness": fitness,
                "alignment_rows": format_alignment(a1, a2),
            }

    return best


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
        return {
            "needleman_wunsch": grid_search_dna(seq1, seq2, nw, "nw"),
            "smith_waterman": grid_search_dna(seq1, seq2, sw, "sw"),
        }
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
            "needleman_wunsch": grid_search_protein(
                seq1, seq2, nw, "nw", matrix_name or "BLOSUM62"
            ),
            "smith_waterman": grid_search_protein(
                seq1, seq2, sw, "sw", matrix_name or "BLOSUM62"
            ),
        }
