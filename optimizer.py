from src.alignment import nw, sw


SEARCH_RANGES = {
    "match": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
    "mismatch": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
    "gap": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
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


def grid_search(seq1: str, seq2: str, alg_func, alg_type: str) -> dict:
    best_fitness = -float("inf")
    best = {
        "weights": {"match": 0, "mismatch": 0, "gap": 0},
        "aligned_s1": "",
        "aligned_s2": "",
        "score": 0,
        "stats": {"matches": 0, "mismatches": 0, "gaps": 0, "length": 0, "identity": 0},
        "fitness": 0,
        "alignment_rows": [],
    }

    for match in SEARCH_RANGES["match"]:
        for mismatch in SEARCH_RANGES["mismatch"]:
            for gap in SEARCH_RANGES["gap"]:
                grid, a1, a2 = alg_func(seq1, seq2, match, mismatch, gap)
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
                        "aligned_s1": a1,
                        "aligned_s2": a2,
                        "score": score,
                        "stats": stats,
                        "fitness": fitness,
                        "alignment_rows": format_alignment(a1, a2),
                    }

    return best


def optimize_both(seq1: str, seq2: str) -> dict:
    return {
        "needleman_wunsch": grid_search(seq1, seq2, nw, "nw"),
        "smith_waterman": grid_search(seq1, seq2, sw, "sw"),
    }
