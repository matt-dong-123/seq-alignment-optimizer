#!/usr/bin/env python3

import argparse
import alignment


def main():
    parser = argparse.ArgumentParser(description="Alignment algorithms")
    parser.add_argument("file1", help="File with 1st sequence")
    parser.add_argument("file2", help="File with 2nd sequence")
    parser.add_argument("--match", type=float, default=1, help="Match score")
    parser.add_argument("--mismatch", type=float, default=1, help="Mismatch penalty")
    parser.add_argument("--gap", type=float, default=2, help="Gap penalty")
    parser.add_argument(
        "--alg",
        type=str,
        default="Needleman-Wunsch",
        help="Choose alignment algorithm:\nNeedleman-Wunsch\nSmith-Waterman",
    )

    args = parser.parse_args()

    with open(args.file1, "r") as f1, open(args.file2, "r") as f2:
        s1 = f1.read().replace("\n", "").strip()
        s2 = f2.read().replace("\n", "").strip()

    score_fn = alignment.dna_score_fn(args.match, args.mismatch)
    a1, a2 = "", ""
    if args.alg == "Needleman-Wunsch":
        _, a1, a2 = alignment.nw(s1, s2, args.gap, score_fn)
    if args.alg == "Smith-Waterman":
        _, a1, a2 = alignment.sw(s1, s2, args.gap, score_fn)

    print(a1)
    print(a2)


if __name__ == "__main__":
    main()
