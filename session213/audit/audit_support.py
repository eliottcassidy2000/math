"""Replay the complete bounded comb support without the author's Gray traversal."""
import json
import sys
from pathlib import Path
from audit_comb import quadratic
from independent_counts import check


def main():
    author = json.loads(Path(sys.argv[1]).read_text())["all_length_comb_support"]
    limit = author["limit"]
    max_m = (limit - 1).bit_length() - 1
    attained, cases = set(), 0
    for m in range(1, max_m + 1):
        for q in range(1 << m):
            for tau in range(2):
                for forward in range(2):
                    h = quadratic(m, q, tau, forward)
                    cases += 1
                    if h <= limit:
                        attained.add(h)
    check(sorted(attained) == author["sorted_support"], "Different complete support")
    print(json.dumps({"status": "PASS", "method": "Naive O(m^2) quadratic evaluation; no Gray-code updates and no author imports",
                      "presentations": cases, "distinct_values": len(attained), "limit": limit,
                      "matching_author_support_sha256": author["sorted_support_sha256"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
