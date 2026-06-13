#!/usr/bin/env python3
"""
OEIS lookup for sequences discovered in opus-2026-05-28-S3.
Uses OEIS API to check which sequences are already catalogued.
"""

import urllib.request
import urllib.parse
import json
import time

def oeis_search(sequence_str, description=""):
    """Search OEIS for a sequence given as comma-separated values."""
    url = "https://oeis.org/search?fmt=json&q=" + urllib.parse.quote(sequence_str)
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        count = data.get('count', 0)
        results = data.get('results', [])
        return count, results
    except Exception as e:
        return -1, []

def search_and_report(label, sequence_list, desc=""):
    """Search for a sequence and print results."""
    # Use first 6-8 terms for search
    terms = ",".join(str(x) for x in sequence_list[:8])
    print(f"\n{label}: {terms}")
    if desc:
        print(f"  ({desc})")
    count, results = oeis_search(terms)
    if count == -1:
        print("  ERROR: OEIS search failed")
    elif count == 0:
        print("  ** NOT IN OEIS ** (new sequence!)")
    else:
        print(f"  Found {count} result(s):")
        for r in results[:3]:
            aid = r.get('number', '?')
            name = r.get('name', '?')[:80]
            print(f"    A{aid:06d}: {name}")
    time.sleep(1)  # be polite to OEIS
    return count, results

def main():
    print("=" * 70)
    print("OEIS LOOKUP — opus-2026-05-28-S3")
    print("=" * 70)

    # -------------------------------------------------------
    # 1. SC tiling sequence (path-fixed SC tournaments)
    # -------------------------------------------------------
    search_and_report(
        "SC tiling sequence",
        [1, 5, 50, 903, 30773, 2032504, 264271477, 68184627441],
        "#{tilings of n-vertex staircase with SC tournament}"
    )

    # -------------------------------------------------------
    # 2. Non-SC tiling sequence
    # -------------------------------------------------------
    search_and_report(
        "Non-SC tiling sequence",
        [1, 3, 14, 121, 1995, 64648, 4163979, 534849295],
        "#{tilings NOT strongly connected}"
    )

    # -------------------------------------------------------
    # 3. Good-cuts triangle by rows
    # -------------------------------------------------------
    # First non-trivial terms: flatten triangle rows
    triangle_flat = [
        # n=3: [1,0,1]  n=4: [1,0,2,5]  n=5: [1,0,3,10,50]
        1, 0, 1, 1, 0, 2, 5, 1, 0, 3, 10, 50, 1, 0, 4, 15, 101, 903
    ]
    search_and_report(
        "Good-cuts triangle (flat)",
        [1, 0, 1, 1, 0, 2, 5, 1, 0, 3],
        "T(n,k) = #{tilings with exactly k good cuts}"
    )

    # -------------------------------------------------------
    # 4. Column k=2 of good-cuts: n-2
    # -------------------------------------------------------
    search_and_report(
        "Column k=2 of good-cuts (= n-2)",
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "trivial: natural numbers"
    )

    # -------------------------------------------------------
    # 5. Column k=3 of good-cuts: 5(n-3)
    # -------------------------------------------------------
    search_and_report(
        "Column k=3: 5(n-3)",
        [0, 0, 5, 10, 15, 20, 25, 30, 35, 40],
        "multiples of 5"
    )

    # -------------------------------------------------------
    # 6. Column k=4: (n-4)(n+95)/2
    # -------------------------------------------------------
    search_and_report(
        "Column k=4: (n-4)(n+95)/2",
        [0, 0, 0, 50, 101, 153, 206, 260, 315, 371],
        "quadratic in n"
    )

    # -------------------------------------------------------
    # 7. King-count: K(n,1) = n*2^{C(n-1,2)}
    # -------------------------------------------------------
    search_and_report(
        "K(n,1) = #{tournaments with 1 king}",
        [2, 6, 32, 320, 6144, 229376],
        "= n * 2^{C(n-1,2)}"
    )

    # -------------------------------------------------------
    # 8. K(n,n) = #{tournaments with all n kings}
    # -------------------------------------------------------
    search_and_report(
        "K(n,n) = #{tournaments with all n kings}",
        [0, 2, 0, 64, 1680, 110048],
        "= 0 for even n>=4, positive for odd n"
    )

    # -------------------------------------------------------
    # 9. K(n,n) for odd n only
    # -------------------------------------------------------
    search_and_report(
        "K(n,n) for odd n: n=3,5,7",
        [2, 64, 1680, 110048],
        "all-king tournaments, odd n"
    )

    # -------------------------------------------------------
    # 10. 3-cycle distribution at n=4
    # -------------------------------------------------------
    search_and_report(
        "3-cycle dist n=4: [24,16,24]",
        [24, 16, 24],
        "too short, try with n=5 row"
    )
    search_and_report(
        "3-cycle dist n=5",
        [120, 120, 240, 240, 280, 24],
        "C3(5,k) for k=0..5"
    )

    # -------------------------------------------------------
    # 11. H-distribution at n=6, H values
    # -------------------------------------------------------
    search_and_report(
        "H-values at n=6",
        [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45],
        "achievable H in labeled 6-vertex tournament"
    )

    # -------------------------------------------------------
    # 12. H-distribution counts at n=6
    # -------------------------------------------------------
    search_and_report(
        "H-dist counts n=6",
        [720, 960, 2160, 2960, 1440, 1440, 2208],
        "#{6-vertex tournaments with H=1,3,5,9,11,13,15}"
    )

    # -------------------------------------------------------
    # 13. #{distinct H values at n}
    # -------------------------------------------------------
    search_and_report(
        "#{distinct H values at n}",
        [1, 2, 3, 7, 19],
        "n=2..6, count of distinct achievable H values"
    )

    # -------------------------------------------------------
    # 14. SC tiling distribution by upward tiles at n=5
    # -------------------------------------------------------
    search_and_report(
        "SC tilings by upward count n=5",
        [1, 9, 18, 15, 6, 1],
        "for j=1..6: #{SC tilings with j upward tiles}"
    )

    # -------------------------------------------------------
    # 15. SC tilings by upward count at n=6
    # -------------------------------------------------------
    search_and_report(
        "SC tilings by upward count n=6",
        [1, 17, 81, 180, 240, 208, 120, 45, 10, 1],
        "for j=1..10: #{SC tilings of 6-vertex with j upward tiles}"
    )

    # -------------------------------------------------------
    # 16. Q(d,1) = SC sequence (same as seq 1)
    # -------------------------------------------------------
    search_and_report(
        "Q(d,2) convolution values",
        [1, 10, 125, 2306, 73076, 4463038, 552760703],
        "[x^d]A(x)^2 where A(x) = sum SC(b+1)*x^b"
    )

    # -------------------------------------------------------
    # 17. Total 3-cycles across all n-vertex tournaments
    # -------------------------------------------------------
    search_and_report(
        "Total directed 3-cycles in all n-vertex tournaments",
        [2, 64, 2560, 163840, 18350080, 3758096384],
        "= 2^{C(n,2)-2} * C(n,3), n=3..8"
    )

    # -------------------------------------------------------
    # 18. K(n,3) = #{tournaments with exactly 3 kings}
    # -------------------------------------------------------
    search_and_report(
        "K(n,3) = #{tournaments with exactly 3 kings}",
        [2, 32, 520, 11600, 402640],
        "for n=3..7"
    )

    # -------------------------------------------------------
    # 19. Good-cuts last column (SC(n) sequence)
    # -------------------------------------------------------
    search_and_report(
        "Last column of good-cuts triangle",
        [1, 1, 5, 50, 903, 30773, 2032504, 264271477],
        "= SC(n): path-fixed SC tournament counts"
    )

    # -------------------------------------------------------
    # 20. Non-SC ratio numerator: nonSC / 2^{m-n+3}
    # -------------------------------------------------------
    # Numerators when written over 2^m:
    search_and_report(
        "Non-SC counts (starting n=3)",
        [1, 3, 14, 121, 1995, 64648],
        "non-SC tiling counts for n=3..8"
    )

    # -------------------------------------------------------
    # 21. K(n,4) for n=5..7
    # -------------------------------------------------------
    search_and_report(
        "K(n,4) for n=5,6,7",
        [120, 7920, 527520],
        "#{n-vertex tournaments with exactly 4 kings}"
    )

    # -------------------------------------------------------
    # 22. K(n,5) for n=5..7
    # -------------------------------------------------------
    search_and_report(
        "K(n,5) for n=5,6,7",
        [64, 5424, 491568],
        "#{n-vertex tournaments with exactly 5 kings}"
    )

    # -------------------------------------------------------
    # 23. Row sums of king triangle (= 2^{C(n,2)} confirmed)
    # -------------------------------------------------------

    # -------------------------------------------------------
    # 24. SC n=5 upward counts: 1,9,18,15,6,1 — check if A-number
    # -------------------------------------------------------

    # -------------------------------------------------------
    # 25. Maximum H at n (A038375) — check our values match
    # -------------------------------------------------------
    search_and_report(
        "Max H (A038375 check)",
        [1, 1, 3, 5, 15, 45],
        "n=1..6: max Hamiltonian paths in n-vertex tournament"
    )

    # -------------------------------------------------------
    # Summary
    # -------------------------------------------------------
    print("\n" + "=" * 70)
    print("OEIS LOOKUP COMPLETE")
    print("=" * 70)
    print("Sequences marked ** NOT IN OEIS ** are candidates for submission.")

if __name__ == "__main__":
    main()
