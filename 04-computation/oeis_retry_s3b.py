#!/usr/bin/env python3
"""
Retry OEIS lookups for sequences from opus-2026-05-28 sessions.
Previous attempt (S3) got all network errors. Try again.
"""
import urllib.request, urllib.parse, json, time

def oeis_search(seq_str):
    url = "https://oeis.org/search?fmt=json&q=" + urllib.parse.quote(seq_str)
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read().decode())
        return data.get('count', 0), data.get('results', [])
    except Exception as e:
        return -1, [str(e)]

def check(label, seq, desc=""):
    terms = ",".join(str(x) for x in seq[:8])
    print(f"\n{label}")
    print(f"  Seq: {terms}")
    if desc: print(f"  Desc: {desc}")
    count, results = oeis_search(terms)
    if count == -1:
        print(f"  ERROR: {results[0] if results else 'unknown'}")
    elif count == 0:
        print(f"  *** NOT IN OEIS *** — new sequence candidate!")
    else:
        print(f"  Found {count} match(es):")
        for r in results[:2]:
            aid = r.get('number', '?')
            name = r.get('name', '?')[:80]
            print(f"    A{aid:06d}: {name}")
    time.sleep(1.5)
    return count, results

print("=" * 65)
print("OEIS LOOKUP RETRY — opus-2026-05-28-S3b")
print("=" * 65)

# ── Core new sequences ──────────────────────────────────────

check("SC tiling sequence (n=3..10)",
    [1, 5, 50, 903, 30773, 2032504, 264271477, 68184627441],
    "Path-fixed SC tournaments: strongly connected tilings of n-vertex staircase")

check("Non-SC tiling sequence (n=3..10)",
    [1, 3, 14, 121, 1995, 64648, 4163979, 534849295],
    "Non-SC tilings: complement of SC tilings in 2^{C(n-1,2)}")

# ── Q(d,k) triangle — new candidate ─────────────────────────
# Flatten by rows: d=2:[1], d=3:[5], d=4:[50,1], d=5:[903,10], ...
q_flat = [1, 5, 50, 1, 903, 10, 30773, 125, 1, 2032504, 2306, 15]
check("Q(d,k) triangle (flat by rows)",
    q_flat[:9],
    "Q(d,k)=[x^d]B(x)^k where B(x)=sum SC(n+1)x^n; k=1..floor(d/2)")

# ── Q(d,2) column — convolution of SC with itself ───────────
check("Q(d,2) column = [x^d]B(x)^2",
    [1, 10, 125, 2306, 73076, 4463038, 552760703],
    "d=4,5,6,7,8,9,10: [x^d](sum_{a>=2} SC(a+1)x^a)^2")

# ── King distribution sequences ─────────────────────────────
check("K(n,1) = #{n-tourn with 1 king}",
    [2, 6, 32, 320, 6144, 229376],
    "= n * 2^{C(n-1,2)}, n=2..7")

check("K(n,n) = #{n-tourn with all kings}",
    [0, 2, 0, 64, 1680, 110048],
    "n=2..7; zero for even n>=4, positive for odd n")

check("K(n,3) for n=3..7",
    [2, 32, 520, 11600, 402640],
    "#{n-vertex tournaments with exactly 3 kings}")

# ── H-spectrum sequences ─────────────────────────────────────
check("#{distinct H values at n}",
    [1, 2, 3, 7, 19],
    "n=2..6, count of achievable H(T) values")

check("H-dist counts at n=6",
    [720, 960, 2160, 2960, 1440, 1440, 2208],
    "n=6: counts for H=1,3,5,9,11,13,15")

# ── SC tilings by upward tile count ─────────────────────────
check("SC tilings n=5 by upward count",
    [1, 9, 18, 15, 6, 1],
    "n=5: #{SC tilings with j=1..6 upward tiles}")

check("SC tilings n=6 by upward count",
    [1, 17, 81, 180, 240, 208, 120, 45, 10, 1],
    "n=6: #{SC tilings with j=1..10 upward tiles}")

# ── Known total directed k-cycle count (should hit OEIS) ────
check("Total 3-cycles in all n-tourn",
    [2, 64, 2560, 163840, 18350080, 3758096384],
    "= 2^{C(n,2)-2}*C(n,3), n=3..8")

# ── Good-cuts triangle column 1 (= n-2) ─────────────────────
# This is trivial: just 0,1,2,3,...
# Skip.

# ── SC tiling upward-count triangle ─────────────────────────
# n=5: [0,1,9,18,15,6,1], n=6: [0,1,17,81,...]
# The first element of each row (j=1) is always 1 (only 1 SC tiling with 1 upward tile)
# The last element (j=m) is always 1 (only 1 SC tiling with all m upward tiles = transitive-complement)
# The j=2 column: 9, 17, ... ?
check("SC upward-count col j=2 (n=5,6,7...)",
    [9, 17],
    "#{SC tilings with exactly 2 upward tiles}: at n=5: 9, n=6: 17")

# ── Q-triangle sub-diagonal Q(2k+1,k) = 5k ─────────────────
check("Q sub-diagonal = 5k",
    [5, 10, 15, 20, 25, 30, 35, 40],
    "Q(2k+1,k) = 5k for k=1,2,...; multiples of 5")

# ── Q-triangle 2nd sub-diagonal Q(2k+2,k) = 25k(k+3)/2 ─────
check("Q 2nd sub-diagonal",
    [50, 125, 225, 350, 500, 675],
    "Q(2k+2,k) = 25k(k+3)/2 for k=1,2,...")

print("\n" + "=" * 65)
print("DONE")
