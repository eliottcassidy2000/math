# AUDIT a16 (independent auditor, 2026-10-21): spectra from the auditor's census a15_census_n*.txt (n = 3..9).
import math
from collections import Counter
for n in range(3, 10):
    R = [tuple(map(int, l.split())) for l in open(f"a15_census_n{n}.txt")]
    f = math.factorial(n)
    lab = sum(f // a for (_, _, _, a) in R)
    assert all(f % a == 0 and h % a == 0 for (h, _, _, a) in R)
    Hs = Counter(h for h, *_ in R); Cs = Counter(c for _, c, _, _ in R); tc = Counter(h // a for h, _, _, a in R)
    sH = Counter(h for h, c, s, a in R if s)
    maxH = max(Hs); maxC = max(Cs)
    oddmiss = [h for h in range(1, maxH + 1, 2) if h not in Hs]
    # 'ceiling': first missing value after 21
    after = [h for h in oddmiss if h > 21]
    ceiling = after[0] if after else None
    print(f"n={n}: classes {len(R)}, sum n!/|Aut| = {lab} == 2^C(n,2): {lab == 2 ** (n * (n - 1) // 2)}; max H {maxH}, max C {maxC}")
    print(f"   odd H missing below the first missing value > 21 (= {ceiling}): {[h for h in oddmiss if ceiling is None or h < ceiling]}; "
          f"#missing odd H in [{ceiling},{maxH}]: {len(after)}")
    if n in (7, 8, 9):
        print(f"   odd H missing (first 12): {oddmiss[:12]}")
    print(f"   tc = H/|Aut| odd values missing below 60: {[t for t in range(1, 60, 2) if t not in tc]}")
    print(f"   C values missing in [0,{maxC}]: {[c for c in range(maxC + 1) if c not in Cs]}")
    smin = min(sH)
    print(f"   strong H min {smin}; odd values missing in [{smin},{smin + 40}]: {[h for h in range(smin, smin + 41, 2) if h not in sH]}")
    if n == 7:
        print("   tournaments with C = 24:", [(h, c, s, a) for h, c, s, a in R if c == 24])
# compare with the session's census files (H, C, strong, |Aut| columns)
import gzip, os
for n in range(3, 10):
    p = f"../census_n{n}.txt"
    fh = open(p) if os.path.exists(p) else gzip.open(p + ".gz", "rt")
    S = Counter((int(a), int(b), int(c), int(d)) for a, b, c, d, *_ in (l.split() for l in fh))
    A = Counter(tuple(map(int, l.split())) for l in open(f"a15_census_n{n}.txt"))
    print(f"n={n}: multiset of (H,C,strong,|Aut|) identical to session census: {S == A}")
