#!/usr/bin/env python3
"""
shadow_lattice_gap_rank_klein_S313.py — klein-2026-07-15-S313 (cont.2)

THE (r, g) SHADOW LATTICE: rank-r Vandermonde truncations of the figurate table,
summed along gap-g diagonals — with the polyhedral/polygonal case as (r, g) = (inf, *)/(1, *).

  Rank-r figurate table:  T_r(k, j) = sum_{i=0..r} C(k-1, i) C(j, i+1)   (col 0 := all ones)
      r = 1  -> polygonal numbers,  r >= k-1 -> simplex numbers (Pascal columns).
  Gap-g sum: D_{r,g}(d) = sum_{k>=0} T_r(k, d-(g+1)k+1)   (g=0: rows, g=1: shallow diagonals).

MASTER THEOREM (one-line proof, verified here):
  full_g(x) = 1/(1-x-x^{g+1})  (g-bonacci: 2^n, Fibonacci, Narayana, ...)
  D_{r,g}(x) = [N_{r,g}(x)] / [(1-x)^{r+2} (1-x^{g+1})^{r+1} (1-x-x^{g+1})],
      N_{r,g}(x) = (1-x)^{r+2} (1-x^{g+1})^{r+1} - x^{(g+2)(r+2)-1},
  and N_{r,g} is DIVISIBLE by the g-bonacci kernel 1-x-x^{g+1}:
      at any kernel root, x^{g+1} = 1-x, so 1-x^{g+1} = x and
      x^{(g+2)(r+2)-1} = (x^{g+1})^{r+2} x^{r+1} = (1-x)^{r+2} x^{r+1}.  QED
  Hence D_{r,g} is a quasi-polynomial of degree 2r+2 with period g+1.

MISSING-REGION LAW: deficit_g(d) = full_g(d) - D_{r,g}(d) is 0 for d < d0 and EXACTLY 1 at
  d0 = (g+2)(r+2)-1.  (r,g)=(1,0): d0 = 5 -> Moser's missing 32nd region.  (1,1): d0 = 8
  (34 vs 33).  Every (r,g) has its own "32nd region" moment, always of deficit exactly 1.

TOURNAMENT SIDE (the +8 climb):
  x(T) = sum d_v^2 = (n^3-n)/3 - 8*c3(T)  where c3 = # cyclic triangles (Kendall).
  So THM-866's k(T) = ((n^3-n)/3 - x)/8 = c3(T): the +8 quantum IS one cyclic triangle,
  and every tie-splitting flip destroys EXACTLY ONE net 3-cycle.  Also: |Aut(T)| is odd
  (exhaustive n<=5), the n=5 floor is the round tournament C5 with c3 = 5 = k_max(5).
"""
from math import comb
from fractions import Fraction
import itertools, random

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

# ---------- the (r,g) lattice ----------
def T_r(r, k, j):
    if k == 0: return 1
    return sum(comb(k - 1, i) * comb(j, i + 1) for i in range(0, min(r, k - 1) + 1))

def D_rg(r, g, d):
    return sum(T_r(r, k, d - (g + 1) * k + 1) for k in range(0, d // (g + 1) + 1))

def full_g(g, N):        # g-bonacci: a(d) = a(d-1) + a(d-g-1), GF 1/(1-x-x^{g+1})
    a = []
    for d in range(N):
        a.append(1 if d == 0 else a[d - 1] + (a[d - g - 1] if d >= g + 1 else 0))
    return a

N = 64
# consistency with part 1 (r=1, g=1 = the D sequence; r=1, g=0 = Moser)
D11 = [D_rg(1, 1, d) for d in range(20)]
check("(r,g)=(1,1) reproduces D = 1,1,2,3,5,8,13,21,33,51,...",
      D11 == [1, 1, 2, 3, 5, 8, 13, 21, 33, 51, 76, 111, 157, 218, 295, 393, 513, 661, 838, 1051])
check("(r,g)=(1,0) reproduces Moser row sums", all(
    D_rg(1, 0, n) == 1 + comb(n + 1, 2) + comb(n + 1, 4) for n in range(30)))

# polyhedral side: r = inf (use r = d) gives the g-bonacci family exactly
check("r=inf gap-g sums = g-bonacci (2^n, Fibonacci, Narayana A000930, A003269, A003520)",
      all([D_rg(d + 2, g, d) for d in range(N)] == full_g(g, N) for g in range(0, 5)))

# ---------- master theorem: kernel divisibility, all small (r,g) ----------
def polymul(A, B):
    C_ = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            C_[i + j] += a * b
    return C_
def polypow(A, e):
    out = [1]
    for _ in range(e): out = polymul(out, A)
    return out
def divisible_by_kernel(num, g):
    # divide by 1 - x - x^{g+1}
    r_ = num[:] + [0] * (g + 2)
    for i in range(len(r_) - (g + 1)):
        c = r_[i]
        if c:
            r_[i] = 0; r_[i + 1] += c; r_[i + g + 1] += c
    return all(v == 0 for v in r_)

ok_div = True
for r in range(1, 5):
    for g in range(0, 6):
        one_x = [1, -1]
        one_xg = [1] + [0] * g + [-1]
        num = polymul(polypow(one_x, r + 2), polypow(one_xg, r + 1))
        e = (g + 2) * (r + 2) - 1
        num = num + [0] * max(0, e + 1 - len(num))
        num[e] -= 1
        if not divisible_by_kernel(num, g): ok_div = False
check("MASTER: (1-x)^{r+2}(1-x^{g+1})^{r+1} - x^{(g+2)(r+2)-1} divisible by 1-x-x^{g+1}, r<=4, g<=5",
      ok_div)

# GF identity: D_{r,g} = full_g - tail, tail = x^{(g+2)(r+2)-1}/[(1-x^{g+1})^{r+1}(1-x)^{r+2}(1-x-x^{g+1})]
def series_from_gf(num, den, L):
    a = [Fraction(0)] * L
    for d in range(L):
        s = Fraction(num[d] if d < len(num) else 0)
        for t in range(1, d + 1):
            if t < len(den):
                s -= den[t] * a[d - t]
        a[d] = s / den[0]
    return a
ok_gf = True
for r in range(1, 4):
    for g in range(0, 4):
        one_x = [1, -1]; one_xg = [1] + [0] * g + [-1]
        ker = [1] + [0] * (g + 1)
        ker[1] -= 1; ker[g + 1] -= 1          # 1 - x - x^{g+1}  (g=0 -> 1 - 2x)
        den = polymul(polymul(polypow(one_x, r + 2), polypow(one_xg, r + 1)), ker)
        e = (g + 2) * (r + 2) - 1
        num = polymul(polypow(one_x, r + 2), polypow(one_xg, r + 1))
        num = num + [0] * max(0, e + 1 - len(num)); num[e] -= 1
        ser = series_from_gf(num, den, 40)
        if [D_rg(r, g, d) for d in range(40)] != [int(v) for v in ser]: ok_gf = False
check("GF identity: D_{r,g}(x) = N_{r,g}/[(1-x)^{r+2}(1-x^{g+1})^{r+1}(1-x-x^{g+1})], r<=3, g<=3", ok_gf)

# ---------- missing-region law ----------
ok_def = True
first_table = []
for r in range(1, 4):
    for g in range(0, 5):
        fg = full_g(g, N)
        defc = [fg[d] - D_rg(r, g, d) for d in range(N)]
        d0 = (g + 2) * (r + 2) - 1
        if not (all(v == 0 for v in defc[:d0]) and defc[d0] == 1): ok_def = False
        if r == 1: first_table.append((g, d0, fg[d0], fg[d0] - 1))
check("MISSING-REGION LAW: deficit 0 before d0=(g+2)(r+2)-1 and EXACTLY 1 at d0 (r<=3, g<=4)", ok_def)
print("   r=1 'missing 32nd region' moments: (g, d0, full, shadow):", first_table)

# quasi-polynomial: along each residue class mod g+1, Delta^{2r+3} vanishes eventually
ok_qp = True
for r in range(1, 4):
    for g in range(0, 4):
        vals = [D_rg(r, g, d) for d in range(N)]
        for res in range(g + 1):
            sub = [vals[d] for d in range(res, N, g + 1)]
            sub = sub[6:]                      # past the pre-period
            for _ in range(2 * r + 3):
                sub = [sub[i + 1] - sub[i] for i in range(len(sub) - 1)]
            if any(v != 0 for v in sub): ok_qp = False
check("D_{r,g} is a quasi-polynomial of degree 2r+2, period g+1 (finite differences)", ok_qp)

# ---------- new sequences (OEIS candidates) ----------
print()
print("NEW-SEQUENCE CANDIDATES (first 24 terms):")
for (r, g, label) in [(1, 2, "D_{1,2} (gap-2 polygonal vs Narayana A000930)"),
                      (1, 3, "D_{1,3} (gap-3 polygonal vs A003269)"),
                      (1, 4, "D_{1,4} (gap-4 polygonal vs A003520)"),
                      (2, 1, "D_{2,1} (rank-2 'polyhedral-gonal' Fibonacci shadow, sextic)"),
                      (2, 2, "D_{2,2} (rank-2 gap-2 shadow)")]:
    print(f"  {label}:")
    print("   ", [D_rg(r, g, d) for d in range(24)])
t1 = [full_g(1, 30)[d] - D_rg(1, 1, d) for d in range(30)]
print("  deficit t_{1,1}(d) = F(d+1) - D(d):", t1[8:24])

# Brown completeness of the new families (terms sorted ascending; they are nondecreasing)
def brown(a):
    s = 0
    for i, t in enumerate(a):
        if i == 0 and t != 1: return False
        if t > 1 + s: return False
        s += t
    return True
check("Brown criterion: D_{1,g} complete for g = 0..4; D_{2,1}, D_{2,2} complete",
      all(brown([D_rg(1, g, d) for d in range(48)]) for g in range(5)) and
      all(brown([D_rg(2, g, d) for d in range(48)]) for g in (1, 2)))

# ---------- integration ladder (pyramidal etc.) ----------
# s-fold integrated polygonal columns: I_s(k, j) = sum over j' <= j, s times, of P(k+1, .)
# closed: I_s(k, j) = C(j+s-1... verified: I_s(k,j) = C(j+s-1, s+1)?? -- computed directly:
def I_s(s, k, j):     # s-fold partial sum of column k of the polygonal table (col 0 stays 1s)
    if k == 0: return 1
    if s == 0: return comb(j, 1) + (k - 1) * comb(j, 2)
    return sum(I_s(s - 1, k, t) for t in range(1, j + 1))
ok_int = True
rows_int = {}
for s in range(0, 3):
    rows = [sum(I_s(s, k, n - k + 1) for k in range(n + 1)) for n in range(20)]
    rows_int[s] = rows
    pred = [1 + comb(n + s + 1, s + 2) + comb(n + s + 1, s + 4) for n in range(20)]
    # column 0 of the integrated table stays all-ones, so adjust: predicted with 1 replaced?
    if rows != pred: ok_int = False
check("INTEGRATION LADDER: s-fold integrated polygonal rows = 1 + C(n+s+1, s+2) + C(n+s+1, s+4)",
      ok_int)
print("   s=1 (pyramidal Moser):", rows_int[1][:12])
print("   s=2 (4D hyper-pyramidal Moser):", rows_int[2][:12])

# ---------- tournament side: x = (n^3-n)/3 - 8 c3 ; k(T) = c3(T) ----------
def scores(T, n): return [sum(T[u][v] for v in range(n) if v != u) for u in range(n)]
def c3_direct(T, n):
    # a triple is cyclic iff every vertex has out-degree 1 within it (no source)
    c = 0
    for a, b, c_ in itertools.combinations(range(n), 3):
        outd = {a: 0, b: 0, c_: 0}
        for u, v in ((a, b), (a, c_), (b, c_)):
            outd[u if T[u][v] else v] += 1
        if sorted(outd.values()) == [1, 1, 1]: c += 1
    return c
def x_of(T, n):
    return sum((2 * s - (n - 1)) ** 2 for s in scores(T, n))

ok_id = True
n = 5
for bits in range(2 ** 10):
    T = [[0] * n for _ in range(n)]
    idx = 0
    for u in range(n):
        for v in range(u + 1, n):
            T[u][v] = (bits >> idx) & 1
            T[v][u] = 1 - T[u][v]
            idx += 1
    if x_of(T, n) != (n ** 3 - n) // 3 - 8 * c3_direct(T, n): ok_id = False
check("x = (n^3-n)/3 - 8*c3 (EXHAUSTIVE n=5, all 1024)", ok_id)

rng = random.Random(313)
ok_id2 = True
for n in (6, 7, 8, 9):
    for _ in range(40):
        T = [[0] * n for _ in range(n)]
        for u in range(n):
            for v in range(u + 1, n):
                T[u][v] = rng.randint(0, 1); T[v][u] = 1 - T[u][v]
        if x_of(T, n) != (n ** 3 - n) // 3 - 8 * c3_direct(T, n): ok_id2 = False
check("x = (n^3-n)/3 - 8*c3 (random n=6..9)", ok_id2)

# tie-splitting walk: each step Delta x = +8 AND Delta c3 = -1; halts transitive in c3 steps
ok_walk = True
for n in (5, 6, 7, 8):
    for _ in range(25):
        T = [[0] * n for _ in range(n)]
        for u in range(n):
            for v in range(u + 1, n):
                T[u][v] = rng.randint(0, 1); T[v][u] = 1 - T[u][v]
        k0 = c3_direct(T, n); steps = 0
        while True:
            s = scores(T, n)
            tied = [(u, v) for u in range(n) for v in range(u + 1, n) if s[u] == s[v]]
            if not tied: break
            u, v = tied[0]
            x0, c0 = x_of(T, n), c3_direct(T, n)
            T[u][v], T[v][u] = T[v][u], T[u][v]
            if x_of(T, n) - x0 != 8 or c3_direct(T, n) - c0 != -1: ok_walk = False
            steps += 1
        if not (steps == k0 and c3_direct(T, n) == 0 and sorted(scores(T, n)) == list(range(n))):
            ok_walk = False
check("tie-split walk: every step Dx=+8 AND Dc3=-1; halts transitive after exactly c3(T) steps",
      ok_walk)

# |Aut| odd, exhaustively n<=5; C5{1,2} floor: c3=5, |Aut|=5, unique max
ok_aut, maxaut, floor_hits = True, {}, []
for n in (3, 4, 5):
    m = n * (n - 1) // 2
    for bits in range(2 ** m):
        T = [[0] * n for _ in range(n)]
        idx = 0
        for u in range(n):
            for v in range(u + 1, n):
                T[u][v] = (bits >> idx) & 1; T[v][u] = 1 - T[u][v]; idx += 1
        a = 0
        for perm in itertools.permutations(range(n)):
            if all(T[perm[u]][perm[v]] == T[u][v] for u in range(n) for v in range(n) if u != v):
                a += 1
        if a % 2 == 0: ok_aut = False
        maxaut[n] = max(maxaut.get(n, 0), a)
        if n == 5 and x_of(T, n) == 0:
            floor_hits.append((a, c3_direct(T, n)))
check("|Aut(T)| odd for ALL tournaments n<=5 (exhaustive; => solvable by Feit-Thompson)", ok_aut)
check("n=5 floor x=0: every floor tournament has c3=5 (=k_max) and the max |Aut|=5 occurs there",
      all(c == 5 for (_, c) in floor_hits) and maxaut[5] == 5 and max(a for a, _ in floor_hits) == 5)
print("   max |Aut| by n:", maxaut, "| n=5 floor count:", len(floor_hits),
      "(labelled copies of the round tournament C5)")

# CD-level table
print()
print("CAYLEY-DICKSON LEVELS n = 2^k+1: ceiling (n^3-n)/3 and k_max = c3_max:")
for kk in range(1, 6):
    n_ = 2 ** kk + 1
    print(f"   n={n_:3d}: ceiling={(n_**3 - n_)//3:6d}   c3_max={(n_**3 - n_)//24:6d}")
print("   (n=9: ceiling 240 = #E8 roots; c3_max 30 = #icosahedron edges — numerology, flagged)")

print()
fails = [nm for nm, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed, {len(fails)} failed ===")
for f in fails: print("FAILED:", f)
