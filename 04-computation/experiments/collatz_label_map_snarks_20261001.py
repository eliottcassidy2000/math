#!/usr/bin/env python3
"""The owner's label map F, a Paley-coded restatement of Collatz, and snarks (opus S15, twentieth note, 2026-10-01).

Checks:
  A. the corrected map F(2N) = 3N, F(2N-1) = 2N-1-K_N (K = 0,2,1,4,2,10,3,9,4,15,5) is the Syracuse map
     S(A) = oddpart(3A+1) in labels M = (A+1)/2; closed form F(M) = (oddpart(3M-1)+1)/2; three-rule
     recursion F(2N) = 3N, F(4j+1) = 3j+1, F(4n-1) = F(n) (the microcosm); K_{2j+1} = j, K_{4m} = 5m-1,
     K_{4m-2} = K_m + 6m - 4; the identity 6K(A)+1 = (2^v - 3) S(A); the multiplicity of every value m in K
     is #{v >= 2 : (2^v - 3) | 6m+1}; mean multiplicity sum_{v>=2} 1/(2^v-3); "two copies + one 0" holds for
     ascents + shallow descents (the units 2^1-3 = -1, 2^2-3 = 1)
  B. Paley restatement: Collatz <=> 7*Phi_T(n) is an integer for every n >= 1, and then 7*Phi_T(n) = -2^sigma(n)
     (mod 7) lies in NQR_7 = {3,5,6} (sigma = number of T-steps to 1); negative integers have period-2, 5, 18
     tails; the q-family x/2, (qx+1)/2: the Paley codes are the rational cycles at q = 1; every integer cycle
     of T1 with period L <= 16 is one of the five known ones; the size bound that settles q = 1 and where it
     fails at q = 3 (near-resonances 2^L ~ 3^k)
  C. snarks: Fano colourings = nowhere-zero Z_2^3-flows; Fano lines = translates of {1,2,4} in Singer
     coordinates; any colouring using one line, two lines, concurrent lines or missing a point is a
     3-edge-colouring; Petersen: all nowhere-zero Z_2^3-flows enumerated (every one uses all 7 points and a
     pencil-plus-line of 4 lines at least); the flower snark J5 and the two 18-vertex dot products P.P (the
     Blanusa snarks) need exactly the 4-line pencil-plus-line type; functional graphs are 3-edge-colourable;
     residue classes mod 2^k with no coefficient descent never vanish (no finite reducible set)
Reproduce: python 04-computation/experiments/collatz_label_map_snarks_20261001.py   (about 1-2 minutes)
"""
from fractions import Fraction
from itertools import product, combinations
from collections import Counter
import math

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x


def v2(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v


def S(A):  # Syracuse map on odd A
    return oddpart(3 * A + 1)


def F(M):  # owner's map in labels M = (A+1)/2
    return (S(2 * M - 1) + 1) // 2


def Kof(A):  # descent of an odd A, (A - S(A))/2  (negative on the ascending branch)
    return (A - S(A)) // 2


def KN(N):  # owner's K_N, indexed by A = 4N - 3
    return Kof(4 * N - 3)


# ---------------------------------------------------------------- A. the label map F
print("=== A. the owner's map F is the Syracuse map in labels M = (A+1)/2 ===")
K_given = [0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5]
check(all(KN(N) == K_given[N - 1] for N in range(1, 12)),
      "K_N = (A - S(A))/2 at A = 4N-3 reproduces the owner's K = 0,2,1,4,2,10,3,9,4,15,5")
check(all(F(2 * N - 1) == 2 * N - 1 - K_given[N - 1] and F(2 * N) == 3 * N for N in range(1, 12)),
      "F(2N) = 3N and F(2N-1) = 2N-1-K_N for N = 1..11 (owner's corrected rule)")
check([S(A) for A in (1, 5, 9, 13)] == [1, 1, 7, 5], "owner's examples: S(1)=1, S(5)=1, S(9)=7, S(13)=5")
X = 10 ** 6
check(all(F(M) == (oddpart(3 * M - 1) + 1) // 2 for M in range(1, X + 1)),
      f"closed form F(M) = (oddpart(3M-1)+1)/2 for M <= {X}")
Frec = [0] * (X + 1)
for M in range(1, X + 1):
    if M % 2 == 0:
        Frec[M] = 3 * M // 2
    elif M % 4 == 1:
        Frec[M] = 3 * (M - 1) // 4 + 1
    else:
        Frec[M] = Frec[(M + 1) // 4]
check(all(Frec[M] == F(M) for M in range(1, X + 1)),
      f"three-rule recursion F(2N)=3N, F(4j+1)=3j+1, F(4n-1)=F(n) defines F for M <= {X}")
check(all(KN(2 * j + 1) == j for j in range(0, X // 2)) and all(KN(4 * m) == 5 * m - 1 for m in range(1, X // 4))
      and all(KN(4 * m - 2) == KN(m) + 6 * m - 4 for m in range(1, X // 4)),
      f"K_(2j+1) = j, K_(4m) = 5m-1, K_(4m-2) = K_m + 6m - 4 (N <= {X})")
check(all(6 * Kof(A) + 1 == (2 ** v2(3 * A + 1) - 3) * S(A) for A in range(1, 4 * X, 2)),
      f"identity 6K(A) + 1 = (2^v - 3) S(A), v = v2(3A+1), for all odd A < {4 * X}")

# multiplicities of K
Xm = 100000
cnt = Counter()
for N in range(1, (8 * Xm + 1 + 3) // 4 + 2):
    k = KN(N)
    if k <= Xm:
        cnt[k] += 1
D = [2 ** v - 3 for v in range(2, 64) if 2 ** v - 3 <= 6 * Xm + 1]


def mult_formula(m):
    return sum(1 for d in D if (6 * m + 1) % d == 0)


check(all(cnt[m] == mult_formula(m) for m in range(0, Xm + 1)),
      f"multiplicity of m in K = #{{v >= 2 : (2^v - 3) | 6m+1}} for every 0 <= m <= {Xm}")
dist = Counter(cnt[m] for m in range(0, Xm + 1))
mean_emp = sum(cnt[m] for m in range(0, Xm + 1)) / (Xm + 1)
mean_th = sum(Fraction(1, 2 ** v - 3) for v in range(2, 80))
print(f"   multiplicity distribution on [0,{Xm}]: {dict(sorted(dist.items()))}; mean {mean_emp:.5f}; "
      f"sum_(v>=2) 1/(2^v-3) = {float(mean_th):.5f}")
check(cnt[0] == 1 and cnt[1] == 1 and cnt[2] == 2 and cnt[24] == 3,
      "0 and 1 occur once, 2 twice (13 = 2^4-3), 24 three times (145 = 5*29): 'two copies of each' is false in K")
check(abs(mean_emp - float(mean_th)) < 0.01, "mean multiplicity of K ~ 1.3437, not 2")
# the 'two copies + one 0' that is true: ascents (v=1, unit -1) + shallow descents (v=2, unit +1)
Xd = 50000
absd = Counter()
for M in range(1, 4 * Xd + 1):
    if M % 4 != 3:
        absd[abs(F(M) - M)] += 1
check(absd[0] == 1 and all(absd[m] == 2 for m in range(1, Xd)),
      f"over labels M != 3 (mod 4), M <= {4 * Xd}: |F(M)-M| takes 0 once and every 1 <= m < {Xd} exactly twice "
      "(one ascent at M=2m, one shallow descent at M=4m+1)")
check(all(F(4 * n - 1) - (4 * n - 1) == (F(n) - n) - 3 * n + 1 for n in range(1, X // 4)),
      "microcosm: Delta(4n-1) = Delta(n) - 3n + 1, Delta(M) = F(M) - M (the quarter M = 3 mod 4 is a rescaled copy)")
check(all(2 ** v - 3 not in (1, -1) for v in range(3, 64)),
      "the two 'copies' are the only units 2^v - 3 = -1, 1 (v = 1, 2): Gersonides again")
print("   2^K_N for N = 1..8:", [2 ** KN(N) for N in range(1, 9)], "(1,4,2 then 16,4,1024: no pattern; NUMEROLOGY)")

# ---------------------------------------------------------------- B. Paley restatement and the q = 1 shadow
print("=== B. Paley restatement; the q-family; where the q = 1 argument fails at q = 3 ===")


def T(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


ok = True
res_count = Counter()
inv7 = None
for n in range(1, 20001):
    bits = []
    x = n
    while x != 1:
        bits.append(x & 1)
        x = T(x)
    sigma = len(bits)
    prefix = sum(b << j for j, b in enumerate(bits))
    num = 7 * prefix - 2 ** sigma          # predicted 7 * Phi_T(n)
    Kb = sigma + 40
    mod = 1 << Kb
    # actual first Kb bits of Phi_T(n): continue the orbit through the trivial cycle
    y, phi = n, 0
    for j in range(Kb):
        phi |= (y & 1) << j
        y = T(y)
    if (7 * phi - num) % mod != 0:
        ok = False
    r = num % 7
    res_count[r] += 1
    if r != (-pow(2, sigma, 7)) % 7:
        ok = False
check(ok, "for n <= 20000: 7*Phi_T(n) = 7*prefix - 2^sigma(n) (checked to 40 bits past the entry into 1)")
check(set(res_count) <= {3, 5, 6}, f"residues 7*Phi_T(n) mod 7 are non-residues: {dict(sorted(res_count.items()))}")
# negative integers: tails of period 2, 5, 18
tails = {}
for n in range(-1, -3001, -1):
    seen = {}
    x = n
    while x not in seen:
        seen[x] = 1
        x = T(x)
    cyc = [x]
    y = T(x)
    while y != x:
        cyc.append(y)
        y = T(y)
    tails[min(cyc, key=abs)] = len(cyc)
check(tails == {-1: 2, -5: 5, -17: 18}, f"negative orbits end in T-cycles of lengths {tails} -> denominators 3, 31, 2^18-1")


def Tq(q, x):
    return x // 2 if x % 2 == 0 else (q * x + 1) // 2


cycles_by_q = {}
for q in (1, 3, 5, 7):
    cyc_set = set()
    for x0 in range(-300, 301):
        x, seen = x0, {}
        for _ in range(2000):
            if x in seen or abs(x) > 10 ** 12:
                break
            seen[x] = 1
            x = Tq(q, x)
        if x in seen:
            c = [x]
            y = Tq(q, x)
            while y != x:
                c.append(y)
                y = Tq(q, y)
            cyc_set.add(min(c, key=abs))
    cycles_by_q[q] = sorted(cyc_set)
print("   integer cycles (element of least |x|) of x/2,(qx+1)/2 from |x0| <= 300:", cycles_by_q)
check(cycles_by_q[1] == [0, 1] and cycles_by_q[3] == [-17, -5, -1, 0, 1],
      "q = 1: only 0 and 1 (Collatz trivially true); q = 3: 1 and the three negative cycles")


def U(x):  # q = 1 map on Z_(2): x/2 or (x+1)/2
    return x / 2 if x.numerator % 2 == 0 else (x + 1) / 2


orb1 = [Fraction(1, 7)]
for _ in range(3):
    orb1.append(U(orb1[-1]))
orb2 = [Fraction(3, 7)]
for _ in range(3):
    orb2.append(U(orb2[-1]))
check(orb1 == [Fraction(1, 7), Fraction(4, 7), Fraction(2, 7), Fraction(1, 7)]
      and orb2 == [Fraction(3, 7), Fraction(5, 7), Fraction(6, 7), Fraction(3, 7)],
      "QR_7/7 and NQR_7/7 are the two rational 3-cycles of the q = 1 map (words 100 and 110)")


def cw(word, q):
    k = sum(word)
    c, ones_after = 0, k
    for t, b in enumerate(word):
        if b:
            ones_after -= 1
            c += q ** ones_after * 2 ** t
    return c


found = set()
for L in range(1, 17):
    for word in product((0, 1), repeat=L):
        k = sum(word)
        den = 2 ** L - 3 ** k
        c = cw(word, 3)
        if c % den == 0:
            x = c // den
            orbit = [x]
            for _ in range(L - 1):
                orbit.append(Tq(3, orbit[-1]))
            found.add(min(orbit, key=abs))
check(found == {0, 1, -1, -5, -17}, f"every integer T1-cycle with period L <= 16 (all 2^L words): least-|x| elements {sorted(found)}")
check(all(cw(w, 1) < 2 ** len(w) - 1 for L in range(1, 15) for w in product((0, 1), repeat=L) if 0 < sum(w) < L),
      "q = 1 size bound: c_w(1) < 2^L - 1 for every non-constant word (no integer cycle but 0, 1)")
rows = []
for L in range(1, 40):
    k = int(L * math.log(2) / math.log(3))
    while 2 ** L <= 3 ** k:
        k -= 1
    lb = Fraction(3 ** k - 2 ** k, 2 ** L - 3 ** k)
    rows.append((L, k, 2 ** L - 3 ** k, float(lb)))
print("   q = 3, positive sheet: (L, k, 2^L - 3^k, min_w c_w(3)/(2^L-3^k)) with ratio > 1 (size bound fails):")
print("   ", [(L, k, d, round(r, 2)) for (L, k, d, r) in rows if r > 1][:12])

# ---------------------------------------------------------------- C. snarks and Fano colourings
print("=== C. snarks: Fano colourings, the Paley line {1,2,4}, and the translations snarks need ===")
# F_8 = F_2[a]/(a^3+a+1), points a^i (i in Z_7) as 3-bit vectors
pw = [1]
for i in range(1, 7):
    x = pw[-1] << 1
    if x & 8:
        x ^= 0b1011
    pw.append(x)
check(len(set(pw)) == 7, "a^0..a^6 are the 7 nonzero vectors of F_2^3 (Singer cycle)")
D = (1, 2, 4)
lines = [frozenset(pw[(d + t) % 7] for d in D) for t in range(7)]
check(all(pw[(1 + t) % 7] ^ pw[(2 + t) % 7] ^ pw[(4 + t) % 7] == 0 for t in range(7)) and len(set(lines)) == 7,
      "the translates D+t of the Paley set D = {1,2,4} are the 7 Fano lines (a^(1+t)+a^(2+t)+a^(4+t) = 0)")
check(all(frozenset(pw[(2 * d) % 7] for d in D) == lines[0] for _ in [0]),
      "the Frobenius x -> x^2 (exponent x2) fixes the line D and rotates its points 1 -> 2 -> 4")


def line_type(Ls):
    Ls = list(Ls)
    l = len(Ls)
    pts = Counter(p for L_ in Ls for p in L_)
    allpts = len(pts) == 7
    if l == 3:
        return "3-concurrent" if max(pts.values()) == 3 else "3-triangle"
    if l == 4:
        return "4-pencil+line" if max(pts.values()) == 3 else "4-quadrilateral"
    return f"{l}"


def petersen():
    E = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    return 10, E


def flower(n):
    a = lambda i: 4 * (i % n)
    b = lambda i: 4 * (i % n) + 1
    c = lambda i: 4 * (i % n) + 2
    d = lambda i: 4 * (i % n) + 3
    E = []
    for i in range(n):
        E += [(a(i), b(i)), (a(i), c(i)), (a(i), d(i)), (b(i), b(i + 1))]
    cyc = [c(i) for i in range(n)] + [d(i) for i in range(n)]
    for i in range(2 * n):
        E.append((cyc[i], cyc[(i + 1) % (2 * n)]))
    return 4 * n, E


def dot_product(G, H, e1, e2, hx):
    """Isaacs dot product: delete independent edges e1=ab, e2=cd of G; delete adjacent x,y of H (edge hx);
    join a,b to x's other neighbours and c,d to y's other neighbours."""
    nG, EG = G
    nH, EH = H
    a, b = e1
    c, d = e2
    x, y = hx
    E = [e for e in EG if set(e) not in (set(e1), set(e2))]
    nbr = {v: [] for v in range(nH)}
    for u, w in EH:
        nbr[u].append(w)
        nbr[w].append(u)
    keep = [v for v in range(nH) if v not in (x, y)]
    idx = {v: nG + i for i, v in enumerate(keep)}
    for u, w in EH:
        if x not in (u, w) and y not in (u, w):
            E.append((idx[u], idx[w]))
    x1, x2 = [w for w in nbr[x] if w != y]
    y1, y2 = [w for w in nbr[y] if w != x]
    E += [(a, idx[x1]), (b, idx[x2]), (c, idx[y1]), (d, idx[y2])]
    return nG + nH - 2, E


def girth(n, E):
    adj = {v: [] for v in range(n)}
    for u, w in E:
        adj[u].append(w)
        adj[w].append(u)
    best = 10 ** 9
    for s in range(n):
        dist = {s: 0}
        par = {s: -1}
        q = [s]
        for u in q:
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    par[w] = u
                    q.append(w)
                elif par[u] != w:
                    best = min(best, dist[u] + dist[w] + 1)
    return best


def fano_colourable(n, E, allowed):
    """Search with unit propagation: label edges by points (3-bit vectors) so that every vertex star is a
    line in `allowed` (two labelled edges at a vertex force the third, their XOR)."""
    from itertools import permutations
    inc = [[] for _ in range(n)]
    for i, (u, w) in enumerate(E):
        inc[u].append(i)
        inc[w].append(i)
    allowed = [frozenset(L_) for L_ in allowed]
    aset = set(allowed)
    lab = [0] * len(E)

    def propagate(trail):
        changed = True
        while changed:
            changed = False
            for v in range(n):
                vals = [lab[i] for i in inc[v]]
                z = vals.count(0)
                if z == 0:
                    if frozenset(vals) not in aset:
                        return False
                elif z == 1:
                    a, b = [x for x in vals if x]
                    c = a ^ b
                    if c == 0 or frozenset((a, b, c)) not in aset:
                        return False
                    i = inc[v][vals.index(0)]
                    lab[i] = c
                    trail.append(i)
                    changed = True
                elif z == 2:
                    a = next(x for x in vals if x)
                    if not any(a in L_ for L_ in allowed):
                        return False
        return True

    def rec():
        best, bz = None, 4
        for v in range(n):
            z = sum(1 for i in inc[v] if lab[i] == 0)
            if 0 < z < bz:
                best, bz = v, z
        if best is None:
            return True
        es = [i for i in inc[best] if lab[i] == 0]
        fixed = {lab[i] for i in inc[best] if lab[i]}
        for L_ in allowed:
            if not fixed <= L_:
                continue
            rest = sorted(L_ - fixed)
            if len(rest) != len(es):
                continue
            for perm in permutations(rest):
                trail = []
                for i, x in zip(es, perm):
                    lab[i] = x
                    trail.append(i)
                if propagate(trail) and rec():
                    return True
                for i in trail:
                    lab[i] = 0
        return False
    return rec()


L0 = lines[0]
type_reps = {
    "1": [lines[0]],
    "2": [lines[0], lines[1]],
}
# concurrent triple: the 3 lines through a point p; triangle: 3 lines with no common point
p0 = pw[0]
pencil = [L_ for L_ in lines if p0 in L_]
nonp = [L_ for L_ in lines if p0 not in L_]
type_reps["3-concurrent"] = pencil
tri = next(T3 for T3 in combinations(lines, 3) if line_type(T3) == "3-triangle")
type_reps["3-triangle"] = list(tri)
type_reps["4-quadrilateral"] = nonp
type_reps["4-pencil+line"] = pencil + [nonp[0]]
type_reps["5"] = lines[:5]
type_reps["7"] = lines
check(line_type(type_reps["4-pencil+line"]) == "4-pencil+line" and line_type(nonp) == "4-quadrilateral",
      "line-set types: 4 lines = pencil+line or quadrilateral (complement of a pencil)")

P10 = petersen()
snarks = {"Petersen": P10, "flower J5": flower(5)}   # J7 (28 vertices) is too slow for this pure-Python search
# Blanusa snarks: dot products P.P with independent edges at distance 1 and 2 in P
nP, EP = P10
adjP = {v: set() for v in range(nP)}
for u, w in EP:
    adjP[u].add(w)
    adjP[w].add(u)
e1 = EP[0]  # (0,1)
cand = [e for e in EP if not set(e) & set(e1)]
d1 = next(e for e in cand if any(w in adjP[u] for u in e1 for w in e))
d2 = next(e for e in cand if not any(w in adjP[u] for u in e1 for w in e))
snarks["Blanusa (dot, dist 1)"] = dot_product(P10, P10, e1, d1, EP[0])
snarks["Blanusa (dot, dist 2)"] = dot_product(P10, P10, e1, d2, EP[0])
for name, (n, E) in snarks.items():
    degs = Counter()
    for u, w in E:
        degs[u] += 1
        degs[w] += 1
    cubic = all(degs[v] == 3 for v in range(n)) and len(E) == 3 * n // 2
    g = girth(n, E)
    res = {t: fano_colourable(n, E, Ls) for t, Ls in type_reps.items()}
    print(f"   {name}: n={n}, cubic={cubic}, girth={g}, colourable by line-set type: {res}")
    check(cubic and g >= 5 and not res["1"] and not res["2"] and not res["3-concurrent"] and not res["3-triangle"]
          and not res["4-quadrilateral"] and res["4-pencil+line"],
          f"{name}: snark (not 3-edge-colourable), needs >= 4 lines, and the 4-line pencil+line type suffices")

# Petersen: all nowhere-zero Z_2^3 flows via the cycle space
n, E = P10
m = len(E)
# spanning tree + fundamental cycles
adj = {v: [] for v in range(n)}
for i, (u, w) in enumerate(E):
    adj[u].append((w, i))
    adj[w].append((u, i))
par = {0: (None, None)}
q = [0]
for u in q:
    for w, i in adj[u]:
        if w not in par:
            par[w] = (u, i)
            q.append(w)
tree = {par[v][1] for v in par if par[v][1] is not None}


def path_edges(u):
    out = []
    while par[u][0] is not None:
        out.append(par[u][1])
        u = par[u][0]
    return out


basis = []
for i, (u, w) in enumerate(E):
    if i in tree:
        continue
    pu, pw_ = set(path_edges(u)), set(path_edges(w))
    basis.append((pu ^ pw_) | {i})
check(len(basis) == m - n + 1 == 6, "Petersen cycle space has dimension 6 (8^6 = 262144 Z_2^3-flows)")
flows = 0
types = Counter()
allseven = True
for coeffs in product(range(8), repeat=6):
    f = [0] * m
    for c, cyc in zip(coeffs, basis):
        if c:
            for i in cyc:
                f[i] ^= c
    if 0 in f:
        continue
    flows += 1
    stars = set()
    for v in range(n):
        es = [i for i, e in enumerate(E) if v in e]
        stars.add(frozenset(f[i] for i in es))
    if len(set(f)) != 7:
        allseven = False
    types[line_type(stars)] += 1
print(f"   Petersen: {flows} nowhere-zero Z_2^3-flows (Fano colourings); line-set types used: {dict(types)}")
check(flows > 0 and allseven and set(types) <= {"4-pencil+line", "5", "6", "7"},
      "every Fano colouring of Petersen uses all 7 points and >= 4 lines; the 4-line ones are pencil+line")
check(flows % 168 == 0, "the flow count is a multiple of |GL(3,2)| = 168 (GL(3,2) acts freely on colourings)")


def functional_graph_class1(N):
    """The Collatz T-graph on [1, N] closed under T: components are trees or unicyclic, max degree 3."""
    verts = set()
    for n0 in range(1, N + 1):
        x = n0
        while x not in verts:
            verts.add(x)
            x = T(x)
    deg = Counter()
    for x in verts:
        deg[x] += 1
        deg[T(x)] += 1
    return max(deg.values()), len(verts)


mx, nv = functional_graph_class1(2000)
check(mx <= 3, f"the Collatz graph on the {nv} vertices reached from [1,2000] has max degree {mx} (functional graph: "
      "each component a tree or unicyclic, so 3-edge-colourable; every non-cycle edge is a bridge)")

# reducibility: residue classes mod 2^k with no coefficient descent within k steps never vanish
dens = []
for k in (5, 10, 20, 40, 80):
    dp = {0: 1}  # number of odd steps so far -> count, restricted to paths with 3^o >= 2^j for all j
    for j in range(1, k + 1):
        nd = Counter()
        for o, cnum in dp.items():
            for b in (0, 1):
                o2 = o + b
                if 3 ** o2 >= 2 ** j:
                    nd[o2] += cnum
        dp = nd
    dens.append((k, sum(dp.values()) / 2 ** k))
print("   density of residues mod 2^k with no coefficient descent in k steps:", [(k, f"{d:.3e}") for k, d in dens])
check(all(d > 0 for _, d in dens), "never zero (the class of -1 mod 2^k ascends for k steps): no finite reducible set")

print()
print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
