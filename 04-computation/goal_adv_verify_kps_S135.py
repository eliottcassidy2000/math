# Adversarial INDEPENDENT verification of the GOAL-theorem proof (H = |Aut| => C3).
# Written from scratch; shares no code with .scratch/goal_thm_crossover_dp_kps_S134.py.
import sys, math
from functools import lru_cache
sys.setrecursionlimit(100000)

def fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r

def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x

def v3(x):
    c = 0
    while x % 3 == 0:
        x //= 3
        c += 1
    return c

# ---------------- tournament utilities (adj[u] = bitmask of out-neighbours) -------------
def is_strong(adj, n):
    full = (1 << n) - 1
    radj = [0] * n
    for u in range(n):
        m = adj[u]
        while m:
            b = m & -m
            v = b.bit_length() - 1
            radj[v] |= 1 << u
            m ^= b
    def reach(ad):
        seen = 1
        frontier = 1
        while True:
            nxt = 0
            m = frontier
            while m:
                b = m & -m
                u = b.bit_length() - 1
                nxt |= ad[u]
                m ^= b
            frontier = nxt & ~seen
            if not frontier:
                return seen
            seen |= frontier
    return reach(adj) == full and reach(radj) == full

def ham_count(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, 1 << n):
        row = dp[s]
        for last in range(n):
            c = row[last]
            if not c:
                continue
            av = adj[last] & ~s
            while av:
                b = av & -av
                v = b.bit_length() - 1
                dp[s | b][v] += c
                av ^= b
    full = (1 << n) - 1
    return sum(dp[full])

def aut_count(adj, n):
    score = [bin(adj[u]).count('1') for u in range(n)]
    cand = [[w for w in range(n) if score[w] == score[i]] for i in range(n)]
    mapping = [-1] * n
    cnt = 0
    def rec(i, used):
        nonlocal cnt
        if i == n:
            cnt += 1
            return
        for w in cand[i]:
            if used >> w & 1:
                continue
            ok = True
            for j in range(i):
                if (adj[i] >> j & 1) != (adj[w] >> mapping[j] & 1):
                    ok = False
                    break
            if ok:
                mapping[i] = w
                rec(i + 1, used | (1 << w))
        mapping[i] = -1
    rec(0, 0)
    return cnt

def find_iso(adjA, adjB, n):
    sA = sorted(bin(adjA[u]).count('1') for u in range(n))
    sB = sorted(bin(adjB[u]).count('1') for u in range(n))
    if sA != sB:
        return None
    scoreA = [bin(adjA[u]).count('1') for u in range(n)]
    scoreB = [bin(adjB[u]).count('1') for u in range(n)]
    cand = [[w for w in range(n) if scoreB[w] == scoreA[i]] for i in range(n)]
    mapping = [-1] * n
    def rec(i, used):
        if i == n:
            return list(mapping)
        for w in cand[i]:
            if used >> w & 1:
                continue
            ok = True
            for j in range(i):
                if (adjA[i] >> j & 1) != (adjB[w] >> mapping[j] & 1):
                    ok = False
                    break
            if ok:
                mapping[i] = w
                r = rec(i + 1, used | (1 << w))
                if r:
                    return r
        mapping[i] = -1
        return None
    return rec(0, 0)

def apply_perm(adj, n, p):
    out = [0] * n
    for u in range(n):
        m = adj[u]
        while m:
            b = m & -m
            v = b.bit_length() - 1
            out[p[u]] |= 1 << p[v]
            m ^= b
    return out

# ---------------- group machinery ----------------
def close_group(gens, n):
    idp = tuple(range(n))
    seen = {idp}
    stack = [idp]
    while stack:
        g = stack.pop()
        for h in gens:
            f = tuple(h[g[i]] for i in range(n))
            if f not in seen:
                seen.add(f)
                stack.append(f)
    return seen

def pair_orbits(gens, n):
    seen = {}
    orbits = []
    for a in range(n):
        for b in range(n):
            if a == b or (a, b) in seen:
                continue
            comp = [(a, b)]
            seen[(a, b)] = len(orbits)
            k = 0
            while k < len(comp):
                x, y = comp[k]
                k += 1
                for g in gens:
                    z = (g[x], g[y])
                    if z not in seen:
                        seen[z] = len(orbits)
                        comp.append(z)
            orbits.append(comp)
    return orbits, seen

def invariant_tournaments(gens, n):
    orbits, oid = pair_orbits(gens, n)
    # pair each orbit with its reversal orbit
    rev = {}
    for i, comp in enumerate(orbits):
        a, b = comp[0]
        j = oid[(b, a)]
        assert j != i, "self-paired orbit => even-order element, impossible for odd group"
        rev[i] = j
    classes = []
    done = set()
    for i in range(len(orbits)):
        if i in done:
            continue
        j = rev[i]
        classes.append((i, j))
        done.add(i)
        done.add(j)
    ts = []
    bitschoices = []
    for mask in range(1 << len(classes)):
        adj = [0] * n
        for c, (i, j) in enumerate(classes):
            chosen = orbits[i] if (mask >> c) & 1 == 0 else orbits[j]
            for (x, y) in chosen:
                adj[x] |= 1 << y
        ts.append(adj)
        bitschoices.append(mask)
    return ts, classes, orbits

# =========================================================================
print("=" * 78)
print("SECTION 1: m = 9  (window [f(9), o(9)] = [75, 81])")
print("=" * 78)

# T9 = C3[C3,C3,C3]
def build_T9():
    adj = [0] * 9
    for g in range(3):
        base = 3 * g
        for t in range(3):
            adj[base + t] |= 1 << (base + (t + 1) % 3)   # inner C3
    for g in range(3):
        for x in range(3):
            for y in range(3):
                adj[3 * g + x] |= 1 << (3 * ((g + 1) % 3) + y)  # block g -> g+1
    return adj

T9 = build_T9()
assert is_strong(T9, 9), "T9 not strong?!"
A9 = aut_count(T9, 9)
H9 = ham_count(T9, 9)
print(f"T9 strong = True, |Aut(T9)| = {A9}, H(T9) = {H9}, tc = {H9 // A9 if H9 % A9 == 0 else 'NOT DIVISIBLE'}")
assert A9 == 81 and H9 == 3159 and H9 % A9 == 0 and H9 // A9 == 39

# Lagrange window
op9 = oddpart(fact(9))
print(f"oddpart(9!) = {op9} (= 3^4*5*7: {op9 == 3**4*5*7}); v3(9!) = {v3(fact(9))}")
win = [d for d in range(75, 82) if d % 2 == 1 and op9 % d == 0]
print(f"odd divisors of oddpart(9!) in [75,81]: {win}")
assert win == [81]
assert oddpart(fact(8)) == 315 and 315 % 27 != 0, "27 | oddpart(8!)??"
print(f"oddpart(8!) = {oddpart(fact(8))}; 27 divides it: {oddpart(fact(8)) % 27 == 0}  (so |Aut|=27 impossible at m=8; o(8) vs f(8)=45)")

# Sylow-3 of S9, standard: <(012),(345),(678),(036)(147)(258)>
def cyc(n, *cycles):
    p = list(range(n))
    for c in cycles:
        for k in range(len(c)):
            p[c[k]] = c[(k + 1) % len(c)]
    return tuple(p)

g1 = cyc(9, (0, 1, 2))
g2 = cyc(9, (3, 4, 5))
g3 = cyc(9, (6, 7, 8))
g4 = cyc(9, (0, 3, 6), (1, 4, 7), (2, 5, 8))
P = close_group([g1, g2, g3, g4], 9)
print(f"|<gens>| = {len(P)} (Sylow-3 of S9 has order 81: {len(P) == 81})")
assert len(P) == 81

ts, classes, orbits = invariant_tournaments([g1, g2, g3, g4], 9)
sizes = sorted(len(o) for o in orbits)
print(f"orbits on ordered pairs: {len(orbits)} with sizes {sizes}; classes = {len(classes)} -> {len(ts)} invariant tournaments")
assert sizes == [9, 9, 27, 27] and len(ts) == 4

allT9 = True
for i, adj in enumerate(ts):
    st = is_strong(adj, 9)
    a = aut_count(adj, 9)
    h = ham_count(adj, 9)
    iso = find_iso(adj, T9, 9)
    ok = iso is not None
    if ok:  # verify the iso really is one
        im = apply_perm(adj, 9, iso)
        ok = (im == T9)
    print(f"  invariant #{i}: strong={st} |Aut|={a} H={h} iso-to-T9(verified)={ok}")
    allT9 = allT9 and st and a == 81 and h == 3159 and ok
assert allT9
print("m=9 CONCLUSION: any T with |Aut| in [75,81] has |Aut|=81 = full Sylow-3 order")
print(" => Aut = a Sylow-3 subgroup => T conjugate to a P-invariant tournament => T iso T9.")
print(" T9 is strong, |Aut|=81, H=3159 != 81 (tc=39). NO strong m=9 counterexample. VERIFIED.")

# =========================================================================
print()
print("=" * 78)
print("SECTION 2: m = 27 (window [f(27), o(27)] = [1171875, 1594323])")
print("=" * 78)

def build_T27():
    adj = [0] * 27
    T9l = build_T9()
    for g in range(3):
        base = 9 * g
        for u in range(9):
            m = T9l[u]
            while m:
                b = m & -m
                v = b.bit_length() - 1
                adj[base + u] |= 1 << (base + v)
                m ^= b
    for g in range(3):
        for x in range(9):
            for y in range(9):
                adj[9 * g + x] |= 1 << (9 * ((g + 1) % 3) + y)
    return adj

T27 = build_T27()
assert is_strong(T27, 27)
w1 = cyc(27, (0, 1, 2))
w2 = cyc(27, (0, 3, 6), (1, 4, 7), (2, 5, 8))
w3 = cyc(27, tuple(range(0, 27, 9)), *[tuple(i + j for j in range(0, 27, 9)) for i in range(1, 9)])
# w3 = (0 9 18)(1 10 19)...(8 17 26)
ts27, classes27, orbits27 = invariant_tournaments([w1, w2, w3], 27)
sizes27 = sorted(len(o) for o in orbits27)
print(f"orbits on ordered pairs: {len(orbits27)} sizes {sizes27}; -> {len(ts27)} invariant tournaments")
assert sizes27 == [27, 27, 81, 81, 243, 243] and len(ts27) == 8

# constructive iso: level flips
fin = list(range(27))
for t in range(0, 27, 3):
    fin[t + 1], fin[t + 2] = t + 2, t + 1
fmid = list(range(27))
for s in range(0, 27, 9):
    for i in range(3):
        fmid[s + 3 + i], fmid[s + 6 + i] = s + 6 + i, s + 3 + i
fout = list(range(27))
for i in range(9):
    fout[9 + i], fout[18 + i] = 18 + i, 9 + i

def compose_perm(p, q):  # p after q
    return [p[q[i]] for i in range(len(p))]

all27 = True
for i, adj in enumerate(ts27):
    st = is_strong(adj, 27)
    b_in = 0 if (adj[0] >> 1) & 1 else 1     # T27 has 0->1
    b_mid = 0 if (adj[0] >> 3) & 1 else 1    # T27 has 0->3
    b_out = 0 if (adj[0] >> 9) & 1 else 1    # T27 has 0->9
    perm = list(range(27))
    if b_in:
        perm = compose_perm(fin, perm)
    if b_mid:
        perm = compose_perm(fmid, perm)
    if b_out:
        perm = compose_perm(fout, perm)
    im = apply_perm(adj, 27, perm)
    ok = (im == T27)
    print(f"  invariant #{i}: strong={st} bits=({b_in},{b_mid},{b_out}) explicit-iso-to-T27 verified={ok}")
    all27 = all27 and st and ok
assert all27
print(f"v3(27!) = {v3(fact(27))} (claimed 13)")
assert v3(fact(27)) == 13
f27 = min(3 ** a * 5 ** b for a in range(0, 14) for b in range(0, 10) if 2 * a + 3 * b == 26)
pows = [3 ** k for k in range(1, 14) if f27 <= 3 ** k <= 3 ** 13]
print(f"f(27) = {f27}; 3-powers in [f(27), 3^13]: {pows}")
assert f27 == 1171875 and pows == [3 ** 13]
# GL(3,3) odd part
gl33 = (27 - 1) * (27 - 3) * (27 - 9)
print(f"|GL(3,3)| = {gl33}, affine primitive cap = 27*oddpart = {27 * oddpart(gl33)} < f(27): {27 * oddpart(gl33) < f27}")
assert 27 * oddpart(gl33) == 9477
# S9 odd subgroup order pinch: divisors of 2835 in (59, 81]
divs = [d for d in range(60, 82) if 2835 % d == 0]
print(f"divisors of oddpart(9!)=2835 in (59,81]: {divs}  (63 excluded via N_S9(Z7): order 84, 63 does not divide 84)")
assert divs == [63, 81] and 84 % 63 != 0
print("m=27 CONCLUSION: all 8 Sylow-invariant tournaments are (explicit-iso) the tower T27 =")
print(" C3[T9,T9,T9]: strong but NON-PRIME, so the prime case is vacuous at m=27. VERIFIED.")

# =========================================================================
print()
print("=" * 78)
print("SECTION 3: independent o(m) upper bounds, f(m), crossovers; TRUE lower bounds")
print("=" * 78)

def divisors(n):
    return [d for d in range(2, n) if n % d == 0]

def prime_power(n):
    if n < 2:
        return None
    p = 2
    while p * p <= n:
        if n % p == 0:
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            return (p, k) if n == 1 else None
        p += 1
    return (n, 1)

def gl_order(p, d):
    o = 1
    q = p ** d
    for i in range(d):
        o *= q - p ** i
    return o

@lru_cache(None)
def Tub(n):  # upper bound for odd transitive subgroup of S_n
    if n == 1:
        return 1
    if n % 2 == 0:
        return 0
    pp = prime_power(n)
    best = 0
    if pp:
        p, d = pp
        best = n * oddpart(n - 1) if d == 1 else n * oddpart(gl_order(p, d))
    for d in divisors(n):
        q = n // d
        cand = Oub(d) ** q * Tub(q)
        best = max(best, cand)
    return best

@lru_cache(None)
def Oub(m):  # upper bound for any odd subgroup of S_m
    if m <= 0:
        return 1
    best = 1 if m == 1 else 0
    for j in range(1, m + 1, 2):
        t = Tub(j)
        if t:
            best = max(best, t * Oub(m - j))
    return best

def fB(m):
    return min(3 ** a * 5 ** b for a in range(0, (m + 1) // 2 + 1)
               for b in range(0, m // 3 + 1) if 2 * a + 3 * b == m - 1)

# TRUE achievable lower bounds: transitive witnesses t_true on parts
# 3-towers 3^k: 3^{(3^k-1)/2}; R5: 5; P7: 21; Z3wrZ5 (15): 1215; P7[C3^7] (21): 45927
t_true = {1: 1, 3: 3, 5: 5, 7: 21, 9: 81, 15: 1215, 21: 45927, 27: 3 ** 13, 81: 3 ** 40}
# also 45 = Z3 wr (Z3wrZ5)? blocks of 3 over 15: 3^15 * 1215 = 3^20*5; 63 = blocks of 3 over 21: 3^21*45927=3^28*21? (careful: 63 = 3*21, kernel 3^21, top t(21)=45927 -> 3^21*3^7*21 = 3^28*21)
t_true[45] = 3 ** 20 * 5
t_true[63] = 3 ** 28 * 21
t_true[243] = 3 ** 121

@lru_cache(None)
def Ltrue(m):  # achievable |Aut| lower bound (products of transitive witnesses)
    if m == 0:
        return 1
    best = 0
    for part, val in t_true.items():
        if part <= m:
            best = max(best, val * Ltrue(m - part))
    return best

hdr = f"{'m':>3} {'f(m)':>16} {'oUB(m)':>16} {'Ltrue(m)':>16}  cross?"
print(hdr)
crossUB, crossTRUE = [], []
otable = {}
for m in range(3, 62):
    fm, om, lm = fB(m), Oub(m), Ltrue(m)
    otable[m] = om
    cu = om >= fm
    ct = lm >= fm
    if cu:
        crossUB.append(m)
    if ct:
        crossTRUE.append(m)
    if m <= 28 or cu or ct or m in (45, 54):
        print(f"{m:>3} {fm:>16} {om:>16} {lm:>16}  {'UBcross' if cu else ''}{' TRUEcross' if ct else ''}")
print(f"crossovers (upper bound) on [3,61]: {crossUB}")
print(f"crossovers (TRUE, achievable group) on [3,61]: {crossTRUE}")
assert crossUB == [3, 9, 27, 54]
assert crossTRUE == [3, 9, 27, 54], f"true crossover mismatch: {crossTRUE}"

# claimed o-table m=3..13 and witness values
claimed = {3: 3, 4: 3, 5: 5, 6: 9, 7: 21, 8: 21, 9: 81, 10: 81, 11: 81, 12: 243, 13: 243,
           14: 441, 15: 1215, 16: 1701, 17: 1701, 18: 6561, 19: 6561, 20: 6561, 21: 45927}
mismatch = {m: (otable[m], c) for m, c in claimed.items() if otable[m] != c}
print(f"o-table m=3..21 matches claimed table: {not mismatch} {mismatch if mismatch else ''}")
assert not mismatch
print(f"near-miss margins: m=21: o/f = {otable[21]/fB(21):.4f}; m=45: {otable[45]/fB(45):.4f}; m=28: {otable[28]/fB(28):.4f}; m=36: {otable[36]/fB(36):.4f}")

# m=54 shape exclusion, independently: max over partitions of 54 into odd parts != (27,27)
def max_excl(m):
    best = [0]
    def rec(rem, maxpart, prod, n27):
        if rem == 0:
            if n27 != 2 or prod != Tub(27) ** 2:
                best[0] = max(best[0], prod)
            return
        p = min(maxpart, rem)
        if p % 2 == 0:
            p -= 1
        while p >= 1:
            t = Tub(p)
            if t:
                rec(rem - p, p, prod * t, n27 + (1 if p == 27 else 0))
            p -= 2
    rec(m, m, 1, 0)
    return best[0]

e54 = max_excl(54)
print(f"m=54: max odd order, orbit shape != (27,27): {e54} = 3^24*5? {e54 == 3**24*5}; < f(54) = {fB(54)}: {e54 < fB(54)}")
assert e54 < fB(54) and Oub(54) == 3 ** 26

# intransitive cap at 27
intr27 = max(Tub(j) * Oub(27 - j) for j in range(1, 27, 2) if Tub(j))
print(f"m=27 intransitive cap: {intr27} = 3^12: {intr27 == 3**12}; < f(27): {intr27 < fB(27)}")
assert intr27 == 3 ** 12

# Burnside at prime degrees 5,7,11,13: p*oddpart(p-1) vs f(p)
for p in [5, 7, 11, 13, 17, 19, 23]:
    bb = p * oddpart(p - 1)
    print(f"  prime p={p}: Burnside cap {bb} vs f(p) = {fB(p)}  ok={bb < fB(p)}")
    assert bb < fB(p) or p == 3

# =========================================================================
print()
print("=" * 78)
print("SECTION 4: primitive/affine inequality m^(1+log3 m) < 5^((m-1)/3), m in 25..5000")
print("=" * 78)
viol = [k for k in range(25, 5001)
        if (1 + math.log(k, 3)) * math.log(k) >= (k - 1) / 3 * math.log(5)]
print(f"violations in [25,5000]: {viol}")
assert not viol
k = 25
lhs = (1 + math.log(k, 3)) * math.log(k)
rhs = (k - 1) / 3 * math.log(5)
print(f"margin at k=25 (tightest): exp(lhs) = {math.exp(lhs):.3e} vs 5^8 = {5**8} (ratio {math.exp(lhs)/5**8:.3f})")

# =========================================================================
print()
print("=" * 78)
print("SECTION 5: Theorem R spot checks (lex composition inequalities) + controls")
print("=" * 78)

def lex(Q, nQ, mods):  # mods = list of (adj, size)
    off = []
    tot = 0
    for (a, s) in mods:
        off.append(tot)
        tot += s
    adj = [0] * tot
    for i, (a, s) in enumerate(mods):
        for u in range(s):
            m = a[u]
            while m:
                b = m & -m
                v = b.bit_length() - 1
                adj[off[i] + u] |= 1 << (off[i] + v)
                m ^= b
    for i in range(nQ):
        for j in range(nQ):
            if i != j and (Q[i] >> j) & 1:
                for u in range(mods[i][1]):
                    for v in range(mods[j][1]):
                        adj[off[i] + u] |= 1 << (off[j] + v)
    return adj, tot

C3 = [0b010, 0b100, 0b001]  # 0->1->2->0
pt = ([0], 1)
R5 = [0] * 5
for i in range(5):
    for d in (1, 2):
        R5[i] |= 1 << ((i + d) % 5)
assert is_strong(R5, 5) and aut_count(R5, 5) == 5
print(f"R5: strong, |Aut| = 5, H = {ham_count(R5, 5)}")

ctrl, n5 = lex(C3, 3, [(C3, 3), ([0], 1), ([0], 1)])
Hc, Ac = ham_count(ctrl, n5), aut_count(ctrl, n5)
print(f"C3[C3,pt,pt]: H = {Hc} (canon 15), |Aut| = {Ac} (canon 3), tc = {Hc//Ac}, strong = {is_strong(ctrl, n5)}")
assert Hc == 15 and Ac == 3

import random
random.seed(134)
def rand_tourn(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
    return adj

print("random Theorem-R inequality checks (Q in {C3,R5}, random modules, total <= 9):")
okR = True
for trial in range(60):
    if trial % 2:
        Q, nQ = C3, 3
        sizes = random.choice([[2, 1, 1], [2, 2, 1], [3, 1, 1], [2, 2, 2], [3, 2, 1], [1, 1, 1]])
    else:
        Q, nQ = R5, 5
        sizes = random.choice([[2, 1, 1, 1, 1], [1, 1, 1, 1, 1], [2, 2, 1, 1, 1], [3, 1, 1, 1, 1]])
    mods = [(rand_tourn(s), s) for s in sizes]
    T, nT = lex(Q, nQ, mods)
    HT, AT = ham_count(T, nT), aut_count(T, nT)
    HQ, AQ = ham_count(Q, nQ), aut_count(Q, nQ)
    prodH = HQ
    prodA = AQ
    for (a, s) in mods:
        prodH *= ham_count(a, s)
        prodA *= aut_count(a, s)
    lower_ok = HT >= prodH
    upper_ok = AT <= prodA
    lem_ok = HT % AT == 0 and HT >= AT
    strict_ok = True
    if Q is C3 and max(sizes) >= 2:
        strict_ok = HT > prodH  # R' strictness clause
    okR = okR and lower_ok and upper_ok and lem_ok and strict_ok
    if not (lower_ok and upper_ok and lem_ok and strict_ok):
        print(f"  FAIL trial {trial}: sizes={sizes} HT={HT} prodH={prodH} AT={AT} prodA={prodA}")
print(f"all 60 random R-checks pass (H(T) >= H(Q)*prod, |Aut(T)| <= |Aut(Q)|*prod, LEM-003+ divisibility, R' strictness): {okR}")
assert okR

# =========================================================================
print()
print("=" * 78)
print("SECTION 6: q >= 5 elements at m <= 11 -- families, exact |Aut| for every member")
print("=" * 78)

def family(name, m, cycles, fexp):
    sigma = cyc(m, *cycles)
    ts, cls, orbs = invariant_tournaments([sigma], m)
    mx = 0
    mxs = 0
    nstr = 0
    for adj in ts:
        a = aut_count(adj, m)
        st = is_strong(adj, m)
        mx = max(mx, a)
        if st:
            nstr += 1
            mxs = max(mxs, a)
    print(f"  {name:<28} members={len(ts):>5}  max|Aut|={mx:>4}  max|Aut| over strong={mxs:>4} (#strong={nstr})  f({m})={fexp}  all < f: {mx < fexp}")
    return mx

mx1 = family("m=11 (5,5,1)", 11, [tuple(range(5)), tuple(range(5, 10))], fB(11))
mx2 = family("m=11 (5,3,3)", 11, [tuple(range(5)), (5, 6, 7), (8, 9, 10)], fB(11))
mx3 = family("m=11 (7,3,1)", 11, [tuple(range(7)), (7, 8, 9)], fB(11))
mx4 = family("m=11 (11) circulants", 11, [tuple(range(11))], fB(11))
mx5 = family("m=10 (5,5)", 10, [tuple(range(5)), tuple(range(5, 10))], fB(10))
mx6 = family("m=10 (7,3)", 10, [tuple(range(7)), (7, 8, 9)], fB(10))
mx7 = family("m=9  (5,3,1)", 9, [tuple(range(5)), (5, 6, 7)], fB(9))
mx8 = family("m=9  (7,1,1)", 9, [tuple(range(7))], fB(9))
mx9 = family("m=8  (5,3)", 8, [tuple(range(5)), (5, 6, 7)], fB(8))
mx10 = family("m=8  (7,1)", 8, [tuple(range(7))], fB(8))
assert all(v < 225 for v in [mx1, mx2, mx3, mx4]) and mx5 < 125 and mx6 < 125 and mx7 < 75 and mx8 < 75 and mx9 < 45 and mx10 < 45
print("NOTE: remaining cycle types with q>=5 at m<=11 ((5,1^k),(9,1^k),(5,3,1^k)...) are covered by")
print("the orbit-product bound: any odd G <= S11 with an orbit of size 5/7/9/11 has")
print(f"|G| <= max(5*Oub(6), 21*Oub(4), 81*Oub(2), 55) = {max(5*Oub(6), 21*Oub(4), 81*Oub(2), 55)} < f(11) = {fB(11)}")

# =========================================================================
print()
print("=" * 78)
print("SECTION 7: BEYOND 61 -- the TRUE crossover landscape (how big is the ML burden?)")
print("=" * 78)
crossT = [m for m in range(3, 401) if Ltrue(m) >= fB(m)]
noncross = [m for m in range(62, 401) if Ltrue(m) < fB(m)]
print(f"TRUE crossovers (achievable odd group order >= Busch floor) in [3,400]:")
print(f"  {crossT}")
print(f"largest NON-crossover m in [62,400]: {max(noncross) if noncross else 'none'}")
print(f"#non-crossovers in [62,400]: {len(noncross)}")
big = [m for m in range(200, 401) if Ltrue(m) < fB(m)]
print(f"non-crossovers in [200,400]: {big}")
print()
print("ALL CHECKS PASSED")
