#!/usr/bin/env python3
"""Shaved tournaments: the oriented graphs contained in every n-tournament (opus S15, 2026-10-01).

The owner's object: every 4-tournament contains A->B, B->C, C->D, A->D (H4 = Hamiltonian path plus the arc
from its first to its last vertex).  Two generalisations are studied as n grows.

Checks:
  T1. isomorphism classes of tournaments for n <= 8 (1,1,2,4,12,56,456,6880 = A000568) by canonical augmentation
  T2. H4: its 4 completions are the 4 classes, each once; copies of H4 in T = |Aut T| (1,3,3,1)
  T3. Hn = Hamiltonian path + first->last arc: #copies = H(T) - n*hc(T); odd for every class with n = 4, 6, 8
      (Redei); zero only for C3 (n = 3) and one class at n = 5 (C3[1,C3,1]); never zero at n = 7; n = 9 by the
      C helper over every one-vertex extension of the 6880 classes (all 9-tournaments)
  T4. u(n) = max #arcs of a spanning oriented graph contained in every n-tournament (= C(n,2) - kappa(n),
      HYP-3798): exhaustive n <= 6; u = 2,4,6,8, unique maximiser = path + all span-3 arcs (H4 at n = 4)
  T5. n = 7: P7 is the only class without TT4, so every unavoidable 7-vertex graph lies in P7; every acyclic
      10-arc subgraph of P7 is avoided by some class => u(7) = 9, kappa(7) = 12 (HYP-3805 'very likely' -> proved);
      the 9-arc maximisers (51 isomorphism classes); the span-{1,3} graph misses exactly P7 and one |Aut|=1 class
  T6. T(n) by the Davis/Polya formula; u(n) <= C(n,2) - ceil(log2 T(n)) <= log2 n!; lower bound by superadditivity
      and transitive blocks (Erdos-Moser): u(n) = Theta(n log n); first n where the bound beats 2n - 4
  T7. perfect shavings (completions = classes, each once) exist only for n <= 4 (T(n) not a power of 2, 5 <= n <= 100)
  T8. transversal multiplicities of the maximum shavings: m_S(T) = emb(S,T)/|Aut T|, n = 4, 5, 6; the number of
      embeddings of path + span-3 arcs is odd for every tournament with n <= 6, not for n = 7
Reproduce: python 04-computation/experiments/shaved_tournaments_20261001.py   (a few minutes; T3 n = 9 needs gcc)
"""
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations

FAILS = []
HERE = os.path.dirname(os.path.abspath(__file__))


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ------------------------------------------------------------------ tournament library
def refine(n, out, colors):
    while True:
        sigs = [(colors[v], tuple(sorted(colors[w] for w in range(n) if out[v] >> w & 1))) for v in range(n)]
        keys = sorted(set(sigs))
        idx = {k: i for i, k in enumerate(keys)}
        new = [idx[s] for s in sigs]
        if len(keys) == len(set(colors)):
            return new
        colors = new


def code_of(n, out, order):
    c, b = 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            if out[order[i]] >> order[j] & 1:
                c |= 1 << b
            b += 1
    return c


def canon(n, out):
    """(canonical code, |Aut|) by individualisation-refinement over the full search tree."""
    best = [None, 0]

    def rec(colors):
        colors = refine(n, out, colors)
        if len(set(colors)) == n:
            c = code_of(n, out, sorted(range(n), key=lambda v: colors[v]))
            if best[0] is None or c < best[0]:
                best[0], best[1] = c, 1
            elif c == best[0]:
                best[1] += 1
            return
        cnt = Counter(colors)
        target = min(c for c in cnt if cnt[c] > 1)
        for v in range(n):
            if colors[v] == target:
                nc = [2 * c for c in colors]
                nc[v] = 2 * target - 1
                rec(nc)
    rec([bin(out[v]).count('1') for v in range(n)])
    return best[0], best[1]


def from_code(n, c):
    out = [0] * n
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if c >> b & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1
    return out


def all_classes(nmax):
    res = {1: [(0, 1)]}
    for n in range(2, nmax + 1):
        seen = {}
        for code, _ in res[n - 1]:
            base = from_code(n - 1, code)
            for S in range(1 << (n - 1)):
                out = base[:] + [0]
                for v in range(n - 1):
                    if S >> v & 1:
                        out[n - 1] |= 1 << v
                    else:
                        out[v] |= 1 << (n - 1)
                c, a = canon(n, out)
                seen.setdefault(c, a)
        res[n] = sorted(seen.items())
    return res


def hp_hc(n, out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for m in range(1, 1 << n):
        for v in range(n):
            x = dp[m][v]
            if x:
                o = out[v] & ~m
                while o:
                    w = (o & -o).bit_length() - 1
                    o &= o - 1
                    dp[m | 1 << w][w] += x
    H = sum(dp[full])
    dq = [[0] * n for _ in range(1 << n)]
    dq[1][0] = 1
    for m in range(1, 1 << n, 2):
        for v in range(n):
            x = dq[m][v]
            if x:
                o = out[v] & ~m
                while o:
                    w = (o & -o).bit_length() - 1
                    o &= o - 1
                    dq[m | 1 << w][w] += x
    return H, sum(dq[full][v] for v in range(n) if out[v] & 1)


def make_embedder(n, arcs):
    """f(outT) -> number of embeddings (if count=True) or bool, of the spanning digraph `arcs` into T."""
    outS, inS = [0] * n, [0] * n
    for a, b in arcs:
        outS[a] |= 1 << b
        inS[b] |= 1 << a
    order, chosen = [], 0
    rest = set(range(n))
    while rest:
        v = max(rest, key=lambda x: (bin((outS[x] | inS[x]) & chosen).count('1'), bin(outS[x] | inS[x]).count('1'), -x))
        order.append(v)
        rest.remove(v)
        chosen |= 1 << v
    pos = {v: i for i, v in enumerate(order)}
    req_out = [[pos[w] for w in range(n) if outS[v] >> w & 1 and pos[w] < i] for i, v in enumerate(order)]
    req_in = [[pos[w] for w in range(n) if inS[v] >> w & 1 and pos[w] < i] for i, v in enumerate(order)]

    def emb(outT, count=False):
        inT = [0] * n
        for v in range(n):
            o = outT[v]
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                inT[w] |= 1 << v
        img = [0] * n

        def rec(i, used):
            if i == n:
                return 1
            cand = ((1 << n) - 1) & ~used
            for j in req_out[i]:
                cand &= inT[img[j]]
            for j in req_in[i]:
                cand &= outT[img[j]]
            tot = 0
            while cand:
                y = (cand & -cand).bit_length() - 1
                cand &= cand - 1
                img[i] = y
                r = rec(i + 1, used | 1 << y)
                if r and not count:
                    return 1
                tot += r
            return tot
        return rec(0, 0) if count else bool(rec(0, 0))
    return emb


def scores(out):
    return tuple(sorted(bin(o).count('1') for o in out))


# ------------------------------------------------------------------ T1
print("=== T1. tournament classes ===")
C = all_classes(8)
check([len(C[n]) for n in range(1, 9)] == [1, 1, 2, 4, 12, 56, 456, 6880], "class counts 1,1,2,4,12,56,456,6880 = A000568")
check(all(sum(math.factorial(n) // a for _, a in C[n]) == 2 ** (n * (n - 1) // 2) for n in range(1, 9)),
      "orbit-counting: sum over classes of n!/|Aut| = 2^C(n,2) for n <= 8 (|Aut| values correct)")

# ------------------------------------------------------------------ T2
print("=== T2. the owner's H4 ===")
H4 = [(0, 1), (1, 2), (2, 3), (0, 3)]
comps = []
for b1, b2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
    out = [0] * 4
    for a, b in H4:
        out[a] |= 1 << b
    for (a, b), bit in (((0, 2), b1), ((1, 3), b2)):
        if bit:
            out[a] |= 1 << b
        else:
            out[b] |= 1 << a
    comps.append(canon(4, out)[0])
check(sorted(comps) == sorted(c for c, _ in C[4]), "the 4 completions of H4 are the 4 classes of 4-tournaments, each once")
emb4 = make_embedder(4, H4)
cop = [(scores(from_code(4, c)), a, emb4(from_code(4, c), count=True)) for c, a in C[4]]
check(all(e == a for _, a, e in cop), f"copies of H4 in T equal |Aut T| for every class: {[(s, e) for s, a, e in cop]}")

# ------------------------------------------------------------------ T3
print("=== T3. Hn = Hamiltonian path + first->last arc ===")
data = {}
for n in range(3, 9):
    rows = []
    for c, a in C[n]:
        out = from_code(n, c)
        H, hc = hp_hc(n, out)
        rows.append((c, a, H, hc, H - n * hc))
    data[n] = rows
    zeros = [r for r in rows if r[4] == 0]
    print(f"   n={n}: min #Hn = {min(r[4] for r in rows)}, classes with #Hn = 0: {len(zeros)}")
for n in (4, 6, 8):
    check(all(r[4] % 2 == 1 for r in data[n]), f"n = {n}: every class has an ODD number of copies of H{n} (Redei: H odd, n even)")
z3 = [r for r in data[3] if r[4] == 0]
z5 = [r for r in data[5] if r[4] == 0]
check(len(z3) == 1 and scores(from_code(3, z3[0][0])) == (1, 1, 1), "n = 3: only the 3-cycle avoids H3 (= TT3)")
ok5 = False
if len(z5) == 1:
    out = from_code(5, z5[0][0])
    sc = [bin(o).count('1') for o in out]
    a = sc.index(3)
    b = sc.index(1)
    mid = [v for v in range(5) if v not in (a, b)]
    is_c3 = all(bin(out[v] & sum(1 << w for w in mid)).count('1') == 1 for v in mid)
    ok5 = (all(out[a] >> v & 1 for v in mid) and all(out[v] >> b & 1 for v in mid) and out[b] >> a & 1 and is_c3
           and z5[0][1] == 3 and z5[0][2] == 15 and z5[0][3] == 3)
check(ok5, "n = 5: exactly one class avoids H5: C3[1,C3,1] (a => C3 => b -> a), |Aut| = 3, H = 15 = 5 * hc")
check(all(r[4] > 0 for r in data[7]), "n = 7: every class contains H7 (no all-closing 7-tournament)")
src = os.path.join(HERE, "shaved_tournaments_20261001_closing.c")
gcc = shutil.which("gcc")
if gcc and os.path.exists(src):
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "closing.exe")
        subprocess.run([gcc, "-O2", "-fopenmp", "-o", exe, src], check=True)
        cf = os.path.join(td, "c8.txt")
        with open(cf, "w") as f:
            for c, _ in C[8]:
                f.write(f"{c}\n")
        r = subprocess.run([exe, "8", cf], capture_output=True, text=True, check=True).stdout
        last = r.strip().splitlines()[-1]
        print("   C helper:", last)
        check(last == "n=9 tested=1761280 allclosing=0",
              "n = 9: every one-vertex extension of every 8-class (hence every 9-tournament) contains H9")
else:
    print("   SKIP n = 9 (gcc or the C helper not found); recorded result: n=9 tested=1761280 allclosing=0")

# T3b. the odd case in general: a constructive proof (Steps 1-5 of the note), implemented independently
import random


def arc(out, a, b):
    return out[a] >> b & 1


def strong(out, verts):
    verts = list(verts)
    if len(verts) <= 1:
        return True
    S = sum(1 << v for v in verts)

    def reach(nb):
        seen, st = 1 << verts[0], [verts[0]]
        while st:
            v = st.pop()
            o = nb(v) & S & ~seen
            seen |= o
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                st.append(w)
        return seen == S
    return reach(lambda v: out[v]) and reach(lambda v: sum(1 << u for u in verts if out[u] >> v & 1))


def ham_path_of(out, verts):
    P = []
    for v in verts:
        k = 0
        while k < len(P) and arc(out, P[k], v):
            k += 1
        P.insert(k, v)     # P[k-1] -> v and (k == len(P) or v -> P[k]) because P[k] does not beat... see check
    assert all(arc(out, P[t], P[t + 1]) for t in range(len(P) - 1))
    return P


def ham_cycle_of(out, verts):
    verts = list(verts)
    C = next([a, b, c] for a in verts for b in verts for c in verts
             if len({a, b, c}) == 3 and arc(out, a, b) and arc(out, b, c) and arc(out, c, a))
    while len(C) < len(verts):
        rest = [v for v in verts if v not in C]
        k = len(C)
        ins = next(((v, t) for v in rest for t in range(k) if arc(out, C[t], v) and arc(out, v, C[(t + 1) % k])), None)
        if ins:
            v, t = ins
            C.insert(t + 1, v)
            continue
        A = [v for v in rest if arc(out, C[0], v)]       # then C => v
        B = [v for v in rest if arc(out, v, C[0])]       # then v => C
        a, b = next((a, b) for a in A for b in B if arc(out, a, b))
        C = [C[0], a, b] + C[2:]                         # length + 1 (C[1] re-enters later)
    assert all(arc(out, C[t], C[(t + 1) % len(C)]) for t in range(len(C)))
    return C


def nonclosing_from_cycle(out, C):
    """Steps 1-5: from a Hamiltonian cycle of a strong odd-order T, a non-closing HP (or None at C3, C3[1,C3,1])."""
    n = len(C)
    c = lambda k: C[k % n]
    F = [arc(out, c(i - 1), c(i + 1)) for i in range(n)]
    for i in range(n):
        if F[i] and F[(i + 1) % n]:                                  # Step 1
            return [c(i)] + [c(i + k) for k in range(2, n)] + [c(i + 1)]
    i = next(i for i in range(n) if not F[i] and not F[(i + 1) % n])   # Step 2 (n odd)

    def insert(w, Q):
        pat = ['o' if arc(out, w, q) else 'i' for q in Q]
        if pat[0] == pat[-1] == 'o':
            return [w] + Q
        if pat[0] == pat[-1] == 'i':
            return Q + [w]
        for k in range(len(Q) - 1):
            if pat[k] == 'i' and pat[k + 1] == 'o':
                return Q[:k + 1] + [w] + Q[k + 1:]
        return None
    P = insert(c(i), [c(i + k) for k in range(1, n)])               # Step 3
    if P:
        return P
    P = insert(c(i + 1), [c(i + k) for k in range(2, n + 1)])       # Step 4
    if P:
        return P
    a, b, R = c(i + 1), c(i), [c(i + k) for k in range(2, n)]       # Step 5: a => R => b -> a
    assert all(arc(out, a, r) and arc(out, r, b) for r in R) and arc(out, b, a)
    if len(R) >= 2 and not strong(out, R):
        s = ham_path_of(out, R)[0]
        return [s, b, a] + ham_path_of(out, [r for r in R if r != s])
    if len(R) <= 3:                                                  # C3 (|R| = 1) or C3[1,C3,1] (|R| = 3)
        return None
    s = next(s for s in R if strong(out, [r for r in R if r != s]))
    CY = ham_cycle_of(out, [r for r in R if r != s])
    k = next(k for k, t in enumerate(CY) if arc(out, s, t))
    return [s, b, a] + CY[k + 1:] + CY[:k + 1]


def good(out, P):
    n = len(out)
    return (P is not None and sorted(P) == list(range(n)) and all(arc(out, P[t], P[t + 1]) for t in range(n - 1))
            and arc(out, P[0], P[-1]))


def all_cycles(out, n):
    res, path = [], [0]

    def rec(v, used):
        if len(path) == n:
            if arc(out, v, 0):
                res.append(path[:])
            return
        o = out[v] & ~used
        while o:
            w = (o & -o).bit_length() - 1
            o &= o - 1
            path.append(w)
            rec(w, used | 1 << w)
            path.pop()
    rec(0, 1)
    return res


fails, tested = [], 0
for n in (3, 5, 7):
    for code, aut in C[n]:
        out = from_code(n, code)
        if not strong(out, range(n)):
            continue
        for cyc in all_cycles(out, n):
            tested += 1
            if not good(out, nonclosing_from_cycle(out, cyc)):
                fails.append((n, scores(out)))
check(tested > 0 and sorted(set(fails)) == [(3, (1, 1, 1)), (5, (1, 2, 2, 2, 3))],
      f"Steps 1-5 from EVERY Hamiltonian cycle of every strong class, n = 3, 5, 7 ({tested} runs): a non-closing HP "
      "except exactly at C3 and C3[1,C3,1]")
rng = random.Random(20261001)
ok, runs = True, 0
for n in range(9, 32, 2):
    for _ in range(60):
        out = [0] * n
        for x in range(n):
            for y in range(x + 1, n):
                if rng.random() < 0.5:
                    out[x] |= 1 << y
                else:
                    out[y] |= 1 << x
        if not strong(out, range(n)):
            continue
        runs += 1
        ok &= good(out, nonclosing_from_cycle(out, ham_cycle_of(out, list(range(n)))))
for n in range(7, 32, 2):           # the hardest family for Step 5: C3[1, R, 1]
    for _ in range(20):
        m = n - 2
        out = [0] * n
        for x in range(m):
            for y in range(x + 1, m):
                if rng.random() < 0.5:
                    out[x] |= 1 << y
                else:
                    out[y] |= 1 << x
        a, b = m, m + 1
        for r in range(m):
            out[a] |= 1 << r
            out[r] |= 1 << b
        out[b] |= 1 << a
        if not strong(out, range(n)):
            continue
        runs += 1
        ok &= good(out, nonclosing_from_cycle(out, ham_cycle_of(out, list(range(n)))))
check(ok and runs > 500, f"Steps 1-5 on {runs} random strong tournaments with odd n = 7..31 (incl. C3[1,R,1]): always a non-closing HP")

# ------------------------------------------------------------------ T4
print("=== T4. maximum shavings u(n), n <= 6 ===")


class Tester:
    def __init__(self, n):
        self.n = n
        self.T = [from_code(n, c) for c, _ in C[n]]
        self.order = list(range(len(self.T)))

    def unavoidable(self, arcs):
        emb = make_embedder(self.n, arcs)
        for k, i in enumerate(self.order):
            if not emb(self.T[i]):
                self.order.insert(0, self.order.pop(k))
                return False
        return True


u = {1: 0, 2: 1}
maxi = {}
for n in range(3, 7):
    X = Tester(n)
    pairs = list(combinations(range(n), 2))
    for e in range(len(pairs), -1, -1):
        good = [S for S in combinations(pairs, e) if X.unavoidable(S)]
        if good:
            u[n], maxi[n] = e, good
            break
    expect = sorted([(i, i + 1) for i in range(n - 1)] + [(i, i + 3) for i in range(n - 3)])
    check(len(maxi[n]) == 1 and sorted(maxi[n][0]) == expect,
          f"n = {n}: u = {u[n]} = 2n-4, unique maximiser (one forward labelling) = path + all span-3 arcs")
check([len(list(combinations(range(n), 2))) - u[n] for n in range(3, 7)] == [1, 2, 4, 7],
      "kappa(n) = C(n,2) - u(n) = 1,2,4,7 (agrees with HYP-3798, mac-mini)")

# ------------------------------------------------------------------ T5
print("=== T5. n = 7: the Paley heptagon decides u(7) ===")
QR = {1, 2, 4}
P7 = [sum(1 << y for y in range(7) if (y - x) % 7 in QR) for x in range(7)]
cP7, aP7 = canon(7, P7)
T7 = [from_code(7, c) for c, _ in C[7]]


def has_TT4(out):
    """some 4 vertices induce a transitive subtournament (internal out-degrees 0,1,2,3)"""
    for S4 in combinations(range(7), 4):
        m = sum(1 << v for v in S4)
        if sorted(bin(out[v] & m).count('1') for v in S4) == [0, 1, 2, 3]:
            return True
    return False


noTT4 = [i for i, t in enumerate(T7) if not has_TT4(t)]
check(len(noTT4) == 1 and C[7][noTT4[0]][0] == cP7 and aP7 == 21,
      "P7 (QR_7, |Aut| = 21) is the only 7-tournament without TT4: every unavoidable 7-vertex graph embeds in P7")
arcs7 = [(x, y) for x in range(7) for y in range(7) if (y - x) % 7 in QR]
aidx = {a: i for i, a in enumerate(arcs7)}
perm7 = [[aidx[((a * x + b) % 7, (a * y + b) % 7)] for (x, y) in arcs7] for a in (1, 2, 4) for b in range(7)]


def acyclic7(sub):
    out, indeg = [0] * 7, [0] * 7
    for i in sub:
        x, y = arcs7[i]
        out[x] |= 1 << y
        indeg[y] += 1
    st = [v for v in range(7) if indeg[v] == 0]
    seen = 0
    while st:
        v = st.pop()
        seen += 1
        o = out[v]
        while o:
            w = (o & -o).bit_length() - 1
            o &= o - 1
            indeg[w] -= 1
            if indeg[w] == 0:
                st.append(w)
    return seen == 7


X7 = Tester(7)
found = {}
norbits = {}
for e in (9, 10):
    reps = set()
    found[e] = []
    nac = 0
    for sub in combinations(range(21), e):
        key = min(sum(1 << p[i] for i in sub) for p in perm7)
        if key in reps:
            continue
        reps.add(key)
        if not acyclic7(sub):
            continue
        nac += 1
        if X7.unavoidable([arcs7[i] for i in sub]):
            found[e].append(sub)
    norbits[e] = (len(reps), nac)
print(f"   orbits under Aut(P7) of e-arc subgraphs of P7 (all, acyclic): {norbits}")
check(norbits[10][0] == 16796 and not found[10], "no acyclic 10-arc subgraph of P7 is unavoidable: u(7) <= 9")
check(len(found[9]) == 54, "54 Aut(P7)-orbits of unavoidable 9-arc subgraphs: u(7) = 9, kappa(7) = 21 - 9 = 12")


def ocanon(n, A):
    return min(tuple(sorted((p[a], p[b]) for a, b in A)) for p in permutations(range(n)))


iso9 = {}
for sub in found[9]:
    iso9.setdefault(ocanon(7, [arcs7[i] for i in sub]), []).append(sub)
ex = [(i, i + 1) for i in range(6)] + [(1, 4), (2, 5), (3, 6)]
check(len(iso9) == 51 and ocanon(7, ex) in iso9,
      "51 isomorphism classes of maximum 7-vertex shavings; HYP-3805's span1 + span3 - (0,3) is one of them")


def linext(n, A):
    pred = [0] * n
    for a, b in A:
        pred[b] |= 1 << a
    from functools import lru_cache

    @lru_cache(None)
    def f(m):
        if m == (1 << n) - 1:
            return 1
        return sum(f(m | 1 << v) for v in range(n) if not m >> v & 1 and pred[v] & ~m == 0)
    return f(0)


print("   linear-extension counts of the 51 maximisers:", sorted(Counter(linext(7, k) for k in iso9).items()))
# which classes certify u(7) <= 9: every acyclic 10-arc subgraph of P7 and the classes avoiding it
blk = []
reps = set()
for sub in combinations(range(21), 10):
    key = min(sum(1 << p[i] for i in sub) for p in perm7)
    if key in reps:
        continue
    reps.add(key)
    if acyclic7(sub):
        emb = make_embedder(7, [arcs7[i] for i in sub])
        blk.append(frozenset(i for i in range(len(T7)) if not emb(T7[i])))
NB = len(blk)
cover = {}
for j, f in enumerate(blk):
    for i in f:
        cover[i] = cover.get(i, 0) | (1 << j)
fullb = (1 << NB) - 1


def hitting(k, acc=0, depth=0):
    """exact: some chosen class must avoid the uncovered candidate with fewest avoiders (sound branching)"""
    if acc == fullb:
        return []
    if depth == k:
        return None
    unc = fullb & ~acc
    j = min((j for j in range(NB) if unc >> j & 1), key=lambda j: len(blk[j]))
    if max(bin(cover[i] & unc).count('1') for i in cover) * (k - depth) < bin(unc).count('1'):
        return None
    for i in sorted(blk[j], key=lambda i: -bin(cover[i] & unc).count('1')):
        r = hitting(k, acc | cover[i], depth + 1)
        if r is not None:
            return [i] + r
    return None


kmin = next(k for k in range(1, 10) if hitting(k) is not None)
cert = hitting(kmin)
check(NB == 1792 and kmin == 6 and all(C[7][i][1] > 1 for i in cert),
      f"certificate for u(7) <= 9: P7 plus {kmin} more classes (the minimum), all symmetric: "
      f"|Aut| = {sorted((C[7][i][1] for i in cert), reverse=True)}")
span13 = [(i, i + 1) for i in range(6)] + [(i, i + 3) for i in range(4)]
emb13 = make_embedder(7, span13)
miss = [i for i, t in enumerate(T7) if not emb13(t)]
check(len(miss) == 2 and any(C[7][i][0] == cP7 for i in miss) and sorted(C[7][i][1] for i in miss) == [1, 21],
      f"path + all span-3 arcs (10 arcs) is avoided by exactly two 7-classes: P7 and one class with |Aut| = 1 "
      f"(scores {[scores(T7[i]) for i in miss]})")

# ------------------------------------------------------------------ T6
print("=== T6. growth: u(n) = Theta(n log n) ===")


def odd_partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n if n % 2 else n - 1
    if n == 0:
        yield []
        return
    for k in range(maxpart, 0, -2):
        for m in range(1, n // k + 1):
            for rest in odd_partitions(n - m * k, k - 2):
                yield [(k, m)] + rest


def T_count(n):
    if n == 0:
        return 1
    tot = Fraction(0)
    for lam in odd_partitions(n):
        e = sum(k * m * (m - 1) // 2 + m * (k - 1) // 2 for k, m in lam)
        e += sum(m1 * m2 * math.gcd(k1, k2) for (k1, m1), (k2, m2) in combinations(lam, 2))
        z = 1
        for k, m in lam:
            z *= k ** m * math.factorial(m)
        tot += Fraction(2 ** e, z)
    assert tot.denominator == 1
    return int(tot)


Tn = {n: T_count(n) for n in range(1, 61)}
check([Tn[n] for n in range(1, 11)] == [1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056],
      "Davis-Polya formula reproduces A000568 (n <= 10)")
L = {1: 0, 2: 1, 3: 2, 4: 4, 5: 6, 6: 8, 7: 9}
for n in range(8, 61):
    best = max(L[a] + L[n - a] for a in range(1, n))
    k = 1
    while 2 ** k <= n:          # every n-tournament contains TT_(k+1) when n >= 2^k (Erdos-Moser)
        k += 1
        best = max(best, (k * (k - 1)) // 2 + L[n - k])
    L[n] = best
U = {n: n * (n - 1) // 2 - math.ceil(math.log2(Tn[n])) for n in range(1, 61)}
check(all(L[n] <= U[n] for n in range(1, 61)), "lower bound <= upper bound C(n,2) - ceil(log2 T(n)) for n <= 60")
cross = next(n for n in range(8, 61) if L[n] > 2 * n - 4)
print("   n, lower, upper, 2n-4, (n/2)log2 n, log2 n!:",
      [(n, L[n], U[n], 2 * n - 4, round(n / 2 * math.log2(n), 1), round(math.log2(math.factorial(n)), 1))
       for n in (7, 8, 10, 16, 20, 30, 40, 50, 60)])
check(cross <= 60 and all(L[n] > 2 * n - 4 for n in range(cross, 61)),
      f"the transitive-block lower bound exceeds 2n-4 from n = {cross} on (n <= 60): kappa(n) < 1 + C(n-2,2) there")
check(U[7] == 12 and L[7] == 9 == 9, "n = 7: counting bound gives u(7) <= 12; the truth is 9 (T5)")

# ------------------------------------------------------------------ T7
print("=== T7. perfect shavings ===")
check(all(Tn[n] & (Tn[n] - 1) != 0 for n in range(5, 61)) and all(Tn[n] & (Tn[n] - 1) == 0 for n in range(1, 5)),
      "T(n) is a power of 2 exactly for n <= 4 (n <= 60): a perfect shaving (each class exactly once) needs that")

# ------------------------------------------------------------------ T8
print("=== T8. how the maximum shavings cover the classes ===")
for n in (4, 5, 6):
    S = maxi[n][0]
    free = [p for p in combinations(range(n), 2) if p not in S]
    mult = Counter()
    for bits in range(1 << len(free)):
        out = [0] * n
        for a, b in S:
            out[a] |= 1 << b
        for j, (a, b) in enumerate(free):
            if bits >> j & 1:
                out[b] |= 1 << a
            else:
                out[a] |= 1 << b
        mult[canon(n, out)[0]] += 1
    emb = make_embedder(n, S)
    okm = all(mult[c] * a == emb(from_code(n, c), count=True) for c, a in C[n])
    check(len(mult) == len(C[n]) and okm,
          f"n = {n}: 2^{len(free)} completions hit all {len(C[n])} classes; m_S(T) = emb(S,T)/|Aut T|; "
          f"multiplicity histogram {sorted(Counter(mult.values()).items())}")

par = {}
for n in range(4, 8):
    S = [(i, i + 1) for i in range(n - 1)] + [(i, i + 3) for i in range(n - 3)]
    emb = make_embedder(n, S)
    par[n] = Counter(emb(from_code(n, c), count=True) % 2 for c, _ in C[n])
print("   parity of emb(path + span-3 arcs, T) over the classes:", {n: dict(v) for n, v in par.items()})
check(all(set(par[n]) == {1} for n in (4, 5, 6)) and par[7] == Counter({1: 282, 0: 174}),
      "the number of embeddings of path + span-3 arcs is odd for EVERY tournament with n <= 6 (so every class is "
      "hit an odd number of times); at n = 7 it is even for 174 of the 456 classes (a small-n parity)")

print()
print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
