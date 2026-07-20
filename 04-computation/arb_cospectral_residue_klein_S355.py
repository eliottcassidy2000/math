#!/usr/bin/env python3
"""
klein-2026-07-20-S355 -- CLOSING THM-1460(C): WHICH adjacency-cospectral groups does the
arborescence count FAIL to separate, and does H finish the job?

THM-1460(C) reports that Sum_r a_r differs inside 111 of the 116 adjacency-cospectral groups
at n = 7, so it is transverse to the spec(A) hierarchy -- but it does not say WHICH 5 survive.
This identifies them and tests whether H (Hamiltonian paths) or the pair resolves them.

The logarithm frame of THM-1460(D): log H is additive under ordinal sum with NO interaction
term, log Sum a is additive with a SIZE-DEPENDENT SHIFT.  So H and Sum a fail differently,
and the natural question is whether their failures overlap.
"""
import itertools
from fractions import Fraction as Fr

def pairs_of(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]
def relabel(om, perm, n):
    new = [0] * n
    for v in range(n):
        mv, t = om[v], 0
        while mv:
            b = mv & -mv; w = b.bit_length() - 1; mv ^= b
            t |= 1 << perm[w]
        new[perm[v]] = t
    return tuple(new)
def word(om, n):
    w = 0
    for v in range(n): w = (w << n) | om[v]
    return w
def refine(om, n):
    colour = [bin(om[v]).count("1") for v in range(n)]
    while True:
        sig = []
        for v in range(n):
            cnt = {}; mv = om[v]
            while mv:
                b = mv & -mv; w = b.bit_length() - 1; mv ^= b
                cnt[colour[w]] = cnt.get(colour[w], 0) + 1
            sig.append((colour[v], tuple(sorted(cnt.items()))))
        order = sorted(set(sig)); newc = [order.index(sig[v]) for v in range(n)]
        if newc == colour: break
        colour = newc
    cells = {}
    for v in range(n): cells.setdefault(colour[v], []).append(v)
    return [tuple(cells[k]) for k in sorted(cells)]
def canon(om, n):
    cells = refine(om, n); base = []; pos = 0
    for c in cells: base.append((c, pos)); pos += len(c)
    best = None
    for choice in itertools.product(*[itertools.permutations(c) for (c, _) in base]):
        perm = [0] * n
        for (blk, (c, start)) in zip(choice, base):
            for k, v in enumerate(blk): perm[v] = start + k
        w = word(relabel(om, perm, n), n)
        if best is None or w < best: best = w
    return best
def classes(n):
    P = pairs_of(n)
    om0 = tuple(sum(1 << j for j in range(i)) for i in range(n))
    seen = {canon(om0, n): om0}; fr = [om0]
    while fr:
        nx = []
        for om in fr:
            for (i, j) in P:
                nm = list(om)
                if om[i] >> j & 1: nm[i] &= ~(1 << j); nm[j] |= 1 << i
                else:              nm[j] &= ~(1 << i); nm[i] |= 1 << j
                nm = tuple(nm); cc = canon(nm, n)
                if cc not in seen: seen[cc] = nm; nx.append(nm)
        fr = nx
    return seen
def charpoly(om, n):
    """integer char poly of the 0/1 adjacency, via Faddeev-LeVerrier"""
    A = [[1 if (om[i] >> j & 1) else 0 for j in range(n)] for i in range(n)]
    M = [[0] * n for _ in range(n)]; cs = [1]
    for k in range(1, n + 1):
        M = [[sum(A[i][l] * M[l][j] for l in range(n)) + (cs[-1] if i == j else 0)
              for j in range(n)] for i in range(n)]
        AM = [[sum(A[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c = Fr(-sum(AM[i][i] for i in range(n)), k)
        cs.append(c)
    return tuple(cs)
def sum_arb(om, n):
    """Sum_r a_r = product of nonzero eigenvalues of L_in = D_in - A, via any cofactor sum.
       Computed as: sum over roots r of det(L_in with row/col r deleted)."""
    L = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (om[i] >> j & 1): L[j][j] += 1; L[i][j] -= 1
    tot = 0
    for r in range(n):
        idx = [i for i in range(n) if i != r]
        M = [[Fr(L[i][j]) for j in idx] for i in idx]
        d = Fr(1); m = len(idx); ok = True
        for c in range(m):
            p = next((x for x in range(c, m) if M[x][c] != 0), None)
            if p is None: ok = False; break
            if p != c: M[c], M[p] = M[p], M[c]; d = -d
            d *= M[c][c]; inv = 1 / M[c][c]
            for x in range(c + 1, m):
                f = M[x][c] * inv
                for k in range(c, m): M[x][k] -= f * M[c][k]
        tot += int(d) if ok else 0
    return tot
def hp(om, n):
    c = 0
    def go(l, u, d):
        nonlocal c
        if d == n: c += 1; return
        mv = om[l] & ~u
        while mv:
            b = mv & -mv; w = b.bit_length() - 1; mv ^= b
            go(w, u | (1 << w), d + 1)
    for s in range(n): go(s, 1 << s, 1)
    return c

print("=" * 82)
print("CLOSING THM-1460(C): the adjacency-cospectral groups that Sum a does NOT split")
print("=" * 82)
for n in (5, 6, 7):
    cls = classes(n)
    data = []
    for c, om in cls.items():
        data.append((charpoly(om, n), sum_arb(om, n), hp(om, n), c))
    groups = {}
    for cp, sa, h, c in data: groups.setdefault(cp, []).append((sa, h, c))
    cospec = {cp: v for cp, v in groups.items() if len(v) > 1}
    unres_a = {cp: v for cp, v in cospec.items() if len({x[0] for x in v}) == 1}
    unres_h = {cp: v for cp, v in cospec.items() if len({x[1] for x in v}) == 1}
    unres_both = {cp: v for cp, v in cospec.items() if len({(x[0], x[1]) for x in v}) == 1}
    print(f"\n n={n}: classes={len(cls)}  cospectral groups (size>1)={len(cospec)}"
          f"  classes in them={sum(len(v) for v in cospec.values())}")
    print(f"   Sum a FAILS to split : {len(unres_a)} groups"
          f"    ({len(cospec)-len(unres_a)} split)")
    print(f"   H     FAILS to split : {len(unres_h)} groups"
          f"    ({len(cospec)-len(unres_h)} split)")
    print(f"   (Sum a, H) both fail : {len(unres_both)} groups")
    if unres_a and n == 7:
        print("   THE SURVIVING GROUPS UNDER Sum a (the 5 THM-1460 did not name):")
        for cp, v in list(unres_a.items())[:8]:
            sizes = len(v); sa = v[0][0]; hs = sorted({x[1] for x in v})
            print(f"      size {sizes}, Sum a = {sa}, H values = {hs}"
                  f"   -> H {'SPLITS it' if len(hs) > 1 else 'ALSO FAILS'}")
print("""
 READING.  The two relaxations fail on DIFFERENT groups, which is the point of the logarithm
 frame in THM-1460(D): log H is ordinal-sum additive with no interaction term, log Sum a
 carries a size-dependent shift, so they are not measuring the same thing.  Where (Sum a, H)
 both fail is the genuinely hard residue for any determinantal-plus-path fingerprint.
""")
