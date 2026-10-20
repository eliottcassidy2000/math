#!/usr/bin/env python3
"""THM-4724 Theorem R(U): for every S_n-submodule U, odd-order P, h in N(P), L = <P,h>:
     |Fix_{T/U}(L)| = sum_{X in (U^perp)^L} sgn_X(h)      (per-h identity)
Checked here for U = 0 (tournaments, U^perp = all graphs) and U = C (cut space; T/C = switching classes,
U^perp = Z = Euler graphs) at n = 5..8, for every odd-order subgroup class P and every h in N(P).
Left side by F_2 linear algebra on the affine action g.v = g(v) + D_g (D_g = reversal set of the transitive
reference tournament), right side by enumerating (U^perp)^L. For n = 6 the class counts are also brute-forced."""
import sys, itertools
from collections import Counter
NMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 5
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8

def cyc(perm_cycles, n):
    g = list(range(n))
    for c in perm_cycles:
        for i in range(len(c)): g[c[i]] = c[(i + 1) % len(c)]
    return tuple(g)
def mul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def inv(a):
    r = [0] * len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)
def closure(gens, n):
    idt = tuple(range(n)); S = {idt}; fr = [idt]
    while fr:
        x = fr.pop()
        for g in gens:
            y = mul(g, x)
            if y not in S: S.add(y); fr.append(y)
    return S
def order(g):
    idt = tuple(range(len(g))); x = g; k = 1
    while x != idt: x = mul(g, x); k += 1
    return k
def odd_partitions(n):
    def parts(n, m):
        if n == 0: yield (); return
        for k in range(min(n, m), 0, -1):
            if k % 2:
                for r in parts(n - k, k): yield (k,) + r
    return [p for p in parts(n, n) if p != tuple([1] * n)]
def perm_of_type(lam, n):
    cycles = []; v = 0
    for L in lam: cycles.append(list(range(v, v + L))); v += L
    return cyc(cycles, n)
def subgroup_reps(n):
    reps = {}
    for lam in odd_partitions(n):
        g = perm_of_type(lam, n); reps[f"Z{order(g)} type {lam}"] = [g]
    if n >= 6: reps["Z3xZ3 <(012),(345)>"] = [cyc([[0,1,2]], n), cyc([[3,4,5]], n)]
    if n >= 7: reps["Z7:Z3"] = [cyc([[0,1,2,3,4,5,6]], n), cyc([[1,2,4],[3,6,5]], n)]
    return reps

# ---- F_2 linear algebra on bitmasks ----
def echelon(vectors):
    """returns list of (pivot_bit, vector) in reduced echelon form"""
    basis = []
    for v in vectors:
        for pb, bv in basis:
            if v >> pb & 1: v ^= bv
        if v:
            pb = v.bit_length() - 1
            # reduce existing
            basis = [(p, (b ^ v) if (b >> pb & 1) else b) for p, b in basis]
            basis.append((pb, v))
    return basis
def reduce(v, basis):
    for pb, bv in basis:
        if v >> pb & 1: v ^= bv
    return v
def solve(rows, rhs, nvar):
    """rows: list of (bitmask over variables, rhs bit). returns (consistent, rank)."""
    piv = {}
    for r, b in zip(rows, rhs):
        for pb in sorted(piv.keys(), reverse=True):
            if r >> pb & 1: r ^= piv[pb][0]; b ^= piv[pb][1]
        if r == 0:
            if b: return False, len(piv)
            continue
        pb = r.bit_length() - 1
        # reduce others
        for q in list(piv.keys()):
            if piv[q][0] >> pb & 1: piv[q] = (piv[q][0] ^ r, piv[q][1] ^ b)
        piv[pb] = (r, b)
    return True, len(piv)
def kernel_basis(rows, nvar):
    """basis of {v : row.v = 0 for all rows} (rows as bitmasks over variables)"""
    piv = {}
    for r in rows:
        for pb in sorted(piv.keys(), reverse=True):
            if r >> pb & 1: r ^= piv[pb]
        if r:
            pb = r.bit_length() - 1
            for q in list(piv.keys()):
                if piv[q] >> pb & 1: piv[q] ^= r
            piv[pb] = r
    free = [k for k in range(nvar) if k not in piv]
    basis = []
    for f in free:
        v = 1 << f
        for pb, r in piv.items():
            if r >> f & 1: v |= 1 << pb
        basis.append(v)
    return basis

for n in range(NMIN, NMAX + 1):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {p: k for k, p in enumerate(pairs)}; m = len(pairs)
    perms = list(itertools.permutations(range(n)))
    def act(g, v):  # image of edge set v under g
        w = 0
        for k, (i, j) in enumerate(pairs):
            if v >> k & 1:
                a, b = g[i], g[j]; w |= 1 << idx[(min(a, b), max(a, b))]
        return w
    def inv_set(g):  # pairs {i<j} with g(i) > g(j)
        w = 0
        for k, (i, j) in enumerate(pairs):
            if g[i] > g[j]: w |= 1 << k
        return w
    def D(g):  # reversal set of T_0 under g: pairs {a<b} with g^{-1}(a) > g^{-1}(b)
        return inv_set(inv(g))
    stars = []
    for v in range(n):
        s = 0
        for k, (i, j) in enumerate(pairs):
            if v in (i, j): s |= 1 << k
        stars.append(s)
    Cbasis = echelon(stars)                    # cut space, dim n-1
    # incidence rows for Euler condition: for vertex v, the row is star(v); Z = kernel
    def count_fixed_classes(gens, quotient):
        """number of L-fixed points of the torsor T/U (U = C if quotient else 0) via the affine system"""
        rows = []; rhs = []
        for g in gens:
            Dg = D(g)
            if quotient: Dg = reduce(Dg, Cbasis)
            # linear map v -> pi((g+1)v): columns
            cols = []
            for k in range(m):
                w = act(g, 1 << k) ^ (1 << k)
                if quotient: w = reduce(w, Cbasis)
                cols.append(w)
            # rows of the matrix: for each coordinate bit t of the codomain
            for t in range(m):
                r = 0
                for k in range(m):
                    if cols[k] >> t & 1: r |= 1 << k
                rows.append(r); rhs.append(Dg >> t & 1)
        ok, rank = solve(rows, rhs, m)
        if not ok: return 0
        sols = 2 ** (m - rank)
        return sols // (2 ** (n - 1)) if quotient else sols
    def signed_sum(gens, h, euler):
        rows = []
        for g in gens:
            for k in range(m):
                pass
        # constraints: (g+1)X = 0 for generators; Euler: star rows
        cons = []
        for g in gens:
            # matrix of g+1: row t gets bit k if (act(g,e_k)^e_k) has bit t
            cols = [act(g, 1 << k) ^ (1 << k) for k in range(m)]
            for t in range(m):
                r = 0
                for k in range(m):
                    if cols[k] >> t & 1: r |= 1 << k
                cons.append(r)
        if euler: cons += stars
        B = kernel_basis(cons, m)
        Ih = inv_set(h); tot = 0
        for bits in range(2 ** len(B)):
            X = 0
            for i, b in enumerate(B):
                if bits >> i & 1: X ^= b
            tot += -1 if bin(X & Ih).count('1') % 2 else 1
        return tot, len(B)
    reps = subgroup_reps(n)
    stats = Counter(); bad = []
    for name, gens in reps.items():
        P = closure(gens, n)
        N = [s for s in perms if frozenset(mul(mul(s, g), inv(s)) for g in P) == P]
        seen_L = set()
        for h in N:
            L = frozenset(closure(gens + [h], n))
            key = (L, h)
            lhs0 = count_fixed_classes(gens + [h], False); rhs0, d0 = signed_sum(gens + [h], h, False)
            lhs1 = count_fixed_classes(gens + [h], True);  rhs1, d1 = signed_sum(gens + [h], h, True)
            stats['pairs'] += 1
            if lhs0 != rhs0: bad.append((n, name, 'U=0', order(h), lhs0, rhs0))
            if lhs1 != rhs1: bad.append((n, name, 'U=C', order(h), lhs1, rhs1))
            if order(h) % 2 == 0: stats['even h'] += 1
            if lhs1 == 0: stats['U=C zero'] += 1
    print(f"n = {n}: {stats['pairs']} pairs (P,h) over {len(reps)} odd subgroup classes, {stats['even h']} with h of even order; "
          f"U=C identity zero on {stats['U=C zero']} pairs; mismatches: {len(bad)}")
    for b in bad[:10]: print("   MISMATCH", b)
    if n == 6:
        # brute force the switching-class fixed-point counts for one P and all h
        name, gens = "Z3 type (3, 1, 1, 1)", [cyc([[0,1,2]], 6)]
        P = closure(gens, 6); N = [s for s in perms if frozenset(mul(mul(s, g), inv(s)) for g in P) == P]
        # classes: tournaments v in F_2^15 modulo C
        Cset = set()
        for bits in range(2 ** len(Cbasis)):
            c = 0
            for i, (pb, bv) in enumerate(Cbasis):
                if bits >> i & 1: c ^= bv
            Cset.add(c)
        def cls(v): return min(v ^ c for c in Cset)
        checks = 0
        for h in N:
            gl = gens + [h]
            fixed = 0
            for v in range(2 ** m):
                if v != cls(v): continue
                if all(cls(act(g, v) ^ D(g)) == v for g in gl): fixed += 1
            assert fixed == count_fixed_classes(gl, True), (h, fixed)
            checks += 1
        print(f"   n = 6 brute-force switching-class fixed-point counts agree with the linear-algebra counts for {checks} elements h")
