# ocf_neg_symfunc_kpo2_verify.py
# ADVERSARIAL VERIFIER, thread C (HYP-2380): claims C2 (omega-invariance of U_T)
# and the computationally checkable parts of C3 (GS Thm 1.39 expansion, Lemma 6.5
# zeta evaluation, the p_1->1, p_odd->-1 specialization giving I(Omega,-2)).
# FRESH code; symmetric functions handled as honest polynomials in m=n variables
# (degree n in n variables determines a symmetric function of degree n uniquely).
#
# Definitions used (to be cross-checked against the GS paper text):
#   U_D = sum over listings w (orderings of V) of L_{Des(w,D),n}
#   Des(w,D) = { i in [n-1] : (w_i, w_{i+1}) is an arc of D }   [convention A]
#       (convention B = complement is ALSO tested; report which one works)
#   L_{S,n} = Gessel fundamental quasisymmetric: sum over 1<=i_1<=...<=i_n<=m
#       with i_j < i_{j+1} whenever j in S.
# Checks per tournament (all 8 at n=3, all 64 at n=4, 5 random at n=5):
#   (1) U_D == sum over collections K of disjoint directed odd cycles of
#       2^{|K|} p_{lambda(K)}  where lambda(K) = cycle lengths + 1^{#uncovered}
#       (GS Thm 1.39 for tournaments)  [under which Des convention?]
#   (2) expansion in p-basis (exact Fraction solve): coefficients supported on
#       all-odd partitions => omega(U)=U since omega(p_lambda)=(-1)^{n-len} p_lambda
#       and n-len is even for all-odd lambda; ALSO verified directly.
#   (3) zeta: p_k -> 1 gives H; principal specialization x=(1,0,0,..) of U gives
#       #{w : Des(w,D) = emptyset} = ham paths of complement.
#   (4) p_1 -> 1, p_{2k+1} -> -1: equals I(Omega,-2) from independent census code.

import itertools, random
from fractions import Fraction
from collections import Counter

random.seed(31337)

def pairs_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_mask(n, mask, pairs):
    adj = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(pairs):
        if (mask >> b) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def odd_cycles_with_support(n, adj):
    res = []
    for k in range(3, n + 1, 2):
        for sub in itertools.combinations(range(n), k):
            v0 = sub[0]
            for perm in itertools.permutations(sub[1:]):
                cyc = (v0,) + perm
                if all(adj[cyc[i]][cyc[(i + 1) % k]] for i in range(k)):
                    m = 0
                    for v in sub:
                        m |= 1 << v
                    res.append((m, k))
    return res

# polynomials in m variables: dict exponent-tuple -> int/Fraction
def pmul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, 0) + ca * cb
    return out

def padd(a, b):
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, 0) + c
        if out[e] == 0:
            del out[e]
    return out

def p_k(k, m):
    out = {}
    for i in range(m):
        e = [0] * m
        e[i] = k
        out[tuple(e)] = 1
    return out

def p_lambda(lam, m):
    out = {tuple([0] * m): 1}
    for k in lam:
        out = pmul(out, p_k(k, m))
    return out

def fundamental_L(S, n, m):
    """L_{S,n} in m variables: sum over 1<=i_1<=..<=i_n<=m, strict at positions in S."""
    out = {}
    def rec(pos, lastvar, e):
        if pos == n:
            t = tuple(e)
            out[t] = out.get(t, 0) + 1
            return
        start = lastvar + 1 if pos > 0 and (pos in S) else max(lastvar, 0)
        if pos == 0:
            start = 0
        for v in range(start, m):
            e[v] += 1
            rec(pos + 1, v, e)
            e[v] -= 1
    rec(0, -1, [0] * m)
    return out

def partitions(n):
    if n == 0:
        yield ()
        return
    def rec(rem, maxp):
        if rem == 0:
            yield ()
            return
        for p in range(min(rem, maxp), 0, -1):
            for rest in rec(rem - p, p):
                yield (p,) + rest
    yield from rec(n, n)

def p_basis_decompose(U, n, m):
    """solve U = sum c_lambda p_lambda exactly; returns dict lambda -> Fraction."""
    lams = list(partitions(n))
    # pick monomial representatives: exponent = partition padded with zeros
    monos = [tuple(list(lam) + [0] * (m - len(lam))) for lam in lams]
    plams = {lam: p_lambda(lam, m) for lam in lams}
    A = [[Fraction(plams[lam].get(mono, 0)) for lam in lams] for mono in monos]
    bvec = [Fraction(U.get(mono, 0)) for mono in monos]
    # gaussian elimination
    k = len(lams)
    for col in range(k):
        piv = next(r for r in range(col, k) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        bvec[col], bvec[piv] = bvec[piv], bvec[col]
        inv = 1 / A[col][col]
        A[col] = [x * inv for x in A[col]]
        bvec[col] *= inv
        for r in range(k):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
                bvec[r] -= f * bvec[col]
    c = {lam: bvec[i] for i, lam in enumerate(lams)}
    # verify reconstruction exactly (also catches non-symmetric U)
    recon = {}
    for lam, coef in c.items():
        if coef != 0:
            recon = padd(recon, {e: coef * v for e, v in plams[lam].items()})
    Uf = {e: Fraction(v) for e, v in U.items() if v != 0}
    assert recon == Uf, "p-basis reconstruction failed (U not symmetric?)"
    return c

def ham_perm(n, adj):
    return sum(1 for p in itertools.permutations(range(n))
               if all(adj[p[i]][p[i + 1]] for i in range(n - 1)))

results = Counter()
conv_match = Counter()

def test_tournament(n, mask):
    m = n
    pairs = pairs_list(n)
    adj = adj_from_mask(n, mask, pairs)
    cycles = odd_cycles_with_support(n, adj)
    # claimed expansion: sum over disjoint collections of 2^{|K|} p_{type}
    full = (1 << n) - 1
    expans = {}
    def rec(idx, used, lens):
        # contribute current collection
        lam = tuple(sorted(lens + [1] * (n - bin(used).count('1')), reverse=True))
        key = (lam, len(lens))
        expans[key] = expans.get(key, 0) + 1
        for i in range(idx, len(cycles)):
            sm, k = cycles[i]
            if sm & used == 0:
                rec(i + 1, used | sm, lens + [k])
    rec(0, 0, [])
    claimed = {}
    for (lam, sz), cnt in expans.items():
        poly = p_lambda(lam, m)
        claimed = padd(claimed, {e: cnt * (2 ** sz) * v for e, v in poly.items()})
    # U_D under both Des conventions
    for conv in ('A:in-D', 'B:not-in-D'):
        U = {}
        for w in itertools.permutations(range(n)):
            if conv == 'A:in-D':
                S = {i + 1 for i in range(n - 1) if adj[w[i]][w[i + 1]]}
            else:
                S = {i + 1 for i in range(n - 1) if not adj[w[i]][w[i + 1]]}
            U = padd(U, fundamental_L(S, n, m))
        if U == claimed:
            conv_match[conv] += 1
            # decompose and run checks
            c = p_basis_decompose(U, n, m)
            allodd = all(all(part % 2 == 1 for part in lam)
                         for lam, coef in c.items() if coef != 0)
            results['thm139_allodd' if allodd else 'thm139_allodd_FAIL'] += 1
            # omega directly
            omega_ok = all(coef == 0 or (n - len(lam)) % 2 == 0
                           for lam, coef in c.items())
            results['omega_invariant' if omega_ok else 'omega_FAIL'] += 1
            # integrality + 2^psi structure: coef of all-odd lambda with t parts>=3
            # should be a nonneg integer divisible by 2^t? (Cor 1.40: N[p1,2p3,...])
            cor140 = True
            for lam, coef in c.items():
                if coef == 0:
                    continue
                t = sum(1 for part in lam if part >= 3)
                if coef.denominator != 1 or coef < 0 or int(coef) % (2 ** t) != 0:
                    cor140 = False
            results['cor140' if cor140 else 'cor140_FAIL'] += 1
            # zeta: p_k -> 1
            zeta = sum(coef for coef in c.values())
            H = ham_perm(n, adj)
            results['zeta==H' if zeta == H else 'zeta==H_FAIL'] += 1
            # principal specialization x=(1,0,...): count monomials e=(n,0,..0)
            e1 = tuple([n] + [0] * (m - 1))
            ps = U.get(e1, 0)
            # = #listings with empty descent set = ham paths of complement digraph
            results['ps1==H' if ps == H else 'ps1==H_FAIL'] += 1
            # negative specialization p_1->1, p_{2k+1}->-1
            spec = 0
            for lam, coef in c.items():
                if coef == 0:
                    continue
                t = sum(1 for part in lam if part >= 3)
                spec += coef * ((-1) ** t)
            a1 = len(cycles)
            a2 = sum(1 for i in range(len(cycles)) for j in range(i + 1, len(cycles))
                     if cycles[i][0] & cycles[j][0] == 0)
            # n<=5 here so alpha_3=0
            Im2 = 1 - 2 * a1 + 4 * a2
            results['negspec==I(-2)' if spec == Im2 else 'negspec==I(-2)_FAIL'] += 1
        else:
            conv_match[conv + '_NO'] += 1

print("symmetric-function verification of GS Thm 1.39 / omega-invariance (fresh code)")
for n in (3, 4):
    pairs = pairs_list(n)
    for mask in range(1 << len(pairs)):
        test_tournament(n, mask)
    print(f"n={n}: done ({1 << len(pairs)} tournaments x 2 Des conventions)")
for _ in range(5):
    test_tournament(5, random.getrandbits(10))
print("n=5: done (5 random)")

print("\nDes convention match counts (which convention satisfies Thm 1.39 expansion):")
for k, v in sorted(conv_match.items()):
    print(f"  {k}: {v}")
print("\ncheck tallies:")
for k, v in sorted(results.items()):
    print(f"  {k}: {v}")
nf = sum(v for k, v in results.items() if k.endswith('_FAIL'))
print(f"\nTOTAL FAILURES: {nf}")
print("done")
