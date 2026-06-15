#!/usr/bin/env python3
"""
permanental_companion_monad.py  (monad-explorer-2026-06-15-S6)  [pure-python, no numpy]

THE PERMANENTAL POLYNOMIAL = the UNSIGNED TWIN of the characteristic polynomial,
and the canonical "all face" of the master cycle-packing polynomial Phi (OPEN-Q-096.1,
reflection `the-master-cycle-packing-polynomial`, THM-505).

BACKGROUND.  For a tournament T on n vertices with 0/1 adjacency A, a *linear subdigraph*
L = a set of vertex-disjoint directed cycles (length >= 3 in a tournament).  The master
polynomial is Phi(T;{y_k}) = sum_L prod_{C in L} y_{|C|}.  Two classical matrix functions
are the all-length, VERTEX-graded faces of Phi, differing ONLY by the cycle-parity sign:

  det(xI - A) = sum_L (-1)^{#cyc(L)} x^{n-|V(L)|}        (SIGNED  -> spectrum / Sachs)
  per(xI + A) = sum_L (+1)^{#cyc(L)} x^{n-|V(L)|}        (UNSIGNED -> permanental poly)

So  [coeff of x^{n-m} in det(xI-A)] = e_m^signed   = sum_{|V(L)|=m} (-1)^{#cyc}   (char poly)
    [coeff of x^{n-m} in per(xI+A)] = e_m^unsigned = sum_{|V(L)|=m} (+1)^{#cyc}   (perm poly)

The determinant is in P (eigenvalues); the permanent is #P-hard (Valiant).  The permanental
polynomial is the "spectrum with the cycle-parity signs stripped" -- it carries exactly the
non-spectral parity information.  This script:

 [1] verifies the PERMANENTAL COEFFICIENT THEOREM   per(xI+A)_coeff == e_m^unsigned
     (Ryser permanent over the formal variable x, vs independent packing enumeration);
 [2] verifies the det/per O/E decomposition: for each m,
        O_m := #{odd-cardinality packings on m vtx} = (e_m^unsigned - e_m^signed)/2
        E_m := #{even-cardinality packings on m vtx} = (e_m^unsigned + e_m^signed)/2 ;
 [3] cospectral probe: does the perm poly SPLIT cospectral classes? at which first n?
 [4] FINGERPRINT completeness (exhaustive n<=6): #distinct char-poly vs #distinct
     (char,perm) vs #distinct (char,perm,H) vs A000568 (#tournaments up to iso);
 [5] does the PAIR (char poly, perm poly) DETERMINE H?  refine cospectral classes by perm
     poly and check whether H is constant -- find the first n where it fails + a witness.

Char poly via Faddeev-LeVerrier (exact int). Perm poly via Ryser over Z[x] (exact int).
"""
import sys, itertools, random
from fractions import Fraction

# ----------------------------------------------------------------------------- tournaments
def random_tournament(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.randint(0, 1): A[i][j] = 1
            else:                 A[j][i] = 1
    return A

def all_labeled_tournaments(n):
    """Generator over all 2^C(n,2) labeled tournaments (n small)."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(1 << len(pairs)):
        A = [[0]*n for _ in range(n)]
        for t, (i, j) in enumerate(pairs):
            if (bits >> t) & 1: A[i][j] = 1
            else:               A[j][i] = 1
        yield A

# ----------------------------------------------------------------------------- char poly
def matmul(A, B, n):
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        Ai = A[i]; Ci = C[i]
        for k in range(n):
            a = Ai[k]
            if a:
                Bk = B[k]
                for j in range(n): Ci[j] += a*Bk[j]
    return C

def charpoly_int(A, n):
    """det(xI - A) = x^n + c1 x^{n-1} + ... ; returns [1,c1,...,cn]; coeff of x^{n-m}=out[m]."""
    M = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    coeffs = [1]
    for k in range(1, n+1):
        AM = matmul(A, M, n)
        tr = sum(AM[i][i] for i in range(n))
        ck = Fraction(-tr, k); assert ck.denominator == 1
        ck = ck.numerator; coeffs.append(ck)
        if k < n:
            M = [[AM[i][j] + (ck if i==j else 0) for j in range(n)] for i in range(n)]
    return coeffs   # e_m^signed = coeffs[m]

# ----------------------------------------------------------------------------- perm poly (Ryser over Z[x])
def polymul(p, q):
    r = [0]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q): r[i+j] += a*b
    return r

def permpoly_int(A, n):
    """per(xI + A) as a list of coeffs [k0,k1,...,kn] with k_t = coeff of x^t.
    Returns the polynomial reindexed so out[m] = coeff of x^{n-m} = e_m^unsigned.
    Ryser:  per(M) = sum_{S subset [n]} (-1)^{n-|S|} prod_i (sum_{j in S} M[i][j]).
    With M = xI + A, row-sum over S is the degree-<=1 poly  ([i in S])*x + a_i(S),
    a_i(S) = #{j in S : A[i][j]=1}.  Product over i has degree |S|."""
    full = [0]*(n+1)   # full[t] = coeff of x^t in per(xI+A)
    rng_n = range(n)
    for Sbits in range(1 << n):
        S = [j for j in rng_n if (Sbits >> j) & 1]
        Sset = set(S)
        prod = [1]
        for i in rng_n:
            ai = sum(1 for j in S if A[i][j])
            term = ([ai, 1] if i in Sset else [ai])   # x*[i in S] + a_i(S)
            prod = polymul(prod, term)
        sign = (-1)**(n - len(S))
        for t, c in enumerate(prod):
            full[t] += sign*c
    # reindex: out[m] = coeff of x^{n-m}
    return [full[n-m] for m in range(n+1)]   # e_m^unsigned = out[m]

# ----------------------------------------------------------------------------- packings (ground truth)
def all_cycles(A, n):
    cycles = []
    for start in range(n):
        path = [start]; visited = {start}
        def dfs(u):
            for w in range(n):
                if w == start:
                    if len(path) >= 3 and A[u][start]:
                        cycles.append((len(path), frozenset(path)))
                elif w > start and w not in visited and A[u][w]:
                    visited.add(w); path.append(w); dfs(w); path.pop(); visited.discard(w)
        dfs(start)
    return cycles

def packing_OE(A, n):
    """From direct packing enumeration: e_signed[m], e_unsigned[m], O[m], E[m]."""
    cycles = all_cycles(A, n)
    vsets = [vs for (_, vs) in cycles]
    e_signed = {}; e_unsigned = {}; O = {}; E = {}
    n_c = len(cycles)
    def rec(start, used, nc, cov):
        m = cov
        e_signed[m]   = e_signed.get(m, 0)   + (-1)**nc
        e_unsigned[m] = e_unsigned.get(m, 0) + 1
        if nc % 2 == 1: O[m] = O.get(m, 0) + 1
        else:           E[m] = E.get(m, 0) + 1
        for i in range(start, n_c):
            if not (vsets[i] & used):
                rec(i+1, used | vsets[i], nc+1, cov + len(vsets[i]))
    rec(0, frozenset(), 0, 0)
    return e_signed, e_unsigned, O, E

def count_ham_paths(A, n):
    """#directed Hamiltonian paths via bitmask DP: O(2^n n^2).  dp[mask][v] = #paths
    covering vertex set `mask` ending at v."""
    full = (1 << n) - 1
    # dp indexed [mask][v]
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        if not any(row): continue
        for v in range(n):
            cv = row[v]
            if not cv: continue
            Av = A[v]
            for w in range(n):
                if not (mask >> w) & 1 and Av[w]:
                    dp[mask | (1 << w)][w] += cv
    return sum(dp[full][v] for v in range(n))

# ----------------------------------------------------------------------------- [1][2] verify theorems
def verify_perm_coeff_theorem(n, n_samples, seed=11):
    rng = random.Random(seed); fail = 0
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        pp = permpoly_int(A, n)                 # e_m^unsigned via Ryser
        es, eu, O, E = packing_OE(A, n)         # via packing enumeration
        for m in range(n+1):
            packing_u = 1 if m == 0 else eu.get(m, 0)
            if pp[m] != packing_u:
                fail += 1
                if fail <= 3: print(f"   [n={n}] PERM-COEFF FAIL m={m}: ryser={pp[m]} packings={packing_u}")
                break
    print(f" [1] permanental coeff thm  per(xI+A)_coeff == e_m^unsigned  n={n}: {n_samples-fail}/{n_samples}")
    return fail == 0

def verify_OE(n, n_samples, seed=12):
    rng = random.Random(seed); fail = 0
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        cp = charpoly_int(A, n)                 # e_m^signed
        pp = permpoly_int(A, n)                 # e_m^unsigned
        es, eu, O, E = packing_OE(A, n)
        for m in range(n+1):
            esig = cp[m]; euns = pp[m]
            if (euns - esig) % 2 or (euns + esig) % 2: fail += 1; break
            Om = (euns - esig)//2; Em = (euns + esig)//2
            Otrue = (1 if m==0 else 0) if m==0 else O.get(m,0)
            Etrue = (1 if m==0 else E.get(m,0))   # empty packing (m=0) is EVEN cardinality (0 cycles)
            if m == 0: Otrue, Etrue = 0, 1
            if Om != Otrue or Em != Etrue:
                fail += 1
                if fail <= 3: print(f"   [n={n}] O/E FAIL m={m}: O={Om}/{Otrue} E={Em}/{Etrue}")
                break
    print(f" [2] det/per O/E split   O_m=(eu-es)/2, E_m=(eu+es)/2 = #odd/even-card packings  n={n}: {n_samples-fail}/{n_samples}")
    return fail == 0

# ----------------------------------------------------------------------------- [3] cospectral probe
def cospectral_perm_probe(n, n_samples, seed=7):
    rng = random.Random(seed); buckets = {}
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        cp = tuple(charpoly_int(A, n))
        pp = tuple(permpoly_int(A, n))
        H  = count_ham_paths(A, n)
        b = buckets.setdefault(cp, {'perm': set(), 'H': set()})
        b['perm'].add(pp); b['H'].add(H)
    perm_split = sum(1 for cp,d in buckets.items() if len(d['perm']) > 1)
    H_split    = sum(1 for cp,d in buckets.items() if len(d['H'])    > 1)
    print(f" [3] cospectral n={n}: {len(buckets)} classes / {n_samples} samples | "
          f"perm-poly splits {perm_split:3d}  |  H splits {H_split:3d}")
    return perm_split, H_split, len(buckets)

# ----------------------------------------------------------------------------- [4] fingerprint completeness (exhaustive)
A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}  # #tournaments up to iso

def fingerprint_completeness(n):
    char_set = set(); cp_set = set(); cph_set = set()
    for A in all_labeled_tournaments(n):
        cp = tuple(charpoly_int(A, n))
        pp = tuple(permpoly_int(A, n))
        H  = count_ham_paths(A, n)
        char_set.add(cp); cp_set.add((cp, pp)); cph_set.add((cp, pp, H))
    iso = A000568[n]
    print(f" [4] n={n}: iso-classes(A000568)={iso:5d} | distinct char-poly={len(char_set):5d} "
          f"| distinct (char,perm)={len(cp_set):5d} | distinct (char,perm,H)={len(cph_set):5d}")
    return len(char_set), len(cp_set), len(cph_set), iso

# ----------------------------------------------------------------------------- [5] does (char,perm) determine H?
def det_perm_determines_H(n, n_samples, seed=23):
    """Refine cospectral classes by perm poly; is H constant on each refined class?"""
    rng = random.Random(seed); refined = {}
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        key = (tuple(charpoly_int(A, n)), tuple(permpoly_int(A, n)))
        H = count_ham_paths(A, n)
        refined.setdefault(key, {}).setdefault(H, []).append(A)
    bad = [(k, d) for k, d in refined.items() if len(d) > 1]
    print(f" [5] (char,perm)->H  n={n}: {len(refined)} (char,perm) classes / {n_samples} samples | "
          f"{len(bad)} classes with SPLIT H  -> "
          f"{'DETERMINES H' if not bad else 'does NOT determine H'}")
    if bad:
        k, d = bad[0]
        Hs = sorted(d.keys())
        print(f"     WITNESS: a (char,perm) class carries H in {Hs}")
        # also report which carriers differ for the witness pair
        A1 = d[Hs[0]][0]; A2 = d[Hs[1]][0]
        es1,eu1,O1,E1 = packing_OE(A1, n); es2,eu2,O2,E2 = packing_OE(A2, n)
        diffsO = {m:(O1.get(m,0),O2.get(m,0)) for m in set(O1)|set(O2) if O1.get(m,0)!=O2.get(m,0)}
        diffsE = {m:(E1.get(m,0),E2.get(m,0)) for m in set(E1)|set(E2) if E1.get(m,0)!=E2.get(m,0)}
        print(f"     (same O_m,E_m by construction; the differing FINE carriers live inside a parity class)")
    return len(bad) == 0

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    print("="*82)
    print(" PERMANENTAL POLYNOMIAL = the UNSIGNED TWIN of the char poly (all-face of Phi)")
    print("="*82)

    print("\n[1] PERMANENTAL COEFFICIENT THEOREM (Ryser per(xI+A) vs packing enumeration)")
    for n in (3,4,5,6,7): verify_perm_coeff_theorem(n, 120 if n < 7 else 60)

    print("\n[2] det/per O/E DECOMPOSITION (parity-of-cardinality packing counts)")
    for n in (3,4,5,6,7): verify_OE(n, 120 if n < 7 else 60)

    print("\n[3] COSPECTRAL PROBE: perm poly vs H, which splits cospectral classes & when")
    for n in (4,5,6,7,8): cospectral_perm_probe(n, 4000 if n <= 6 else 8000)

    print("\n[4] FINGERPRINT COMPLETENESS (exhaustive; vs A000568 = #tournaments up to iso)")
    for n in (3,4,5,6): fingerprint_completeness(n)

    print("\n[5] DOES (char poly, perm poly) DETERMINE H?  (refine cospectral by perm; first fail)")
    for n in (6,7,8,9): det_perm_determines_H(n, 5000 if n <= 7 else 9000)
    print("\nDONE.")
