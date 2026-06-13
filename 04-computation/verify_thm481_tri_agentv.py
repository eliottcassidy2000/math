#!/usr/bin/env python3
"""
verify_thm481_tri_agentv.py — adversarial independent verification (agent "tri_agentv", 2026-06-11)

Verifies, with code written FRESH from definitions (no imports from repo scripts):

(1) THE ORDER-32 PALEY CLAIM (THM-481 part C / THM-480 contrast claim):
    H = bordered Paley_31, indices {inf} u Z_31, with
        H[inf,inf] = 1, H[inf,j] = +1, H[i,inf] = -1, H[i,i] = +1, H[i,j] = chi(j-i) mod 31
    (chi = quadratic-residue character; p = 31 == 3 mod 4 so chi(-x) = -chi(x), H skew-type).
    Claim: C(H) = F2-row-span of (J-H)/2 is a [32,16] self-dual doubly-even (Type II) code
    with weight distribution
        A0=1, A8=620, A12=13888, A16=36518, A20=13888, A24=620, A32=1   (A4 = 0).

(3) THE TOWER CONTRAST (THM-480 consistency):
    Tower: H_8 = border(Paley_7) (same gauge), H_{2m} = [[H, H], [-H^T, H^T]] (THM-447 doubling).
    Claim: C(H_32^tower) is [32,16] Type II with the d32+ enumerator
        A0=1, A4=120, A8=1820, A12=8008, A16=45638, A20=8008, A24=1820, A28=120, A32=1.
    Different enumerators ==> the two CODES are inequivalent. To upgrade to HADAMARD
    inequivalence of the two order-32 +-1 matrices we additionally compute the 4-profile
    (multiset of |sum_j H[a,j]H[b,j]H[c,j]H[d,j]| over all 4-subsets of rows), which IS a
    full Hadamard-equivalence invariant (invariant under row/col permutations AND row/col
    negations), plus the same on columns to cover transposition.

GAUGE CAVEAT (checked empirically below): the binary code C(H) itself is NOT a priori a
Hadamard-equivalence invariant. Row permutations fix C; column permutations send C to an
equivalent code; row negations fix C provided the all-ones vector handling works out
(we verify directly: complementing any single row leaves the span unchanged for both
matrices, which holds because 1 is in C and C is spanned without any single row);
column negations send the generators b_i -> b_i + e_j and CAN change the code (we test
all 32 column negations for both matrices and report whether the span moves).
Hence: enumerator difference alone certifies inequivalence only under the subgroup
<row perms, col perms, row negations>; the 4-profile closes the column-negation gap.

Everything is exact integer / F2 arithmetic. Runtime ~ a few seconds.
"""

import itertools
import sys
from collections import Counter

# ----------------------------------------------------------------------------
# basic constructions
# ----------------------------------------------------------------------------

def chi_fn(p):
    qr = {pow(x, 2, p) for x in range(1, p)}
    def chi(a):
        a %= p
        if a == 0:
            return 0
        return 1 if a in qr else -1
    return chi


def bordered_paley(p):
    """Bordered Paley matrix, order p+1, index 0 = inf, index i+1 = residue i.
    H[inf,*]=+1 (incl. diagonal), H[*,inf]=-1, H[i,i]=+1, H[i,j]=chi(j-i)."""
    chi = chi_fn(p)
    n = p + 1
    H = [[0] * n for _ in range(n)]
    H[0][0] = 1
    for j in range(1, n):
        H[0][j] = 1
        H[j][0] = -1
    for i in range(1, n):
        for j in range(1, n):
            H[i][j] = 1 if i == j else chi((j - 1) - (i - 1))
    return H


def double_skew(H):
    """THM-447 doubling: H_{2m} = [[H, H], [-H^T, H^T]]."""
    m = len(H)
    HT = [[H[j][i] for j in range(m)] for i in range(m)]
    out = []
    for i in range(m):
        out.append(H[i] + H[i])
    for i in range(m):
        out.append([-x for x in HT[i]] + HT[i])
    return out


def check_hadamard(H):
    n = len(H)
    assert all(abs(x) == 1 for row in H for x in row), "entries not +-1"
    for i in range(n):
        for k in range(i, n):
            d = sum(H[i][j] * H[k][j] for j in range(n))
            expect = n if i == k else 0
            if d != expect:
                return False
    return True


def check_skew_type(H):
    n = len(H)
    for i in range(n):
        for j in range(n):
            s = H[i][j] + H[j][i]
            if (i == j and s != 2) or (i != j and s != 0):
                return False
    return True

# ----------------------------------------------------------------------------
# F2 linear algebra on bitmask ints (bit j = coordinate j)
# ----------------------------------------------------------------------------

def binarize_rows(H):
    """rows of (J - H)/2 as bitmasks: bit j set iff H[i][j] == -1."""
    n = len(H)
    rows = []
    for i in range(n):
        v = 0
        for j in range(n):
            if H[i][j] == -1:
                v |= 1 << j
        rows.append(v)
    return rows


def rref_basis(vectors):
    """Canonical reduced-echelon basis (sorted tuple is a canonical form of the span).
    Forward elimination on leading bits, then full back-substitution in increasing
    pivot order (correct: once bit p is cleared from all higher rows it never returns)."""
    piv = {}  # pivot bit -> vector
    for v in vectors:
        while v:
            p = v.bit_length() - 1
            if p in piv:
                v ^= piv[p]
            else:
                piv[p] = v
                break
    for p in sorted(piv):
        for q in piv:
            if q > p and (piv[q] >> p) & 1:
                piv[q] ^= piv[p]
    return tuple(sorted(piv.values()))


def span_words(basis):
    words = [0]
    for b in basis:
        words += [w ^ b for w in words]
    return words


def weight_distribution(basis):
    return Counter(bin(w).count("1") for w in span_words(basis))


def is_self_orthogonal(basis):
    return all(bin(a & b).count("1") % 2 == 0
               for i, a in enumerate(basis) for b in basis[i:])

# ----------------------------------------------------------------------------
# Hadamard-equivalence invariant: the 4-profile
# ----------------------------------------------------------------------------

def four_profile(H, on_columns=False):
    """multiset of |sum_j prod_{r in Q} H[r,j]| over 4-subsets Q of rows (or columns).
    Invariant under H -> P H Q for signed permutation matrices P, Q:
      * column negation flips 4 factors of each product (no change);
      * row negation flips the sign of the whole inner sum (killed by |.|);
      * permutations permute the multiset.  (Symmetrically for the column profile.)"""
    n = len(H)
    M = [[H[j][i] for j in range(n)] for i in range(n)] if on_columns else H
    masks = binarize_rows(M)  # product over Q at col j is -1 iff parity of -1s odd
    prof = Counter()
    for a, b, c, d in itertools.combinations(range(n), 4):
        x = masks[a] ^ masks[b] ^ masks[c] ^ masks[d]
        prof[abs(n - 2 * bin(x).count("1"))] += 1
    return prof

# ----------------------------------------------------------------------------
# gauge-orbit probes (row / column negation effect on the code)
# ----------------------------------------------------------------------------

def gauge_probes(rows, n, canon):
    full = (1 << n) - 1
    row_neg_fixed = all(
        rref_basis([r ^ full if i == k else r for k, r in enumerate(rows)]) == canon
        for i in range(n))
    col_changed, col_dims = [], set()
    for j in range(n):
        c2 = rref_basis([r ^ (1 << j) for r in rows])
        col_dims.add(len(c2))
        if c2 != canon:
            col_changed.append(j)
    return row_neg_fixed, col_changed, sorted(col_dims)

# ----------------------------------------------------------------------------
# expected enumerators
# ----------------------------------------------------------------------------

CLAIM_PALEY32 = {0: 1, 8: 620, 12: 13888, 16: 36518, 20: 13888, 24: 620, 32: 1}
CLAIM_TOWER32 = {0: 1, 4: 120, 8: 1820, 12: 8008, 16: 45638, 20: 8008, 24: 1820, 28: 120, 32: 1}
# d32+ closed form: pair-doubled words on 16 coordinate pairs — choosing j pairs gives
# weight 2j, doubly-even forces j even, so A_{2j} = C(16, j) for even j (dim-15 "d" part);
# the glue coset adds 2^15 words, all of weight 16.
D32PLUS = {2 * j: __import__("math").comb(16, j) for j in range(0, 17, 2)}
D32PLUS[16] += 2 ** 15

CLAIM_TOWER16 = {0: 1, 4: 28, 8: 198, 12: 28, 16: 1}
CLAIM_E8 = {0: 1, 4: 14, 8: 1}


def analyze(name, H, claimed=None):
    n = len(H)
    print(f"\n=== {name} (order {n}) ===")
    had = check_hadamard(H)
    skw = check_skew_type(H)
    print(f"  Hadamard (H H^T = {n} I): {had}")
    print(f"  skew-type (H + H^T = 2I): {skw}")
    rows = binarize_rows(H)
    basis = rref_basis(rows)
    canon = basis
    dim = len(basis)
    print(f"  dim C(H) = rank_F2 (J-H)/2 = {dim}  (n/2 = {n // 2})")
    so = is_self_orthogonal(list(basis))
    wd = weight_distribution(basis)
    de = all(w % 4 == 0 for w in wd)
    selfdual = so and dim == n // 2
    print(f"  self-orthogonal: {so}; self-dual (dim=n/2 & C<=C^perp): {selfdual}")
    print(f"  doubly-even (all weights = 0 mod 4): {de}  => Type II: {selfdual and de}")
    allones = (1 << n) - 1
    has1 = rref_basis(list(basis) + [allones]) == canon
    print(f"  all-ones vector in C: {has1}")
    dist = sorted(wd.items())
    print(f"  weight distribution: {dist}")
    nonzero = [w for w in sorted(wd) if w > 0]
    print(f"  minimum distance: {nonzero[0]}")
    if claimed is not None:
        ok = dict(wd) == claimed
        print(f"  MATCHES claimed enumerator: {ok}")
        if not ok:
            print(f"    claimed: {sorted(claimed.items())}")
    return {"H": H, "rows": rows, "basis": basis, "canon": canon, "wd": dict(wd),
            "selfdual": selfdual, "typeII": selfdual and de, "had": had, "skew": skw}


def main():
    print("verify_thm481_tri_agentv.py — fresh adversarial recomputation")
    print("=" * 78)

    # ---------------- Part 1: bordered Paley_31 ----------------
    HP32 = bordered_paley(31)
    P = analyze("bordered Paley_31, tournament gauge", HP32, CLAIM_PALEY32)

    # ---------------- Part 3: the THM-447 tower ----------------
    H8 = bordered_paley(7)
    T8 = analyze("tower H_8 = border(Paley_7)", H8, CLAIM_E8)
    H16 = double_skew(H8)
    T16 = analyze("tower H_16 = double(H_8)", H16, CLAIM_TOWER16)
    H32 = double_skew(H16)
    T32 = analyze("tower H_32 = double(H_16)", H32, CLAIM_TOWER32)

    print("\n=== d32+ closed-form cross-check ===")
    print(f"  d32+ enumerator (pair-doubling C(16,j) at weight 2j, j even, + 2^15 glue at 16): "
          f"{sorted((k, v) for k, v in D32PLUS.items() if v)}")
    print(f"  tower-32 distribution == d32+ closed form: {T32['wd'] == {k: v for k, v in D32PLUS.items() if v}}")
    s1 = sum(P['wd'].values()); s2 = sum(T32['wd'].values())
    print(f"  codeword-count sanity: Paley32 {s1} == 2^16: {s1 == 65536}; tower32 {s2} == 2^16: {s2 == 65536}")

    # ---------------- code inequivalence ----------------
    print("\n=== code comparison (order 32) ===")
    print(f"  Paley_31 code  A4 = {P['wd'].get(4, 0)}, A8 = {P['wd'].get(8, 0)}")
    print(f"  tower-32 code  A4 = {T32['wd'].get(4, 0)}, A8 = {T32['wd'].get(8, 0)}")
    diff = P["wd"] != T32["wd"]
    print(f"  enumerators differ: {diff}  => codes INEQUIVALENT (enumerator is a code-equivalence invariant)")

    # ---------------- gauge caveats, empirically ----------------
    print("\n=== gauge probes (does the code move inside the Hadamard class?) ===")
    for nm, R in (("Paley32", P), ("tower32", T32)):
        rn, cc, cd = gauge_probes(R["rows"], 32, R["canon"])
        print(f"  {nm}: span fixed under complementing any single row: {rn}; "
              f"columns whose negation CHANGES the span: {len(cc)} of 32; "
              f"dims of column-negated spans: {cd}")

    # ---------------- full Hadamard-equivalence invariant ----------------
    print("\n=== 4-profile (genuine Hadamard-equivalence invariant) ===")
    pr_P = four_profile(HP32)
    pr_T = four_profile(H32)
    pc_P = four_profile(HP32, on_columns=True)
    pc_T = four_profile(H32, on_columns=True)
    print(f"  Paley32 row 4-profile:  {sorted(pr_P.items())}")
    print(f"  tower32 row 4-profile:  {sorted(pr_T.items())}")
    print(f"  Paley32 col 4-profile:  {sorted(pc_P.items())}")
    print(f"  tower32 col 4-profile:  {sorted(pc_T.items())}")
    neq = pr_P != pr_T
    neq_t = pr_P != pc_T  # also rule out Paley32 ~ tower32^T
    print(f"  row-profiles differ: {neq}  => H matrices Hadamard-INEQUIVALENT (PHQ, signed perms)")
    print(f"  Paley32 row-profile != tower32 column-profile: {neq_t}  => also not equivalent to the TRANSPOSE")

    # ---------------- verdict ----------------
    print("\n" + "=" * 78)
    okP = (P["typeII"] and P["had"] and P["skew"] and len(P["basis"]) == 16
           and P["wd"] == CLAIM_PALEY32)
    okT = (T32["typeII"] and T32["had"] and T32["skew"] and len(T32["basis"]) == 16
           and T32["wd"] == CLAIM_TOWER32 and T32["wd"] == {k: v for k, v in D32PLUS.items() if v})
    print(f"VERDICT part 1 (Paley_31 [32,16] Type II, claimed enumerator): {'CONFIRMED' if okP else 'FAILED'}")
    print(f"VERDICT part 3 (tower-32 Type II, d32+ enumerator A4=120):     {'CONFIRMED' if okT else 'FAILED'}")
    print(f"VERDICT codes different: {'CONFIRMED' if diff else 'FAILED'}; "
          f"Hadamard matrices inequivalent (4-profile): {'CONFIRMED' if (neq and neq_t) else 'NOT SEPARATED BY 4-PROFILE'}")
    return 0 if (okP and okT and diff) else 1


if __name__ == "__main__":
    sys.exit(main())
