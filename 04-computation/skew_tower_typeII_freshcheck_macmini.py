#!/usr/bin/env python3
"""
FRESH, self-contained verification (no reuse of project scripts) that the
skew-Walsh tower's tournament-gauge row code C(H_{2^k}) is a doubly-even
self-dual (Type II) [2^k, 2^{k-1}, 4] code, for k = 3,4,5,6 (orders 8,16,32,64).

Construction:
  - Base H_1 = [[1]]  (order 1, skew-type: H+H^T = 2I trivially).
  - Doubling: H_{2m} = [[ H ,  H ],
                        [-H^T, H^T]]   (skew-type preserved: H+H^T=2I).
  - Tournament-gauge binary code: C(H) = rowspan over F_2 of B = (J - H)/2,
    where J is all-ones. (entry -1 -> bit 1, entry +1 -> bit 0).

We implement our own F_2 linear algebra from scratch (Gaussian elimination over
GF(2) using Python ints as bit-vectors). No numpy linalg, no project code.

Checks per order n=2^k:
  (a) dim C = n/2
  (b) self-orthogonal: every pair of rows has even F_2 inner product
  (c) self-dual: C == C^perp  (equivalently dim = n/2 AND self-orthogonal)
  (d) doubly-even (Type II): every codeword has weight divisible by 4
  (e) minimum distance = 4 (min nonzero weight)
  (f) full weight enumerator (for orders up to 32; order 64 enumerated by
      basis-walk since 2^32 codewords is too many -> we instead verify
      doubly-even and self-dual structurally + min weight by a bounded search).
"""

# ----------------------------------------------------------------------
# Build the skew tower as plain Python lists of lists of ints (+1/-1).
# ----------------------------------------------------------------------
def transpose(M):
    n = len(M); m = len(M[0])
    return [[M[i][j] for i in range(n)] for j in range(m)]

def double(H):
    """H_{2m} = [[H, H],[-H^T, H^T]]."""
    m = len(H)
    HT = transpose(H)
    top = [H[i] + H[i] for i in range(m)]                      # [H | H]
    bot = [[-HT[i][j] for j in range(m)] + HT[i][:] for i in range(m)]  # [-H^T | H^T]
    return top + bot

def skew_tower(k):
    """Return H of order 2^k starting from H_1=[[1]]."""
    H = [[1]]
    for _ in range(k):
        H = double(H)
    return H

# ----------------------------------------------------------------------
# Sanity: skew-Hadamard property at the +-1 level (independent of F_2).
# ----------------------------------------------------------------------
def matmul(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def check_skew_hadamard(H):
    n = len(H)
    HT = transpose(H)
    # H + H^T = 2I ?
    skew_ok = all((H[i][j] + HT[i][j]) == (2 if i==j else 0) for i in range(n) for j in range(n))
    # H H^T = n I ?
    P = matmul(H, HT)
    orth_ok = all(P[i][j] == (n if i==j else 0) for i in range(n) for j in range(n))
    return skew_ok, orth_ok

# ----------------------------------------------------------------------
# Binarize: B[i] = (J - H)/2 as an integer bitmask (bit j set <=> -1 at (i,j)).
# ----------------------------------------------------------------------
def binarize_rows(H):
    n = len(H)
    rows = []
    for i in range(n):
        v = 0
        for j in range(n):
            if H[i][j] == -1:
                v |= (1 << j)
            elif H[i][j] == 1:
                pass
            else:
                raise ValueError("non +-1 entry")
        rows.append(v)
    return rows  # list of int bitmasks, length n

def popcount(x):
    return bin(x).count("1")

# ----------------------------------------------------------------------
# GF(2) row reduction on list of int bitmasks -> reduced basis + rank.
# ----------------------------------------------------------------------
def gf2_basis(vectors, nbits):
    """Return a list of basis bitmasks (row-echelon by pivot bit) spanning the
    same F_2 space. Length of returned list = dimension."""
    basis = []          # pivot_bit -> vector, kept as dict
    pivots = {}
    for v in vectors:
        x = v
        # reduce against current pivots (scan from high bit to low for determinism)
        for pb in sorted(pivots.keys(), reverse=True):
            if (x >> pb) & 1:
                x ^= pivots[pb]
        if x != 0:
            pb = x.bit_length() - 1   # highest set bit = pivot
            pivots[pb] = x
    return list(pivots.values())

def gf2_rank(vectors, nbits):
    return len(gf2_basis(vectors, nbits))

def in_span(x, basis_pivots):
    """basis_pivots: dict pivot_bit->vec. Return True if x in span."""
    for pb in sorted(basis_pivots.keys(), reverse=True):
        if (x >> pb) & 1:
            x ^= basis_pivots[pb]
    return x == 0

def basis_as_pivot_dict(basis):
    d = {}
    for v in basis:
        pb = v.bit_length() - 1
        d[pb] = v
    return d

# ----------------------------------------------------------------------
# Dual code: C^perp = { x : <x, c> = 0 mod 2 for all c in basis of C }.
# We build C^perp's basis via the parity-check (the basis of C is the
# parity-check rows). Solve null space over GF(2).
# ----------------------------------------------------------------------
def gf2_nullspace_basis(rows, nbits):
    """rows: list of int bitmasks (the generator/check matrix). Return a basis
    (list of int bitmasks) of { x in F_2^nbits : rows[i] . x = 0 for all i }."""
    # Build augmented system; do elimination tracking which columns are pivots.
    M = [r for r in rows]
    # Gaussian elimination, record pivot columns
    pivot_cols = []
    row_idx = 0
    cols_order = list(range(nbits))  # eliminate column 0,1,2,... (bit positions)
    M = M[:]  # copy
    used = [False]*len(M)
    pivot_for_col = {}
    r = 0
    R = len(M)
    for c in range(nbits):
        # find a row with bit c set, at or below r
        sel = -1
        for i in range(r, R):
            if (M[i] >> c) & 1:
                sel = i; break
        if sel == -1:
            continue
        M[r], M[sel] = M[sel], M[r]
        # eliminate bit c from all other rows
        for i in range(R):
            if i != r and ((M[i] >> c) & 1):
                M[i] ^= M[r]
        pivot_for_col[c] = r
        pivot_cols.append(c)
        r += 1
        if r == R:
            break
    free_cols = [c for c in range(nbits) if c not in pivot_for_col]
    # For each free column, build a null vector
    null_basis = []
    for fc in free_cols:
        x = (1 << fc)
        # for each pivot column, set its bit so that constraints satisfied:
        # row pivot_for_col[c] currently has bit c plus some free bits; the
        # value of x at pivot column c = sum over free bits in that row.
        for c, pr in pivot_for_col.items():
            # bit c of x = (M[pr] AND x_free) parity, excluding pivot bits
            # M[pr] has bit c set; the rest of its set bits are free columns.
            rowmask = M[pr] & ~(1 << c)
            if popcount(rowmask & x) & 1:
                x |= (1 << c)
            # else leave 0
        null_basis.append(x)
    return null_basis

# ----------------------------------------------------------------------
# Enumerate all codewords from a basis (only feasible for dim <= ~22).
# ----------------------------------------------------------------------
def enumerate_weights(basis):
    dim = len(basis)
    weights = {}
    # iterate over all 2^dim combinations via Gray-ish incremental xor
    # simple version: subset-sum
    cur = 0
    # use binary counting; for dim up to ~20 fine
    for mask in range(1 << dim):
        # incremental: compute via lowest changed bit
        x = 0
        m = mask
        idx = 0
        while m:
            if m & 1:
                x ^= basis[idx]
            m >>= 1
            idx += 1
        w = popcount(x)
        weights[w] = weights.get(w, 0) + 1
    return weights

def enumerate_weights_fast(basis):
    """Gray-code incremental enumeration of all codewords' weights."""
    dim = len(basis)
    weights = {}
    x = 0
    weights[0] = 1
    # Gray code: at step i (1..2^dim-1), flip basis[trailing_zeros(i)]
    prev = 0
    for i in range(1, 1 << dim):
        bit = (i & -i).bit_length() - 1
        x ^= basis[bit]
        w = popcount(x)
        weights[w] = weights.get(w, 0) + 1
    return weights

# ----------------------------------------------------------------------
# Min weight by bounded codeword search (for order 64 where full enum is huge).
# We search all codewords of weight < 4 by checking low-order combinations
# AND confirm there is a weight-4 codeword. Since the code is small-dim only
# for n<=32 we full-enumerate; for n=64 (dim 32) we use a randomized + greedy
# bounded check plus the structural guarantee (doubly-even => min weight in
# {4,8,...}; we exhibit a weight-4 word).
# ----------------------------------------------------------------------
def min_weight_bruteforce(basis):
    """Full enumeration min nonzero weight (only call when dim small)."""
    dim = len(basis)
    best = None
    x = 0
    for i in range(1, 1 << dim):
        bit = (i & -i).bit_length() - 1
        x ^= basis[bit]
        w = popcount(x)
        if w > 0 and (best is None or w < best):
            best = w
    return best

def min_weight_meet_in_middle(basis):
    """Min nonzero weight via meet-in-the-middle: split basis in halves.
    Feasible for dim up to ~40 (here dim=32 -> 2^16 each side)."""
    dim = len(basis)
    h = dim // 2
    A = basis[:h]
    Bb = basis[h:]
    # Enumerate all combos of A -> dict mapping vector -> weight (keep min per vec)
    from itertools import product
    # Build all A-combos
    a_vecs = []
    x = 0
    a_vecs.append(0)
    for i in range(1, 1 << len(A)):
        bit = (i & -i).bit_length() - 1
        x ^= A[bit]
        a_vecs.append(x)
    # store min over A side as list aligned to gray order; we want, for the
    # full space, min popcount(a xor b). Brute O(2^h * 2^(dim-h)) is 2^32 -> too big.
    # Instead: min weight of a self-dual doubly-even code is >=4; we just need to
    # confirm min weight == 4 by EXHIBITING a weight-4 codeword and ruling out
    # weight<4. Doubly-even already rules out weights 1,2,3 (not divisible by 4,
    # except weight 0). So if doubly-even is verified, min nonzero weight is a
    # multiple of 4, hence >=4. Then exhibit one weight-4 word.
    raise NotImplementedError

def exhibit_weight4(basis):
    """Find a codeword of weight exactly 4 by searching pairs/singles/triples
    of basis vectors (bounded). Returns the word or None."""
    dim = len(basis)
    # singles
    for i in range(dim):
        if popcount(basis[i]) == 4:
            return basis[i]
    # pairs
    for i in range(dim):
        for j in range(i+1, dim):
            x = basis[i] ^ basis[j]
            if popcount(x) == 4:
                return x
    # triples
    for i in range(dim):
        for j in range(i+1, dim):
            xij = basis[i] ^ basis[j]
            for l in range(j+1, dim):
                x = xij ^ basis[l]
                if popcount(x) == 4:
                    return x
    return None

# ----------------------------------------------------------------------
# Main verification loop.
# ----------------------------------------------------------------------
def verify_order(k, full_enum_max_dim=20):
    n = 1 << k
    H = skew_tower(k)
    assert len(H) == n
    skew_ok, orth_ok = check_skew_hadamard(H)

    rows = binarize_rows(H)          # n bitmasks
    basis = gf2_basis(rows, n)
    dim = len(basis)

    # self-orthogonal: check all pairs of GENERATING rows (suffices, since
    # inner product is bilinear; but to be safe check basis pairs AND row pairs)
    def self_orthogonal(vecs):
        ok = True
        L = len(vecs)
        for i in range(L):
            if popcount(vecs[i] & vecs[i]) % 2 != 0:   # self inner product = weight mod 2
                ok = False
            for j in range(i+1, L):
                if popcount(vecs[i] & vecs[j]) % 2 != 0:
                    ok = False
        return ok
    so_rows = self_orthogonal(rows)
    so_basis = self_orthogonal(basis)

    # self-dual: C == C^perp.  C^perp basis = nullspace of generator rows.
    perp_basis = gf2_nullspace_basis(basis, n)
    perp_dim = len(perp_basis)
    # check C subset C^perp (self-orthogonal) AND dims equal
    bp = basis_as_pivot_dict(basis)
    # every perp basis vector in C ?  and every C basis vector in C^perp ?
    # self-dual iff C == C^perp. Test: dim C == dim C^perp == n/2 and C ⊆ C^perp.
    C_in_perp = so_basis  # C ⊆ C^perp iff self-orthogonal
    perp_in_C = all(in_span(v, bp) for v in perp_basis)
    self_dual = (dim == perp_dim == n//2) and C_in_perp and perp_in_C

    # doubly-even + weight enumerator
    result = {
        "order": n, "k": k, "skew_type_ok": skew_ok, "hadamard_orth_ok": orth_ok,
        "dim": dim, "n_over_2": n//2, "dim_ok": dim == n//2,
        "self_orth_rows": so_rows, "self_orth_basis": so_basis,
        "perp_dim": perp_dim, "self_dual": self_dual,
    }

    if dim <= full_enum_max_dim:
        weights = enumerate_weights_fast(basis)
        doubly_even = all((w % 4 == 0) for w in weights)
        nz = [w for w in weights if w > 0]
        min_d = min(nz) if nz else None
        result["weights"] = dict(sorted(weights.items()))
        result["doubly_even"] = doubly_even
        result["min_distance"] = min_d
        result["num_codewords"] = sum(weights.values())
    else:
        # order 64: dim=32, cannot enumerate 2^32. Verify doubly-even on basis
        # + all pairwise XORs is NOT sufficient for full doubly-even; instead
        # use the standard fact: a binary code with (i) every basis vector
        # weight divisible by 4 and (ii) self-orthogonal (every pair of basis
        # vectors has inner product 0 mod 2) is doubly-even. (wt(a+b)=wt(a)+wt(b)
        # -2<a,b>; mod 4: if wt(a),wt(b) ≡0 mod4 and <a,b> even then wt(a+b)≡0 mod4.)
        basis_all_div4 = all(popcount(v) % 4 == 0 for v in basis)
        doubly_even = basis_all_div4 and so_basis
        result["basis_weights_mod4"] = [popcount(v) % 4 for v in basis]
        result["basis_all_div4"] = basis_all_div4
        result["doubly_even_structural"] = doubly_even
        # min distance: doubly-even => min nonzero weight is multiple of 4 => >=4.
        # exhibit a weight-4 codeword:
        w4 = exhibit_weight4(basis)
        result["weight4_witness_found"] = (w4 is not None)
        result["min_distance"] = 4 if (w4 is not None and doubly_even) else None
    return result

if __name__ == "__main__":
    import json
    for k in (3, 4, 5, 6):
        r = verify_order(k, full_enum_max_dim=20)  # dim 4,8,16 enumerated; dim32 structural
        print(f"=== order {r['order']} (k={k}) ===")
        for key, val in r.items():
            if key in ("order", "k"): continue
            print(f"  {key}: {val}")
        print()
