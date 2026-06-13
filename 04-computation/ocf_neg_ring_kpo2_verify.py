# ocf_neg_ring_kpo2_verify.py
# ADVERSARIAL VERIFIER, thread C (HYP-2380), claims C1 (det quotient) and C4 (path GF
# reciprocity). FRESH code, independent of the worker's.
#
# Ring R = Z[x_1..x_n]/(x_i^2): elements are dicts {vertex-bitmask: int coeff}.
# det via Leibniz directly (each permutation contributes a single +-monomial since
# entry (i,j) of I+XA is delta_ij + x_i A_ij). Inverse via Neumann series with an
# explicit final assertion u * u^{-1} == 1. No division anywhere.
#
# Claims tested as FULL RING EQUALITIES (stronger than the eps-functional claim):
#   (i)   det(I+XA) * inv(det(I-XA)) == G_2   (collection GF at t=2)
#   (ii)  det(I-XA) * inv(det(I+XA)) == G_{-2}
#   (iii) [quotient]^2 == G_4 and [reciprocal]^2 == G_{-4}   (general even-t claim)
#   (iv)  eps of (i)/(ii) == I(Omega,2)=H and I(Omega,-2) from an independent census
#   (v)   P(x) * P(-x) == 1 where P = simple-directed-path GF with vertex marks
#         (equivalent to W_T(z) W_T(-z) = 1 in R[z], since z^k always pairs with a
#         degree-k squarefree monomial)
#   (vi)  edge cases n=1,2
# where G_t = sum over collections K of pairwise vertex-disjoint DIRECTED odd cycles
# of t^{|K|} x^{union of supports} (MISTAKE-023 guard: directed cycles, not supports).
#
# Scope (runtime guard): ALL iso-class representatives n=3..5, all 56 reps n=6,
# plus 10 random labeled tournaments per n in 3..7.

import itertools, random
from collections import Counter

random.seed(98765)

def pairs_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_mask(n, mask, pairs):
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    for b, (i, j) in enumerate(pairs):
        if (mask >> b) & 1:
            adj[i][j] = 1; out[i] |= 1 << j
        else:
            adj[j][i] = 1; out[j] |= 1 << i
    return adj, out

# ---------- ring ops ----------
def rmul(a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            if ma & mb == 0:
                k = ma | mb
                out[k] = out.get(k, 0) + ca * cb
    return {k: v for k, v in out.items() if v != 0}

def radd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}

def rneg(a):
    return {k: -v for k, v in a.items()}

ONE = {0: 1}

def rinv(u):
    """inverse of u with constant term 1, by Neumann series; asserts exactness."""
    assert u.get(0, 0) == 1
    v = {k: c for k, c in u.items() if k != 0}
    inv = dict(ONE)
    term = dict(ONE)
    while True:
        term = rneg(rmul(term, v))
        if not term:
            break
        inv = radd(inv, term)
    assert rmul(u, inv) == ONE, "Neumann inverse failed exactness"
    return inv

def eps(a):
    return sum(a.values())

def det_IpXA(n, adj, sign_arc):
    """det(I + sign_arc * X A) by Leibniz; returns ring element."""
    out = {}
    for p in itertools.permutations(range(n)):
        # sign of permutation
        seen = [False] * n
        sgn = 1
        mono = 0
        ok = True
        arcs = 0
        for i in range(n):
            if p[i] != i:
                if not adj[i][p[i]]:
                    ok = False
                    break
                mono |= 1 << i
                arcs += 1
        if not ok:
            continue
        # permutation sign
        for i in range(n):
            if not seen[i]:
                clen = 0
                j = i
                while not seen[j]:
                    seen[j] = True
                    j = p[j]
                    clen += 1
                if clen % 2 == 0:
                    sgn = -sgn
        coeff = sgn * (sign_arc ** arcs)
        out[mono] = out.get(mono, 0) + coeff
    return {k: v for k, v in out.items() if v != 0}

# ---------- combinatorial sides ----------
def odd_cycle_supports(n, adj):
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
                    res.append(m)
    return res

def collection_GF(cycles, t):
    """sum over collections of pairwise disjoint directed odd cycles of t^{|K|} x^{supp}."""
    out = {0: 1}
    L = sorted(range(len(cycles)))
    def rec(idx, mask, size):
        for i in range(idx, len(cycles)):
            if cycles[i] & mask == 0:
                m2 = mask | cycles[i]
                out[m2] = out.get(m2, 0) + t ** (size + 1)
                rec(i + 1, m2, size + 1)
    rec(0, 0, 0)
    return {k: v for k, v in out.items() if v != 0}

def path_GF(n, adj, negate):
    """sum over simple directed paths (incl. empty) of prod x_v, with optional
    per-vertex sign (-1) (i.e. P(-x)). Brute force over injective sequences."""
    out = {0: 1}
    for k in range(1, n + 1):
        for seq in itertools.permutations(range(n), k):
            if all(adj[seq[i]][seq[i + 1]] for i in range(k - 1)):
                m = 0
                for v in seq:
                    m |= 1 << v
                c = (-1) ** k if negate else 1
                out[m] = out.get(m, 0) + c
    return {k: v for k, v in out.items() if v != 0}

def ham_dfs(n, out_):
    full = (1 << n) - 1
    total = 0
    stack = [(v, 1 << v) for v in range(n)]
    while stack:
        v, vis = stack.pop()
        if vis == full:
            total += 1
            continue
        avail = out_[v] & ~vis
        while avail:
            b = avail & -avail
            stack.append((b.bit_length() - 1, vis | b))
            avail ^= b
    return total

def canon_reps(n, pairs):
    """iso-class representatives by brute canonical minimization (n<=6)."""
    seen = set()
    reps = []
    perms = list(itertools.permutations(range(n)))
    for mask in range(1 << len(pairs)):
        adj, _ = adj_from_mask(n, mask, pairs)
        best = None
        scores = [sum(adj[v]) for v in range(n)]
        groups = {}
        for v in range(n):
            groups.setdefault(scores[v], []).append(v)
        keys = sorted(groups)
        for parts in itertools.product(*[itertools.permutations(groups[k]) for k in keys]):
            perm = [v for part in parts for v in part]
            m = 0; b = 0
            for i in range(n):
                pi = perm[i]
                for j in range(i + 1, n):
                    if adj[pi][perm[j]]:
                        m |= 1 << b
                    b += 1
            if best is None or m < best:
                best = m
        if best not in seen:
            seen.add(best)
            reps.append(mask)
    return reps

print("=" * 78)
print("RING-IDENTITY VERIFICATION (C1 det quotient, C4 path reciprocity) -- fresh code")
print("=" * 78)

# edge cases n=1,2
for n in (1, 2):
    pairs = pairs_list(n)
    for mask in range(1 << len(pairs)):
        adj, out_ = adj_from_mask(n, mask, pairs)
        dp = det_IpXA(n, adj, +1); dm = det_IpXA(n, adj, -1)
        q = rmul(dp, rinv(dm)); qr = rmul(dm, rinv(dp))
        assert q == ONE and qr == ONE
        assert eps(q) == 1 == ham_dfs(n, out_)
        P = path_GF(n, adj, False); Pm = path_GF(n, adj, True)
        assert rmul(P, Pm) == ONE
print("n=1,2 edge cases: det quotient == 1 == I(Omega,t) (Omega empty), H=1, "
      "path reciprocity holds")

total = Counter()
for n in range(3, 8):
    pairs = pairs_list(n)
    if n <= 6:
        masks = canon_reps(n, pairs)
    else:
        masks = []
    masks = list(masks) + [random.getrandbits(len(pairs)) for _ in range(10)]
    nfail = 0
    for mask in masks:
        adj, out_ = adj_from_mask(n, mask, pairs)
        cycles = odd_cycle_supports(n, adj)
        dp = det_IpXA(n, adj, +1)
        dm = det_IpXA(n, adj, -1)
        inv_dm = rinv(dm)
        inv_dp = rinv(dp)
        Q = rmul(dp, inv_dm)       # claimed = G_2
        QR = rmul(dm, inv_dp)      # claimed = G_{-2}
        G2 = collection_GF(cycles, 2)
        Gm2 = collection_GF(cycles, -2)
        G4 = collection_GF(cycles, 4)
        Gm4 = collection_GF(cycles, -4)
        ok_i = (Q == G2)
        ok_ii = (QR == Gm2)
        ok_iii = (rmul(Q, Q) == G4) and (rmul(QR, QR) == Gm4)
        H = ham_dfs(n, out_)
        a1 = len(cycles)
        a2 = sum(1 for i in range(a1) for j in range(i + 1, a1)
                 if cycles[i] & cycles[j] == 0)
        a3 = 0
        for i in range(a1):
            for j in range(i + 1, a1):
                if cycles[i] & cycles[j] == 0:
                    u = cycles[i] | cycles[j]
                    for k in range(j + 1, a1):
                        if cycles[k] & u == 0:
                            a3 += 1
        Im2 = 1 - 2 * a1 + 4 * a2 - 8 * a3
        ok_iv = (eps(Q) == H == 1 + 2 * a1 + 4 * a2 + 8 * a3) and (eps(QR) == Im2)
        P = path_GF(n, adj, False)
        Pm = path_GF(n, adj, True)
        ok_v = (rmul(P, Pm) == ONE)
        for name, ok in (('i:Q==G2', ok_i), ('ii:QR==Gm2', ok_ii),
                         ('iii:t=+-4', ok_iii), ('iv:eps==H,I(-2)', ok_iv),
                         ('v:path-reciprocity', ok_v)):
            total[name + ('_PASS' if ok else '_FAIL')] += 1
            if not ok:
                nfail += 1
                print(f"  FAIL n={n} mask={mask} check {name}")
    print(f"n={n}: {len(masks)} tournaments "
          f"({'all iso reps + ' if n<=6 else ''}10 random), failures: {nfail}")

print("\nsummary:", dict(sorted(total.items())))
nf = sum(v for k, v in total.items() if k.endswith('_FAIL'))
print(f"TOTAL FAILURES: {nf}")
print("RING CHECKS:", "ALL PASS" if nf == 0 else "FAILURES PRESENT")
print("done")
