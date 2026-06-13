#!/usr/bin/env python3
"""
ocf_digit_tower_small_kpo2_verify.py  -- ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S2)

Independent verification of THM-466 claims (thread B-2adic-digit-tower) for n <= 6
(full labeled censuses) plus the explicit n=5 score-class counterexample pair.

FRESH METHODS (deliberately different from the worker's Held-Karp DP script):
  * H computed by INCLUSION-EXCLUSION over vertex subsets (walk counting):
        H = sum_{S subseteq V} (-1)^{n-|S|} * (total # walks of length n-1 inside S)
    vectorized in numpy int64 (exact; entries far below 2^63).
  * H cross-checked by PURE brute force over all n! permutations (pure Python)
    on random codes and on the counterexample pair.
  * Directed odd cycles enumerated as DIRECTED cycles: for each odd-size vertex
    subset, every cyclic order with the minimum element fixed first is tested
    (so a vertex set CAN host several directed cycles -- repo mistake guard).
  * alpha_2 = number of unordered pairs of vertex-disjoint directed odd cycles,
    summed over disjoint subset pairs of per-subset directed-Hamiltonian-cycle counts.

Encoding (MINE, chosen independently): for u < v, bit p(u,v) (lexicographic pair
order) is 1 iff u -> v.  T^op = bitwise complement of the code, so full-census
reversal invariance <=> the alpha arrays are PALINDROMES in code order.

Assertions (all full labeled censuses):
  n <= 5 : H == 1 + 2*alpha1            (alpha2 == 0)
  n == 6 : H == 1 + 2*alpha1 + 4*alpha2 (only (3,3) pairs exist)
  Sum H == n! * 2^(m-n+1); Sum c3 == C(n,3)*2^(m-2); Sum c5 == C(n,5)*24*2^(m-5)
  #(H==1) == n!  (transitive); H always odd (Redei); Kendall c3 == C(n,3)-sum C(s_i,2)
  v2(H-1) distributions; score-class constancy of alpha1/c3/c5 parities;
  reversal palindromes; the explicit counterexample pair.
"""
import sys, math, itertools
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
FAIL = 0

def chk(cond, msg):
    global FAIL
    if cond:
        print("  PASS  " + msg)
    else:
        FAIL += 1
        print("  *** FAIL *** " + msg)

# ---------- bit encoding ----------
def pair_index(n):
    idx = {}
    p = 0
    for u in range(n):
        for v in range(u + 1, n):
            idx[(u, v)] = p
            p += 1
    return idx, p

def adj_from_code(code, n, idx):
    A = [[0] * n for _ in range(n)]
    for (u, v), p in idx.items():
        if (code >> p) & 1:
            A[u][v] = 1
        else:
            A[v][u] = 1
    return A

# ---------- pure-Python reference engines (brute force) ----------
def ham_paths_perm(A, n):
    cnt = 0
    for perm in itertools.permutations(range(n)):
        ok = True
        for i in range(n - 1):
            if not A[perm[i]][perm[i + 1]]:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt

def dir_ham_cycles_subset(A, sub):
    """# directed Hamiltonian cycles of induced subtournament on tuple sub."""
    first = sub[0]
    cnt = 0
    for perm in itertools.permutations(sub[1:]):
        seq = (first,) + perm
        ok = all(A[seq[i]][seq[(i + 1) % len(seq)]] for i in range(len(seq)))
        if ok:
            cnt += 1
    return cnt

def alphas_brute(A, n):
    """List all directed odd cycles, then count alpha1, alpha2 by explicit
    disjointness of vertex sets (cycles carried as (frozenset, id))."""
    cycles = []
    for s in (3, 5, 7):
        if s > n:
            break
        for sub in itertools.combinations(range(n), s):
            first = sub[0]
            for perm in itertools.permutations(sub[1:]):
                seq = (first,) + perm
                if all(A[seq[i]][seq[(i + 1) % s]] for i in range(s)):
                    cycles.append(frozenset(sub))
    a1 = len(cycles)
    a2 = 0
    for i in range(a1):
        for j in range(i + 1, a1):
            if not (cycles[i] & cycles[j]):
                a2 += 1
    return a1, a2

# ---------- vectorized census engines ----------
def build_B(codes, n, idx):
    B = {}
    for (u, v), p in idx.items():
        b = ((codes >> p) & 1).astype(np.uint8)
        B[(u, v)] = b
        B[(v, u)] = (1 - b).astype(np.uint8)
    return B

def census_c3(B, n):
    """per-triple cyclic indicators (dict) and total c3 array."""
    N = len(next(iter(B.values())))
    cyc3 = {}
    c3 = np.zeros(N, dtype=np.int32)
    for t in itertools.combinations(range(n), 3):
        i, j, k = t
        f = B[(i, j)] & B[(j, k)] & B[(k, i)]
        g = B[(i, k)] & B[(k, j)] & B[(j, i)]
        cyc3[t] = f | g
        c3 += cyc3[t]
    return cyc3, c3

def census_cycles_size(B, n, s):
    """per-subset directed-Hamiltonian-cycle COUNT arrays for subsets of size s."""
    N = len(next(iter(B.values())))
    out = {}
    for sub in itertools.combinations(range(n), s):
        acc = np.zeros(N, dtype=np.uint16)
        first = sub[0]
        for perm in itertools.permutations(sub[1:]):
            seq = (first,) + perm
            ind = B[(seq[0], seq[1])].copy()
            for i in range(1, s):
                ind &= B[(seq[i], seq[(i + 1) % s])]
            acc += ind
        out[sub] = acc
    return out

def batched_power(M, e):
    """M: (N,s,s) int64, return M^e (e>=1) by binary exponentiation."""
    result = None
    base = M
    while e:
        if e & 1:
            result = base if result is None else np.matmul(result, base)
        e >>= 1
        if e:
            base = np.matmul(base, base)
    return result

def H_inclusion_exclusion(codes, n, idx):
    N = len(codes)
    A = np.zeros((N, n, n), dtype=np.int64)
    for (u, v), p in idx.items():
        b = (codes >> p) & 1
        A[:, u, v] = b
        A[:, v, u] = 1 - b
    H = np.zeros(N, dtype=np.int64)
    for r in range(1, n + 1):
        sign = (-1) ** (n - r)
        for S in itertools.combinations(range(n), r):
            Ssub = A[:, S, :][:, :, S]
            if n == 1:
                W = np.ones(N, dtype=np.int64)
            else:
                P = batched_power(Ssub, n - 1)
                W = P.sum(axis=(1, 2))
            H += sign * W
    return H

def v2_of(x):
    """2-adic valuation array for positive ints (int64)."""
    lo = (x & -x).astype(np.float64)
    return np.round(np.log2(lo)).astype(np.int64)

def score_class_analysis(scores_sorted_key, parity_vectors, labels):
    uniq, inv = np.unique(scores_sorted_key, return_inverse=True)
    cnt = np.bincount(inv).astype(np.float64)
    res = {}
    for lab, par in zip(labels, parity_vectors):
        odd = np.bincount(inv, weights=par.astype(np.float64))
        noncon = int(((odd > 0.5) & (odd < cnt - 0.5)).sum())
        res[lab] = noncon
    return len(uniq), res

# =====================================================================
print("=" * 78)
print("PART 0: edge cases n = 1, 2, 3 (pure Python, definitions only)")
print("=" * 78)
for n in (1, 2, 3):
    idx, m = pair_index(n)
    bad = 0
    for code in range(2 ** m):
        A = adj_from_code(code, n, idx)
        H = ham_paths_perm(A, n) if n > 1 else 1
        a1, a2 = alphas_brute(A, n)
        if H != 1 + 2 * a1 + 4 * a2:
            bad += 1
    chk(bad == 0, f"n={n}: H == 1+2*a1+4*a2 on all {2**m} tournaments (a2==0 trivially)")

results_small = {}
for n in (4, 5, 6):
    print("=" * 78)
    print(f"PART n={n}: FULL labeled census, 2^{n*(n-1)//2} tournaments")
    print("=" * 78)
    idx, m = pair_index(n)
    N = 2 ** m
    codes = np.arange(N, dtype=np.int64)
    B = build_B(codes, n, idx)

    cyc3, c3 = census_c3(B, n)
    c5 = np.zeros(N, dtype=np.int32)
    cyc5 = {}
    if n >= 5:
        cyc5 = census_cycles_size(B, n, 5)
        for sub, arr in cyc5.items():
            c5 += arr
    alpha1 = c3 + c5
    alpha2 = np.zeros(N, dtype=np.int32)
    if n >= 6:
        trips = list(cyc3.keys())
        for i in range(len(trips)):
            for j in range(i + 1, len(trips)):
                if not (set(trips[i]) & set(trips[j])):
                    alpha2 += (cyc3[trips[i]] & cyc3[trips[j]]).astype(np.int32)

    H = H_inclusion_exclusion(codes, n, idx)

    # central identity
    rhs = 1 + 2 * alpha1.astype(np.int64) + 4 * alpha2.astype(np.int64)
    nf = int((H != rhs).sum())
    chk(nf == 0, f"H == 1 + 2*alpha1 + 4*alpha2 EXACTLY on all {N} tournaments (failures={nf})")
    if n <= 5:
        chk(int(alpha2.max()) == 0, "alpha2 == 0 everywhere (n<=5, so H == 1+2*alpha1 exactly)")

    # Redei / aggregates
    chk(bool((H % 2 == 1).all()) and int(H.min()) >= 1, "H always odd and >= 1 (Redei; H=0 impossible)")
    chk(int(H.sum()) == math.factorial(n) * 2 ** (m - n + 1),
        f"Sum H == n!*2^(m-n+1) == {math.factorial(n)*2**(m-n+1)}")
    chk(int((H == 1).sum()) == math.factorial(n), f"#(H==1) == n! == {math.factorial(n)} (transitive)")
    chk(int(c3.sum()) == math.comb(n, 3) * 2 ** (m - 2),
        f"Sum c3 == C(n,3)*2^(m-2) == {math.comb(n,3)*2**(m-2)}")
    if n >= 5:
        chk(int(c5.sum()) == math.comb(n, 5) * 24 * 2 ** (m - 5),
            f"Sum c5 == C(n,5)*24*2^(m-5) == {math.comb(n,5)*24*2**(m-5)}")

    # Kendall: c3 score-determined exactly
    scores = np.zeros((N, n), dtype=np.int64)
    for u in range(n):
        for v in range(n):
            if u != v:
                scores[:, u] += B[(u, v)]
    kendall = math.comb(n, 3) - (scores * (scores - 1) // 2).sum(axis=1)
    chk(bool((kendall == c3).all()), "c3 == C(n,3) - sum C(s_i,2) (Kendall) everywhere")

    # reversal invariance: code of T^op is bitwise complement -> palindrome test
    chk(bool(np.array_equal(alpha1, alpha1[::-1])) and bool(np.array_equal(alpha2, alpha2[::-1]))
        and bool(np.array_equal(H, H[::-1])),
        "FULL-census reversal invariance: alpha1, alpha2, H are palindromes under T -> T^op")

    # v2(H-1) distribution (from actual H, independent of the identity)
    Hm1 = H - 1
    inf_cnt = int((Hm1 == 0).sum())
    pos = Hm1[Hm1 > 0]
    vals = v2_of(pos)
    dist = {int(v): int(c) for v, c in zip(*np.unique(vals, return_counts=True))}
    print(f"  v2(H-1) distribution n={n}: " +
          ", ".join(f"v2={k}: {v}" for k, v in sorted(dist.items())) + f", inf: {inf_cnt}")
    p_odd = int((alpha1 % 2 == 1).sum())
    print(f"  P(alpha1 odd) = {p_odd}/{N} = {p_odd/N:.6f}")
    chk(inf_cnt == math.factorial(n), "v2(H-1)=inf exactly on the n! transitive tournaments")
    chk(dist.get(1, 0) == p_odd, "#(v2(H-1)==1) == #(alpha1 odd)  [digit-lemma corollary]")
    c5odd = int((c5 % 2 == 1).sum())
    print(f"  #(c5 odd) = {c5odd}/{N}")

    # score-class analysis
    skey = np.sort(scores, axis=1)
    key = np.zeros(N, dtype=np.int64)
    for i in range(n):
        key = key * (n + 1) + skey[:, i]
    ncls, noncon = score_class_analysis(
        key, [(alpha1 % 2).astype(np.int64), (c3 % 2).astype(np.int64), (c5 % 2).astype(np.int64),
              (H % 4).astype(np.int64) // 2],
        ["alpha1_parity", "c3_parity", "c5_parity", "Hmod4_bit"])
    print(f"  score classes: {ncls}; non-constant alpha1-parity classes: {noncon['alpha1_parity']}; "
          f"c5-parity: {noncon['c5_parity']}; c3-parity: {noncon['c3_parity']}")
    chk(noncon["c3_parity"] == 0, "c3 parity constant on every score class (it is score-determined)")
    chk(noncon["alpha1_parity"] == noncon["Hmod4_bit"],
        "alpha1 parity and H mod 4 agree class-by-class (digit lemma)")
    if n == 4:
        chk(noncon["alpha1_parity"] == 0, "n=4: alpha1 parity (hence H) IS score-determined")
    results_small[n] = dict(ncls=ncls, noncon=noncon, p_odd=p_odd, dist=dist, c5odd=c5odd)

    # spot cross-validation against pure brute force
    rng = np.random.default_rng(271828 + n)
    sample = rng.integers(0, N, size=12)
    ok = True
    for code in sample:
        A = adj_from_code(int(code), n, idx)
        hp = ham_paths_perm(A, n)
        a1b, a2b = alphas_brute(A, n)
        if hp != int(H[code]) or a1b != int(alpha1[code]) or a2b != int(alpha2[code]):
            ok = False
    chk(ok, "12 random codes: numpy census == pure-Python permutation brute force (H, alpha1, alpha2)")

print("=" * 78)
print("PART CE: the claimed n=5 counterexample pair (explicit arcs, brute force only)")
print("=" * 78)
n = 5
arcs1 = [(0, 3), (0, 4), (1, 0), (2, 0), (2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)]
A1 = [[0] * n for _ in range(n)]
for u, v in arcs1:
    A1[u][v] = 1
# completeness check
comp = all((A1[u][v] + A1[v][u] == 1) for u in range(n) for v in range(u + 1, n))
chk(comp, "T1 arc list is a tournament (every pair oriented exactly once)")
# T2 = T1 with the 3-cycle on {0,1,3} reversed
c013 = A1[0][3] and A1[3][1] and A1[1][0]
chk(bool(c013), "T1 has the directed 3-cycle 0->3->1->0 on {0,1,3}")
A2 = [row[:] for row in A1]
for (u, v) in [(0, 3), (3, 1), (1, 0)]:
    A2[u][v], A2[v][u] = 0, 1
s1 = sorted(sum(A1[u]) for u in range(n))
s2 = sorted(sum(A2[u]) for u in range(n))
chk(s1 == [1, 2, 2, 2, 3] and s2 == [1, 2, 2, 2, 3],
    f"both have sorted score sequence (1,2,2,2,3); got {s1} and {s2}")
def census_one(A):
    c3 = sum(dir_ham_cycles_subset(A, t) for t in itertools.combinations(range(5), 3))
    c5 = sum(dir_ham_cycles_subset(A, t) for t in itertools.combinations(range(5), 5))
    return c3, c5, ham_paths_perm(A, 5)
c3a, c5a, Ha = census_one(A1)
c3b, c5b, Hb = census_one(A2)
print(f"  T1: c3={c3a}, c5={c5a}, H={Ha} (H mod 4 = {Ha%4});  T2: c3={c3b}, c5={c5b}, H={Hb} (H mod 4 = {Hb%4})")
chk((c3a, c5a, Ha) == (4, 1, 11), "T1 matches claim: c3=4, c5=1, H=11 == 3 (mod 4)")
chk((c3b, c5b, Hb) == (4, 2, 13), "T2 matches claim: c3=4, c5=2, H=13 == 1 (mod 4)")
chk(Ha == 1 + 2 * (c3a + c5a) and Hb == 1 + 2 * (c3b + c5b), "both satisfy H == 1+2*alpha1 (n=5)")
chk(Ha % 4 != Hb % 4, "H mod 4 differs on the same score sequence => alpha1 parity NOT score-determined")

print("=" * 78)
print(f"SUMMARY (small): {'ALL CHECKS PASSED' if FAIL == 0 else str(FAIL) + ' FAILURES'}")
print("=" * 78)
sys.exit(0 if FAIL == 0 else 1)
