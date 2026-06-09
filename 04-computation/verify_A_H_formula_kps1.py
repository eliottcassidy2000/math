#!/usr/bin/env python3
"""
verify_A_H_formula_kps1.py — kind-pasteur-2026-06-09-S1
ADVERSARIAL VERIFICATION of branch A-H-formula (sibling agent claims).

Independent recomputation sharing NO code with skew_doubling_core_kps1 /
ip_doubling_hunt_kps1 / k2_spectrum_functional_kps1 / c7_lift_law_verify_kps1:
  - own tournament construction from upper-triangle bits
  - own Hamiltonian-path counters: (a) recursive backtracking, (b) fresh bitmask DP,
    (c) raw permutation scan (10-vertex gold standard for the headline CE numbers)
  - own doubling constructions with INTERLEAVED vertex order (i, i' -> 2i, 2i+1),
    different from the sibling's block order [copy1 | copy2] (also tests label invariance)
  - own cycle counting: permutation method (bases + cross-check) and fresh subset-DP
  - own strong components (BFS reachability), pair features, conflict-graph IP

Targets:
  T0  anchors: H(C3)=3, H(regT5)=15, H(D(C3))=45, H(D(trans3))=13, H(SC(trans3))=41
  T1  CE pair n=5 idx4 (bits 0001000000) vs idx6 (bits 0101000000):
      same IP & typed IP, H(K2) = 3225 vs 2785                       [claim 1]
  T2  74-class sweep: spectrum (c3,c4,c5,c6) -> H(K2): 0/32 broken    [claim 2]
  T3  product law H(T[K2]) = prod_C H(C[K2]) 74/74; fails for D/SC except strong [claim 3]
  T4  lift laws c3'=8c3, c5'=32c5+32c4+6c3 (74/74);
      c7' = 128c7+192c6+80c5+8c4+64p331+48p332+64p341+32p342 (74/74);
      p343 == 2*p332 (74/74)                                          [claims 4,5]
  T5  H(K2) == 2H(T)-1 mod 8 74/74, mod-16 census {0:39, 8:35};
      H(SC) == 1 mod 4 74/74                                          [claim 6]
  T6  D op-asymmetry 50/74 with CE 189 vs 333; n=5: H(D)==1 mod 8 12/12,
      n=3: ==5 mod 8                                                  [claim 7]
  T7  transitive towers H(D)=13,95,1033,15611,313285;
      H(SC)=41,629,14937,513669,24104937                              [claim 11]
  T8  spot recheck of profile/IP equality for spec-(6,7,6,3) strong pair idx42/idx45

Output: 05-knowledge/results/verify_A_H_formula_kps1.out
"""
import itertools, sys, time

sys.setrecursionlimit(10000)
OUT = open('05-knowledge/results/verify_A_H_formula_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

# ---------------- construction ----------------

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def from_bits(n, bits):
    """bits: int (bit k = pair k) or '01' string (char k = pair k); 1 -> i beats j."""
    P = pairs_of(n)
    A = [[0] * n for _ in range(n)]
    if isinstance(bits, str):
        vals = [int(ch) for ch in bits]
    else:
        vals = [(bits >> k) & 1 for k in range(len(P))]
    for k, (i, j) in enumerate(P):
        if vals[k]:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def bits_str(A):
    n = len(A)
    return ''.join(str(A[i][j]) for i, j in pairs_of(n))

def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]

def scores_mine(A):
    return tuple(sum(r) for r in A)

def iso_classes_mine(n):
    """Replicates sibling enumeration order (bits ascending, orbit sweep) so idx aligns."""
    P = pairs_of(n)
    perms = list(itertools.permutations(range(n)))
    seen = set()
    reps = []
    for b in range(1 << len(P)):
        A = from_bits(n, b)
        key = tuple(tuple(r) for r in A)
        if key in seen:
            continue
        reps.append(A)
        for p in perms:
            seen.add(tuple(tuple(A[p[i]][p[j]] for j in range(n)) for i in range(n)))
    return reps

# ---------------- Hamiltonian path counters (3 independent methods) ----------------

def H_backtrack(A):
    n = len(A)
    adj = [[j for j in range(n) if A[i][j]] for i in range(n)]
    count = 0
    def rec(v, visited, depth):
        nonlocal count
        if depth == n:
            count += 1
            return
        for u in adj[v]:
            b = 1 << u
            if not (visited & b):
                rec(u, visited | b, depth + 1)
    for s in range(n):
        rec(s, 1 << s, 1)
    return count

def H_dp(A):
    n = len(A)
    om = []
    for i in range(n):
        m = 0
        for j in range(n):
            if A[i][j]:
                m |= 1 << j
        om.append(m)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, size):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c:
                free = om[last] & ~mask
                while free:
                    lb = free & (-free)
                    j = lb.bit_length() - 1
                    dp[mask | lb][j] += c
                    free ^= lb
    return sum(dp[size - 1])

def H_perm(A):
    """Raw gold standard: scan all n! permutations."""
    n = len(A)
    cnt = 0
    for p in itertools.permutations(range(n)):
        for k in range(n - 1):
            if not A[p[k]][p[k + 1]]:
                break
        else:
            cnt += 1
    return cnt

# ---------------- doublings, interleaved vertex order ----------------

def K2_mine(A):
    """Lex product T[K2]: twin arc i->i'; all arcs between twin-classes follow T."""
    n = len(A)
    N = 2 * n
    B = [[0] * N for _ in range(N)]
    for i in range(n):
        B[2 * i][2 * i + 1] = 1
        for j in range(n):
            if i != j and A[i][j]:
                for a in (0, 1):
                    for b in (0, 1):
                        B[2 * i + a][2 * j + b] = 1
    return B

def D_mine(A):
    """Skew double: copy1 = T, copy2 = T^op, BOTH cross directions follow T, twin i->i'."""
    n = len(A)
    N = 2 * n
    B = [[0] * N for _ in range(N)]
    for i in range(n):
        B[2 * i][2 * i + 1] = 1
        for j in range(n):
            if i != j and A[i][j]:
                B[2 * i][2 * j] = 1
                B[2 * j + 1][2 * i + 1] = 1
                B[2 * i][2 * j + 1] = 1
                B[2 * i + 1][2 * j] = 1
    return B

def SC_mine(A):
    """SC blowup: copy1 = copy2 = T, BOTH cross directions follow T^op, twin i->i'."""
    n = len(A)
    N = 2 * n
    B = [[0] * N for _ in range(N)]
    for i in range(n):
        B[2 * i][2 * i + 1] = 1
        for j in range(n):
            if i != j and A[i][j]:
                B[2 * i][2 * j] = 1
                B[2 * i + 1][2 * j + 1] = 1
                B[2 * j][2 * i + 1] = 1
                B[2 * j + 1][2 * i] = 1
    return B

def check_tournament(B):
    N = len(B)
    for i in range(N):
        if B[i][i]:
            return False
        for j in range(i + 1, N):
            if B[i][j] + B[j][i] != 1:
                return False
    return True

# ---------------- cycle counting ----------------

def cycles_subset_perm(A, S):
    k = len(S)
    first = S[0]
    c = 0
    for p in itertools.permutations(S[1:]):
        if not A[first][p[0]]:
            continue
        for t in range(k - 2):
            if not A[p[t]][p[t + 1]]:
                break
        else:
            if A[p[-1]][first]:
                c += 1
    return c

def adjbits_of(A):
    n = len(A)
    out = []
    for i in range(n):
        m = 0
        for j in range(n):
            if A[i][j]:
                m |= 1 << j
        out.append(m)
    return out

def cycles_subset_dp(adjb, S):
    k = len(S)
    loc = []
    for v in S:
        m = 0
        av = adjb[v]
        for jj, u in enumerate(S):
            if (av >> u) & 1:
                m |= 1 << jj
        loc.append(m)
    dp = [[0] * k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k, 2):
        row = dp[mask]
        for last in range(k):
            c = row[last]
            if not c:
                continue
            free = loc[last] & ~mask
            while free:
                lb = free & (-free)
                nxt = lb.bit_length() - 1
                dp[mask | lb][nxt] += c
                free ^= lb
    full = (1 << k) - 1
    return sum(dp[full][l] for l in range(1, k) if (loc[l] & 1))

def carriers_perm(A, lengths):
    n = len(A)
    out = {}
    for k in lengths:
        if 3 <= k <= n:
            for S in itertools.combinations(range(n), k):
                c = cycles_subset_perm(A, S)
                if c:
                    out[S] = c
    return out

def ck_double_dp(B, k):
    adjb = adjbits_of(B)
    return sum(cycles_subset_dp(adjb, S) for S in itertools.combinations(range(len(B)), k))

def pair_features_mine(carr):
    items = sorted(carr.items())
    pf = {}
    for a in range(len(items)):
        S1, c1 = items[a]
        l1 = len(S1)
        if c1 >= 2:
            kk = (l1, l1, l1)
            pf[kk] = pf.get(kk, 0) + c1 * (c1 - 1) // 2
        s1 = set(S1)
        for b in range(a + 1, len(items)):
            S2, c2 = items[b]
            sh = len(s1 & set(S2))
            lo, hi = sorted((l1, len(S2)))
            kk = (lo, hi, sh)
            pf[kk] = pf.get(kk, 0) + c1 * c2
    return pf

def conflict_ip(carr):
    items = [(frozenset(S), c) for S, c in carr.items()]
    coeffs = {0: 1}
    def rec(i, used, weight, size):
        for j in range(i, len(items)):
            S, c = items[j]
            if used & S:
                continue
            coeffs[size + 1] = coeffs.get(size + 1, 0) + weight * c
            rec(j + 1, used | S, weight * c, size + 1)
    rec(0, frozenset(), 1, 0)
    return [coeffs.get(t, 0) for t in range(max(coeffs) + 1)]

def strong_comps_mine(A):
    n = len(A)
    reach = []
    for s in range(n):
        seen = 1 << s
        stack = [s]
        while stack:
            v = stack.pop()
            for u in range(n):
                if A[v][u] and not (seen >> u) & 1:
                    seen |= 1 << u
                    stack.append(u)
        reach.append(seen)
    comps = []
    assigned = 0
    for v in range(n):
        if (assigned >> v) & 1:
            continue
        comp = [u for u in range(n) if (reach[v] >> u) & 1 and (reach[u] >> v) & 1]
        for u in comp:
            assigned |= 1 << u
        comps.append(comp)
    return comps

def induced(A, comp):
    return [[A[i][j] for j in comp] for i in comp]

def transitive_mine(n):
    return [[1 if j > i else 0 for j in range(n)] for i in range(n)]

# =====================================================================
t0 = time.time()
w("=== verify_A_H_formula_kps1 — adversarial verification of branch A-H-formula ===")
w("All constructions/counters written independently; doublings use INTERLEAVED order.")
w("")

# ---------------- T0: anchors ----------------
w("--- T0: anchors ---")
C3 = from_bits(3, '101')          # 0->1, 1->2, 2->0
assert scores_mine(C3) == (1, 1, 1), scores_mine(C3)
regT5 = [[0] * 5 for _ in range(5)]
for i in range(5):
    for d in (1, 2):
        regT5[i][(i + d) % 5] = 1
tr3 = transitive_mine(3)
anchors = [
    ("H(C3) backtrack", H_backtrack(C3), 3),
    ("H(C3) dp", H_dp(C3), 3),
    ("H(C3) perm", H_perm(C3), 3),
    ("H(regular T5) dp", H_dp(regT5), 15),
    ("H(regular T5) perm", H_perm(regT5), 15),
    ("H(trans5) dp", H_dp(transitive_mine(5)), 1),
    ("H(D(C3)) dp", H_dp(D_mine(C3)), 45),
    ("H(D(C3)) backtrack", H_backtrack(D_mine(C3)), 45),
    ("H(D(trans3)) dp", H_dp(D_mine(tr3)), 13),
    ("H(SC(trans3)) dp", H_dp(SC_mine(tr3)), 41),
    ("H(K2(C3)) dp", H_dp(K2_mine(C3)), 45),
    ("H(SC(C3)) dp", H_dp(SC_mine(C3)), 45),
]
ok0 = True
for name, got, exp in anchors:
    ok = (got == exp)
    ok0 &= ok
    w(f"  {name:24s} = {got:>6}  expect {exp:>6}  {'OK' if ok else '*** MISMATCH ***'}")
for name, B in (("D(C3)", D_mine(C3)), ("K2(C3)", K2_mine(C3)), ("SC(C3)", SC_mine(C3))):
    w(f"  {name} is a tournament: {check_tournament(B)}")
w(f"T0 verdict: {'PASS' if ok0 else 'FAIL'}")
w("")

# ---------------- T1: the CE pair ----------------
w("--- T1: CE pair n=5 idx4 (0001000000) vs idx6 (0101000000) [claim 1] ---")
A4 = from_bits(5, '0001000000')
A6 = from_bits(5, '0101000000')
w(f"  idx4 scores={scores_mine(A4)} (claimed (1,1,2,3,3));"
  f"  idx6 scores={scores_mine(A6)} (claimed (2,1,1,3,3))")
for label, A in (('idx4', A4), ('idx6', A6)):
    carr = carriers_perm(A, range(3, 6))
    spec = tuple(sum(c for S, c in carr.items() if len(S) == k) for k in (3, 4, 5, 6))
    odd = {S: c for S, c in carr.items() if len(S) % 2 == 1}
    ip = conflict_ip(odd)
    lens = sorted(len(S) for S, c in odd.items() for _ in range(c))
    # any disjoint odd pair?
    disj = any(not (set(S1) & set(S2)) for S1 in odd for S2 in odd if S1 < S2)
    hT_b, hT_d, hT_p = H_backtrack(A), H_dp(A), H_perm(A)
    K = K2_mine(A)
    assert check_tournament(K)
    hK_b, hK_d = H_backtrack(K), H_dp(K)
    w(f"  {label}: spec={spec}  odd-cycle lengths={lens}  disjoint-odd-pair={disj}"
      f"  IP={ip}  H(T) bt/dp/perm={hT_b}/{hT_d}/{hT_p}")
    w(f"        H(K2) backtrack={hK_b}  dp={hK_d}")
w("  raw permutation scan (10! = 3628800) on both doubles:")
t1 = time.time()
hK4_perm = H_perm(K2_mine(A4))
hK6_perm = H_perm(K2_mine(A6))
w(f"    H(K2(idx4)) perm-scan = {hK4_perm}   (claimed 3225)   [{time.time()-t1:.1f}s]")
w(f"    H(K2(idx6)) perm-scan = {hK6_perm}   (claimed 2785)")
# non-isomorphism: c4 differs (3 vs 2) is an iso invariant; also brute check
niso = True
for p in itertools.permutations(range(5)):
    if all(A4[p[i]][p[j]] == A6[i][j] for i in range(5) for j in range(5)):
        niso = False
        break
w(f"  idx4 iso idx6? {not niso}  (must be False; c4 differs 3 vs 2)")
ce_ok = (hK4_perm == 3225 and hK6_perm == 2785 and niso)
w(f"T1 verdict: {'CONFIRMED — IP (and typed IP) does NOT determine H(K2)' if ce_ok else '*** REFUTATION OF CLAIM 1 ***'}")
w("")

# ---------------- T2/T3/T4/T5/T6: 74-class sweep ----------------
w("--- T2-T6: independent 74-class sweep (n=3..6) ---")
t1 = time.time()
expected_small = {  # (n, idx) -> (hT, hD, hK, hS) from sibling P1 table
    (3, 0): (1, 13, 1, 41), (3, 1): (3, 45, 45, 45),
    (4, 0): (1, 95, 1, 629), (4, 1): (3, 189, 45, 633),
    (4, 2): (5, 523, 393, 653), (4, 3): (3, 333, 45, 633),
    (5, 0): (1, 1033, 1, 14937), (5, 1): (3, 1809, 45, 14961),
    (5, 2): (5, 2817, 393, 15109), (5, 3): (3, 1809, 45, 14973),
    (5, 4): (9, 8137, 3225, 15313), (5, 5): (5, 5289, 393, 15109),
    (5, 6): (9, 8145, 2785, 15201), (5, 7): (3, 3561, 45, 14961),
    (5, 8): (11, 11017, 6069, 15461), (5, 9): (15, 12129, 10773, 15261),
    (5, 10): (13, 11625, 8097, 15305), (5, 11): (15, 15505, 15565, 15565),
}
expected_n6_hK = {26: 2025, 46: 356977, 54: 421425, 55: 411513, 30: 210945,
                  51: 386445, 52: 386445, 9: 27513, 42: 118113, 45: 118113}
DATA = []
mismatch_small = 0
mismatch_n6 = 0
for n in (3, 4, 5, 6):
    reps = iso_classes_mine(n)
    w(f"  n={n}: {len(reps)} iso classes")
    for idx, A in enumerate(reps):
        carr = carriers_perm(A, range(3, n + 1))
        spec = tuple(sum(c for S, c in carr.items() if len(S) == k) for k in (3, 4, 5, 6))
        pf = pair_features_mine(carr)
        hT = H_dp(A)
        K = K2_mine(A)
        hK = H_dp(K)
        hD = H_dp(D_mine(A))
        hDop = H_dp(D_mine(transpose(A)))
        hS = H_dp(SC_mine(A))
        comps = strong_comps_mine(A)
        prodK = 1
        prodD = 1
        prodS = 1
        for c in comps:
            if len(c) >= 2:
                sub = induced(A, c)
                prodK *= H_dp(K2_mine(sub))
                prodD *= H_dp(D_mine(sub))
                prodS *= H_dp(SC_mine(sub))
        adjbK = adjbits_of(K)
        c3p = sum(cycles_subset_dp(adjbK, S) for S in itertools.combinations(range(2 * n), 3))
        c5p = sum(cycles_subset_dp(adjbK, S) for S in itertools.combinations(range(2 * n), 5))
        c7p = sum(cycles_subset_dp(adjbK, S) for S in itertools.combinations(range(2 * n), 7)) \
            if 2 * n >= 7 else 0
        DATA.append(dict(n=n, idx=idx, A=A, spec=spec, pf=pf, hT=hT, hK=hK, hD=hD,
                         hDop=hDop, hS=hS, csz=tuple(sorted(len(c) for c in comps)),
                         prodK=prodK, prodD=prodD, prodS=prodS,
                         c3p=c3p, c5p=c5p, c7p=c7p))
        if (n, idx) in expected_small and expected_small[(n, idx)] != (hT, hD, hK, hS):
            mismatch_small += 1
            w(f"  *** P1-TABLE MISMATCH n={n} idx={idx}: mine {(hT, hD, hK, hS)} "
              f"theirs {expected_small[(n, idx)]}")
        if n == 6 and idx in expected_n6_hK and expected_n6_hK[idx] != hK:
            mismatch_n6 += 1
            w(f"  *** n=6 hK MISMATCH idx={idx}: mine {hK} theirs {expected_n6_hK[idx]}")
w(f"  sweep done [{time.time()-t1:.1f}s]; "
  f"P1-table spot mismatches: {mismatch_small} (n<=5 full) + {mismatch_n6} (n=6 spots)")
w("")

# T2: spectrum determination
groups = {}
for d in DATA:
    groups.setdefault(d['spec'], []).append(d)
multi = [g for g in groups.values() if len(g) > 1]
def broken(valkey):
    return [g for g in groups.values() if len({d[valkey] for d in g}) > 1]
bK, bD, bS = broken('hK'), broken('hD'), broken('hS')
w(f"T2: spectrum groups: {len(groups)} distinct (claimed 32), {len(multi)} multi (claimed 22)")
w(f"    spectrum -> H(K2): {len(bK)} broken (claimed 0)"
  + ("" if not bK else "  *** BROKEN: " + str([(d['n'], d['idx'], d['hK']) for g in bK for d in g])))
w(f"    spectrum -> H(D): {len(bD)} broken (claimed 22); spectrum -> H(SC): {len(bS)} broken (claimed 9)")
w(f"T2 verdict: {'CONFIRMED' if len(groups)==32 and len(multi)==22 and not bK and len(bD)==22 and len(bS)==9 else 'DISCREPANCY'}")
w("")

# T3: product law
badK = [d for d in DATA if d['prodK'] != d['hK']]
okD = sum(1 for d in DATA if d['prodD'] == d['hD'])
okS = sum(1 for d in DATA if d['prodS'] == d['hS'])
nstrong = sum(1 for d in DATA if len(d['csz']) == 1)
idx26 = [d for d in DATA if d['n'] == 6 and d['idx'] == 26][0]
w(f"T3: product law H(K2): {len(DATA)-len(badK)}/{len(DATA)} (claimed 74/74)"
  + ("" if not badK else "  *** FAILS: " + str([(d['n'], d['idx']) for d in badK])))
w(f"    product law H(D): {okD}/74 (claimed 43);  H(SC): {okS}/74 (claimed 43); "
  f"strong classes: {nstrong} (claimed 43)")
w(f"    n=6 idx26 comps={idx26['csz']} H(K2)={idx26['hK']} (claimed 2025=45^2)")
w(f"T3 verdict: {'CONFIRMED' if not badK and okD==43 and okS==43 and nstrong==43 and idx26['hK']==2025 else 'DISCREPANCY'}")
w("")

# T4: lift laws
ok3 = sum(1 for d in DATA if d['c3p'] == 8 * d['spec'][0])
ok5 = sum(1 for d in DATA if d['c5p'] == 32 * d['spec'][2] + 32 * d['spec'][1] + 6 * d['spec'][0])
ok7 = 0
n_p341 = 0
fails7 = []
okI = 0
for d in DATA:
    c3, c4, c5, c6 = d['spec']
    pf = d['pf']
    p331, p332 = pf.get((3, 3, 1), 0), pf.get((3, 3, 2), 0)
    p341, p342, p343 = pf.get((3, 4, 1), 0), pf.get((3, 4, 2), 0), pf.get((3, 4, 3), 0)
    if p341 > 0:
        n_p341 += 1
    pred = (192 * c6 + 80 * c5 + 8 * c4 + 64 * p331 + 48 * p332 + 64 * p341 + 32 * p342)
    if d['c7p'] == pred:
        ok7 += 1
    else:
        fails7.append((d['n'], d['idx'], d['c7p'], pred))
    if p343 == 2 * p332:
        okI += 1
w(f"T4: c3'(K2) == 8 c3: {ok3}/74;   c5'(K2) == 32c5+32c4+6c3: {ok5}/74")
w(f"    c7'(K2) == 128c7+192c6+80c5+8c4+64p331+48p332+64p341+32p342: {ok7}/74"
  + ("" if not fails7 else "  *** FAILS: " + str(fails7[:5])))
w(f"    NOTE: c7 of base == 0 for ALL n<=6, so the 128c7 coefficient is UNTESTED here.")
w(f"    classes with p341 > 0 (testing the 64*p341 coefficient): {n_p341}")
w(f"    p(3,4,3) == 2*p(3,3,2): {okI}/74")
w(f"T4 verdict: {'CONFIRMED' if ok3==74 and ok5==74 and ok7==74 and okI==74 else 'DISCREPANCY'}")
w("")

# T5: congruences
okm8 = sum(1 for d in DATA if (d['hK'] - (2 * d['hT'] - 1)) % 8 == 0)
cen16 = {}
for d in DATA:
    r = (d['hK'] - (2 * d['hT'] - 1)) % 16
    cen16[r] = cen16.get(r, 0) + 1
okSC4 = sum(1 for d in DATA if d['hS'] % 4 == 1)
w(f"T5: H(K2) == 2H(T)-1 (mod 8): {okm8}/74 (claimed 74)")
w(f"    (H(K2)-(2H(T)-1)) mod 16 census: {dict(sorted(cen16.items()))} (claimed {{0:39, 8:35}})")
w(f"    H(SC) == 1 (mod 4): {okSC4}/74 (claimed 74)")
w(f"T5 verdict: {'CONFIRMED' if okm8==74 and cen16=={0:39,8:35} and okSC4==74 else 'DISCREPANCY'}")
w("")

# T6: D op-asymmetry + odd-n mod laws
asym = [d for d in DATA if d['hD'] != d['hDop']]
d41 = [d for d in DATA if d['n'] == 4 and d['idx'] == 1][0]
n5mod8 = [d['hD'] % 8 for d in DATA if d['n'] == 5]
n3mod8 = [d['hD'] % 8 for d in DATA if d['n'] == 3]
w(f"T6: H(D(T)) != H(D(T^op)): {len(asym)}/74 (claimed 50)")
w(f"    n=4 idx1: H(D)={d41['hD']} H(D(op))={d41['hDop']} (claimed 189 vs 333)")
w(f"    n=5 H(D) mod 8 values: {sorted(set(n5mod8))} (claimed all 1);"
  f"  n=3 mod 8: {sorted(set(n3mod8))} (claimed all 5)")
ok6 = (len(asym) == 50 and d41['hD'] == 189 and d41['hDop'] == 333
       and set(n5mod8) == {1} and set(n3mod8) == {5})
w(f"T6 verdict: {'CONFIRMED' if ok6 else 'DISCREPANCY'}")
w("")

# ---------------- T7: transitive towers ----------------
w("--- T7: transitive towers (independent constructions + fresh DP) ---")
seqD, seqS = [], []
for n in range(3, 8):
    tr = transitive_mine(n)
    t2 = time.time()
    hD = H_dp(D_mine(tr))
    hS = H_dp(SC_mine(tr))
    hK = H_dp(K2_mine(tr))
    seqD.append(hD); seqS.append(hS)
    w(f"  n={n}: H(D)={hD}  H(SC)={hS}  H(K2)={hK}   [{time.time()-t2:.1f}s]")
expD = [13, 95, 1033, 15611, 313285]
expS = [41, 629, 14937, 513669, 24104937]
w(f"  H(D(trans_n)) n=3..7: {seqD} (claimed {expD})")
w(f"  H(SC(trans_n)) n=3..7: {seqS} (claimed {expS})")
w(f"T7 verdict: {'CONFIRMED' if seqD==expD and seqS==expS else 'DISCREPANCY'}")
w("")

# ---------------- T8: profile/IP spot check on strong pair idx42/idx45 ----------------
w("--- T8: spec (6,7,6,3) strong pair n=6 idx42 vs idx45: full odd profile + IP of Omega(K2) ---")
w("    claimed: both profile'=(48,452,3016,10148,9856), IP'=(1,23520,14664,1520,16), H(K2)=118113")
sigs = []
for idx in (42, 45):
    d = [x for x in DATA if x['n'] == 6 and x['idx'] == idx][0]
    B = K2_mine(d['A'])
    adjb = adjbits_of(B)
    t2 = time.time()
    carr = []
    prof = []
    for k in (3, 5, 7, 9, 11):
        tot = 0
        for S in itertools.combinations(range(12), k):
            c = cycles_subset_dp(adjb, S)
            if c:
                m = 0
                for v in S:
                    m |= 1 << v
                carr.append((m, c, k))
                tot += c
        prof.append(tot)
    i1 = sum(c for _, c, _ in carr)
    i2 = 0
    for a in range(len(carr)):
        ma, ca, _ = carr[a]
        for b in range(a + 1, len(carr)):
            mb, cb, _ = carr[b]
            if not (ma & mb):
                i2 += ca * cb
    T3c = [(m, c) for m, c, k in carr if k == 3]
    T5c = [(m, c) for m, c, k in carr if k == 5]
    i3 = i4 = 0
    for a in range(len(T3c)):
        ma, ca = T3c[a]
        for b in range(a + 1, len(T3c)):
            mb, cb = T3c[b]
            if ma & mb:
                continue
            mab = ma | mb
            for c_ in range(b + 1, len(T3c)):
                mc, cc = T3c[c_]
                if mab & mc:
                    continue
                i3 += ca * cb * cc
                mabc = mab | mc
                for e in range(c_ + 1, len(T3c)):
                    me, ce = T3c[e]
                    if not (mabc & me):
                        i4 += ca * cb * cc * ce
            for mf, cf in T5c:
                if not (mab & mf):
                    i3 += ca * cb * cf
    ocf = 1 + 2 * i1 + 4 * i2 + 8 * i3 + 16 * i4
    sigs.append((tuple(prof), (i1, i2, i3, i4)))
    w(f"  idx{idx}: profile'={tuple(prof)} IP'=(1,{i1},{i2},{i3},{i4}) "
      f"I'(2)={ocf} H(K2)={d['hK']} match={ocf == d['hK']}   [{time.time()-t2:.1f}s]")
exp_prof = (48, 452, 3016, 10148, 9856)
exp_ip = (23520, 14664, 1520, 16)
ok8 = (sigs[0] == sigs[1] == (exp_prof, exp_ip))
w(f"T8 verdict: {'CONFIRMED (identical, matches claim)' if ok8 else 'DISCREPANCY: ' + str(sigs)}")
w("")

# cross-validation of the two cycle counters on a 10-vertex double
w("--- cross-validation: perm-method vs subset-DP cycle counts on K2(idx4), K2(idx6) ---")
for label, A in (('idx4', A4), ('idx6', A6)):
    B = K2_mine(A)
    adjb = adjbits_of(B)
    for k in (3, 5, 7, 9):
        dpv = sum(cycles_subset_dp(adjb, S) for S in itertools.combinations(range(10), k))
        pmv = sum(cycles_subset_perm(B, S) for S in itertools.combinations(range(10), k))
        w(f"  {label} c{k}'(K2): dp={dpv} perm={pmv} {'OK' if dpv == pmv else '*** MISMATCH ***'}")
w("")
w(f"=== done in {time.time()-t0:.1f}s ===")
OUT.close()
