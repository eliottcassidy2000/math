# ocf_neg_census_kpo2_verify.py
# ADVERSARIAL VERIFIER for session kind-pasteur-2026-06-10-S2, thread C (HYP-2380).
# FRESH code, written independently of the worker's scripts.
# Methods deliberately DIFFERENT from the worker:
#   * H counted by DFS path extension (worker used Held-Karp DP), cross-checked
#     against plain brute force over all n! permutations at n<=5 (all) and n=6 (sample).
#   * Directed odd cycles enumerated by (subset) x (cyclic order with fixed start = min
#     vertex of the subset), i.e. brute force over (k-1)! cyclic permutations -- NOT a
#     DFS/anchored-walk enumeration, and NOT vertex-set counting (MISTAKE-023 guard:
#     each DIRECTED cycle counted separately; both directions of a support count).
# Exact integer arithmetic throughout.
#
# Checks (full labeled census n=3..6):
#   I(Omega,2) == H  (THM-002 / OCF)
#   I(-2) odd (never 0); H + I(-2) == 2 mod 8; H - I(-2) == 4*alpha_1 mod 16;
#   I(-2) == H - 4*alpha_1 exactly (n<=8, alpha_3=0); H == 1 + 2*alpha_1 mod 4 (GS 7.1);
#   disc = alpha_1^2 - 4*alpha_2 >= 0; alpha_3 == 0 for n<=6;
#   sign tallies of I(-2); I(-1) value distribution (min/max, count of -2, count of +1);
#   per-iso-class table (canonical form via score-stabilizer minimization),
#   class counts must be 2,4,12,56 (A000568); SC = self-converse flag.
#   MISTAKE-023 regression: n=5 bits=40 decoded under both bit conventions.

import itertools, random, sys
from collections import Counter

random.seed(20260610)

def pairs_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_mask(n, mask, pairs, bit_means_ij=True):
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    for b, (i, j) in enumerate(pairs):
        if ((mask >> b) & 1) == (1 if bit_means_ij else 0):
            adj[i][j] = 1; out[i] |= 1 << j
        else:
            adj[j][i] = 1; out[j] |= 1 << i
    return adj, out

def odd_cycle_supports(n, adj):
    """All DIRECTED odd cycles, as support bitmasks (one entry per directed cycle)."""
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

def alpha_vector(cycles):
    """[alpha_1, alpha_2, alpha_3] by brute pair/triple disjointness."""
    a1 = len(cycles); a2 = 0; a3 = 0
    L = cycles
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if L[i] & L[j] == 0:
                a2 += 1
                u = L[i] | L[j]
                for k in range(j + 1, len(L)):
                    if L[k] & u == 0:
                        a3 += 1
    return a1, a2, a3

def I_eval(a, x):
    a1, a2, a3 = a
    return 1 + a1 * x + a2 * x * x + a3 * x ** 3

def ham_dfs(n, out):
    full = (1 << n) - 1
    total = 0
    stack = [(v, 1 << v) for v in range(n)]
    pop = stack.pop; push = stack.append
    while stack:
        v, vis = pop()
        if vis == full:
            total += 1
            continue
        avail = out[v] & ~vis
        while avail:
            b = avail & -avail
            push((b.bit_length() - 1, vis | b))
            avail ^= b
    return total

def ham_perm(n, adj):
    c = 0
    for p in itertools.permutations(range(n)):
        if all(adj[p[i]][p[i + 1]] for i in range(n - 1)):
            c += 1
    return c

def canon(n, adj, pairs):
    """Canonical labeled mask: min over relabelings that sort the score sequence."""
    scores = [sum(adj[v]) for v in range(n)]
    groups = {}
    for v in range(n):
        groups.setdefault(scores[v], []).append(v)
    keys = sorted(groups)
    best = None
    for parts in itertools.product(*[itertools.permutations(groups[k]) for k in keys]):
        perm = [v for part in parts for v in part]  # new position -> old vertex
        m = 0
        b = 0
        for i in range(n):
            pi = perm[i]
            for j in range(i + 1, n):
                if adj[pi][perm[j]]:
                    m |= 1 << b
                b += 1
        if best is None or m < best:
            best = m
    return best

def transpose(adj):
    n = len(adj)
    return [[adj[j][i] for j in range(n)] for i in range(n)]

print("=" * 78)
print("FRESH CENSUS VERIFICATION (verifier kpo2) -- full labeled n=3..6")
print("bit convention: bit b set => i->j for b-th pair (i,j), i<j, lex order")
print("=" * 78)

EXPECTED_CLASSES = {3: 2, 4: 4, 5: 12, 6: 56}

grand_fail = 0

for n in range(3, 7):
    pairs = pairs_list(n)
    nb = len(pairs)
    N = 1 << nb
    fails = Counter()
    sign_neg = sign_zero = sign_pos = 0
    Im1_counter = Counter()
    classes = {}   # canon -> [orbit, rep_mask, data]
    sample_perm_check = set()
    if n == 6:
        sample_perm_check = set(random.sample(range(N), 300))
    for mask in range(N):
        adj, out = adj_from_mask(n, mask, pairs)
        cyc = odd_cycle_supports(n, adj)
        a = alpha_vector(cyc)
        H = ham_dfs(n, out)
        if n <= 5 or mask in sample_perm_check:
            if ham_perm(n, adj) != H:
                fails['dfs_vs_perm'] += 1
        I2 = I_eval(a, 2); Im2 = I_eval(a, -2); Im1 = I_eval(a, -1); I1 = I_eval(a, 1)
        if I2 != H: fails['OCF I(2)==H'] += 1
        if Im2 % 2 == 0: fails['I(-2) odd'] += 1
        if (H + Im2) % 8 != 2: fails['H+I(-2)==2 mod 8'] += 1
        if (H - Im2 - 4 * a[0]) % 16 != 0: fails['H-I(-2)==4a1 mod 16'] += 1
        if Im2 != H - 4 * a[0]: fails['I(-2)==H-4a1 exact (n<=8)'] += 1
        if (H - 1 - 2 * a[0]) % 4 != 0: fails['H==1+2a1 mod 4 (GS7.1)'] += 1
        if a[0] * a[0] - 4 * a[1] < 0: fails['disc>=0'] += 1
        if a[2] != 0: fails['alpha_3==0 at n<=6'] += 1
        if Im2 < 0: sign_neg += 1
        elif Im2 == 0: sign_zero += 1
        else: sign_pos += 1
        Im1_counter[Im1] += 1
        c = canon(n, adj, pairs)
        if c not in classes:
            cT = canon(n, transpose(adj), pairs)
            sc = tuple(sorted(sum(adj[v]) for v in range(n)))
            c3 = sum(1 for m in cyc if bin(m).count('1') == 3)
            c5 = sum(1 for m in cyc if bin(m).count('1') == 5)
            classes[c] = [0, mask, sc, c3, c5, a, H, Im2, Im1, I1, (c == cT)]
        classes[c][0] += 1
    tot = N
    print(f"\n--- n={n}: {N} labeled tournaments, {len(classes)} iso classes "
          f"(expected {EXPECTED_CLASSES[n]}) ---")
    if len(classes) != EXPECTED_CLASSES[n]:
        fails['class count'] += 1
    orbsum = sum(v[0] for v in classes.values())
    if orbsum != N:
        fails['orbit sizes sum'] += 1
    for k, v in sorted(fails.items()):
        print(f"  CHECK FAIL {k}: {v} cases")
        grand_fail += v
    if not fails:
        print("  ALL CHECKS PASS (OCF, 2-adic mirror, parity, GS7.1 congruence, disc>=0,"
              " dfs-vs-perm H, class count)")
    print(f"  I(-2) sign census: neg {sign_neg} ({100.0*sign_neg/tot:.1f}%), "
          f"zero {sign_zero}, pos {sign_pos} ({100.0*sign_pos/tot:.1f}%)")
    print(f"  I(-1): min {min(Im1_counter)}, max {max(Im1_counter)}, "
          f"#(I(-1)=-2) {Im1_counter.get(-2,0)}, #(I(-1)=+1) {Im1_counter.get(1,0)}, "
          f"#(I(-1) in {{0,+-1}}) {sum(c for v,c in Im1_counter.items() if v in (0,1,-1))}")
    if n <= 6:
        print(f"  I(-1) full distribution: {dict(sorted(Im1_counter.items()))}")
    print(f"  per-iso-class table (canon, orbit, scores, c3, c5, alpha, H, I(-2), I(-1), I(1), SC):")
    for c, v in sorted(classes.items()):
        orbit, rep, sc, c3, c5, a, H, Im2, Im1, I1, scflag = v
        print(f"    canon={c:5d} orbit={orbit:5d} scores={sc} c3={c3} c5={c5} "
              f"alpha=({a[0]},{a[1]}) H={H:4d} I(-2)={Im2:5d} I(-1)={Im1:4d} I(1)={I1:4d} "
              f"SC={'Y' if scflag else 'N'}")

# MISTAKE-023 regression: n=5 bits=40, both conventions
print("\n--- MISTAKE-023 regression (n=5, bits=40) ---")
pairs = pairs_list(5)
for conv in (True, False):
    adj, out = adj_from_mask(5, 40, pairs, bit_means_ij=conv)
    cyc = odd_cycle_supports(5, adj)
    a = alpha_vector(cyc)
    H = ham_dfs(5, out)
    full = sum(1 for m in cyc if bin(m).count('1') == 5)
    print(f"  convention bit=>{'i->j' if conv else 'j->i'}: alpha_1={a[0]}, "
          f"directed 5-cycles on full set={full}, H={H}, I(2)={I_eval(a,2)}")
print("  expected (MISTAKES.md): alpha_1=7, three directed 5-cycles on one 5-set, H=15")

print(f"\nGRAND TOTAL CHECK FAILURES: {grand_fail}")
print("done")
