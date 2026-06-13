# ocf_neg_structured_kpo2_verify.py
# ADVERSARIAL VERIFIER, thread C (HYP-2380): claims C5 (strong-component
# factorization), C8 (n=9 triple-triangle counterexample to I(-2)=H-4a1),
# C9 (Paley_7 / circulant census rows), C10 (sign anatomy n=7/8). FRESH code.
# Reduced sample sizes vs worker (runtime guard): 300 random n=7, 80 random n=8.
# H by DFS path extension (validated vs permutation brute force in census script;
# revalidated here on n=7 randoms x20 and the n=9 construction by perm brute force).

import itertools, random
from collections import Counter

random.seed(424242)

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

def out_from_adj(adj):
    n = len(adj)
    return [sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)]

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

def alpha_list(cycles, n):
    """full alpha vector [1, a1, a2, ...] by recursive disjoint-collection count."""
    kmax = n // 3
    a = [0] * (kmax + 1)
    a[0] = 1
    L = cycles
    def rec(idx, mask, size):
        for i in range(idx, len(L)):
            if L[i] & mask == 0:
                a[size + 1] += 1
                if size + 1 < kmax:
                    rec(i + 1, mask | L[i], size + 1)
    rec(0, 0, 0)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def I_eval(a, x):
    return sum(c * x ** k for k, c in enumerate(a))

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

def ham_perm(n, adj):
    c = 0
    for p in itertools.permutations(range(n)):
        if all(adj[p[i]][p[i + 1]] for i in range(n - 1)):
            c += 1
    return c

def strong_components(n, adj):
    """SCCs via reachability bitmask closure (fresh, no Tarjan)."""
    reach = [sum(1 << j for j in range(n) if adj[i][j]) | (1 << i) for i in range(n)]
    changed = True
    while changed:
        changed = False
        for i in range(n):
            r = reach[i]
            m = r
            while m:
                b = m & -m
                r |= reach[b.bit_length() - 1]
                m ^= b
            if r != reach[i]:
                reach[i] = r
                changed = True
    comps = []
    assigned = [False] * n
    for i in range(n):
        if assigned[i]:
            continue
        comp = [j for j in range(n) if (reach[i] >> j) & 1 and (reach[j] >> i) & 1]
        for j in comp:
            assigned[j] = True
        comps.append(comp)
    return comps

def induced(adj, verts):
    k = len(verts)
    return [[adj[verts[i]][verts[j]] for j in range(k)] for i in range(k)]

def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

fails = 0

# ---------------- PART A: strong-component factorization (C5) ----------------
print("PART A: I(Omega,x) factorization over strong components (fresh check)")
nontrivial = 0
tested = 0
for n, count in ((7, 300), (8, 80)):
    pairs = pairs_list(n)
    for _ in range(count):
        mask = random.getrandbits(len(pairs))
        adj, out_ = adj_from_mask(n, mask, pairs)
        cycles = odd_cycle_supports(n, adj)
        a_direct = alpha_list(cycles, n)
        comps = strong_components(n, adj)
        poly = [1]
        for comp in comps:
            sub = induced(adj, comp)
            k = len(comp)
            ac = alpha_list(odd_cycle_supports(k, sub), k) if k >= 3 else [1]
            poly = polymul(poly, ac)
        if poly != a_direct:
            print(f"  FAIL n={n} mask={mask}: direct {a_direct} vs product {poly}")
            fails += 1
        if len(comps) > 1:
            nontrivial += 1
        tested += 1
        # 2-adic mirror + OCF at n=7,8 with full alpha vector
        H = ham_dfs(n, out_)
        if I_eval(a_direct, 2) != H:
            print(f"  OCF FAIL n={n} mask={mask}"); fails += 1
        if I_eval(a_direct, -2) != H - 4 * a_direct[1]:
            print(f"  exact-mirror FAIL n={n} mask={mask}"); fails += 1
        if (H + I_eval(a_direct, -2)) % 8 != 2:
            print(f"  mod8 FAIL n={n} mask={mask}"); fails += 1
print(f"  tested {tested} random tournaments (300 n=7 + 80 n=8); "
      f"{nontrivial} with nontrivial condensation; failures so far: {fails}")

# ---------------- PART B: n=7 sign anatomy (C10) ----------------
print("\nPART B: n=7 sign anatomy on fresh 600-sample (worker used 2000)")
pairs7 = pairs_list(7)
posvecs = Counter(); neg = 0; pos = 0; zero = 0
permcheck = 0
for t in range(600):
    mask = random.getrandbits(21)
    adj, out_ = adj_from_mask(7, mask, pairs7)
    a = alpha_list(odd_cycle_supports(7, adj), 7)
    Im2 = I_eval(a, -2)
    H = ham_dfs(7, out_)
    if I_eval(a, 2) != H:
        print(f"  OCF FAIL mask={mask}"); fails += 1
    if t < 20 and ham_perm(7, adj) != H:
        print(f"  DFS-vs-perm FAIL mask={mask}"); fails += 1
    if Im2 < 0: neg += 1
    elif Im2 == 0: zero += 1; print(f"  ZERO?! mask={mask}"); fails += 1
    else:
        pos += 1
        posvecs[(tuple(a), Im2)] += 1
print(f"  n=7 sample: neg {neg}/600 ({100*neg/600:.1f}%), zero {zero}, pos {pos}")
print(f"  positive cases (alpha, I(-2)) -> count: {dict(posvecs)}")
print("  worker claims positives only alpha in {[1],[1,2,1],[1,3,2]} with I(-2) in {1,1,3}")

# ---------------- PART C: structured families (C8, C9, C10) ----------------
print("\nPART C: structured families")

def block_tournament(blocks):
    """blocks: list of adjacency matrices; all arcs earlier-block -> later-block."""
    sizes = [len(b) for b in blocks]
    n = sum(sizes)
    adj = [[0] * n for _ in range(n)]
    off = []
    s = 0
    for sz in sizes:
        off.append(s); s += sz
    for bi, b in enumerate(blocks):
        for i in range(len(b)):
            for j in range(len(b)):
                adj[off[bi] + i][off[bi] + j] = b[i][j]
        for bj in range(bi + 1, len(blocks)):
            for i in range(sizes[bi]):
                for j in range(sizes[bj]):
                    adj[off[bi] + i][off[bj] + j] = 1
    return adj

C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
# unique strong 4-tournament: triangle 0->1->2->0, 3->0, 1->3, 2->3
S4 = [[0, 1, 0, 0], [0, 0, 1, 1], [1, 0, 0, 1], [1, 0, 0, 0]]

# n=9 triple triangle (C8)
adj9 = block_tournament([C3, C3, C3])
out9 = out_from_adj(adj9)
cyc9 = odd_cycle_supports(9, adj9)
a9 = alpha_list(cyc9, 9)
H9 = ham_dfs(9, out9)
H9p = ham_perm(9, adj9)
Im2_9 = I_eval(a9, -2)
print(f"  C3+C3+C3 (n=9): alpha={a9}, H(dfs)={H9}, H(perm)={H9p}, "
      f"I(2)={I_eval(a9,2)}, I(-2)={Im2_9}, H-4a1={H9-4*a9[1]}")
ok9 = (a9 == [1, 3, 3, 1] and H9 == 27 == H9p == I_eval(a9, 2) and Im2_9 == -1
       and H9 - 4 * a9[1] == 15)
print(f"  C8 counterexample check: {'CONFIRMED' if ok9 else 'FAIL'}"
      f"  ((1+x)^3 -> I(-2)=-1 vs H-4a1=15)")
if not ok9: fails += 1

# S4 alone, S4+S4 (n=8), S4+C3 (n=7), C3+C3 (n=6)
for name, blocks, expect_alpha in (("S4", [S4], [1, 2]),
                                   ("S4+S4 (n=8)", [S4, S4], [1, 4, 4]),
                                   ("S4+C3 (n=7)", [S4, C3], [1, 3, 2]),
                                   ("C3+C3 (n=6)", [C3, C3], [1, 2, 1])):
    adjB = block_tournament(blocks)
    nB = len(adjB)
    aB = alpha_list(odd_cycle_supports(nB, adjB), nB)
    HB = ham_dfs(nB, out_from_adj(adjB))
    disc = aB[1] ** 2 - 4 * (aB[2] if len(aB) > 2 else 0)
    print(f"  {name}: alpha={aB} (expect {expect_alpha}), H={HB}, I(2)={I_eval(aB,2)}, "
          f"I(-2)={I_eval(aB,-2)}, I(-1)={I_eval(aB,-1)}, disc={disc}")
    if aB != expect_alpha or I_eval(aB, 2) != HB:
        print("  ^^ FAIL"); fails += 1

# n=7 circulants (C9): S subset of {1..6}, one of each pair {k,7-k}
print("\n  n=7 circulants T_S (i->i+s mod 7 for s in S):")
for s1 in (1, 6):
    for s2 in (2, 5):
        for s3 in (3, 4):
            S = (s1, s2, s3)
            adjC = [[0] * 7 for _ in range(7)]
            for i in range(7):
                for s in S:
                    adjC[i][(i + s) % 7] = 1
            aC = alpha_list(odd_cycle_supports(7, adjC), 7)
            HC = ham_dfs(7, out_from_adj(adjC))
            paley = set(S) in ({1, 2, 4}, {3, 5, 6})
            print(f"    S={S}: alpha={aC}, H={HC}, I(2)={I_eval(aC,2)}, "
                  f"I(-2)={I_eval(aC,-2)}, {'PALEY' if paley else 'non-Paley'}")
            if I_eval(aC, 2) != HC:
                print("    ^^ OCF FAIL"); fails += 1

# transitive n=7
adjT = [[1 if i < j else 0 for j in range(7)] for i in range(7)]
aT = alpha_list(odd_cycle_supports(7, adjT), 7)
HT = ham_dfs(7, out_from_adj(adjT))
print(f"  transitive n=7: alpha={aT}, H={HT}, I(-2)={I_eval(aT,-2)}")

print(f"\nTOTAL FAILURES: {fails}")
print("done")
