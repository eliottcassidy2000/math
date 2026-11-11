"""a1_visit_matrix.py (AUDIT) - items 1, 2 of the audit brief.

1. Visit matrix by plain depth-first enumeration of walks for ALL classes of order <= 5 and all ordered port pairs,
   and for a sample of order 6;  H(K[A, B]) = sum N_ij P_i(A) P_j(B) on random modules (H of the substituted
   tournament by an independent subset DP in C, P_i by a DP cross-checked against a brute force over permutations).
2. Bound N_ij = 0 unless i, j <= m - 1 (and whether the bound is attained).
3. Two-visit kernels: the factorisation N22 H = L(A) L'(B) + det G H(A) H(B) and the (al, be, ga, de) form.
4. Proposition 2.3: for random pairs (B0, O) with state (r, M) computed from the DEFINITIONS of THM-5130
   (r = H(B0) - W(O), M = U(O)), the kernel K = C(B0, O; x => y) has at most two visits, G with
   (al, be, ga, de) = (r, 0, 0, M), det G = r M; includes O empty (then M = 0).
5. The Gaussian kernel and the first hyperbolic kernel."""
import random, sys, itertools
from aud_common import *

rng = random.Random(20261111)
out_lines = []


def P(msg):
    print(msg)
    out_lines.append(msg)
    sys.stdout.flush()


# ---- self-tests of the tools
for trial in range(40):
    n = rng.randint(1, 6)
    T = random_tournament(n, rng)
    pk = P_tuple(T, 4)
    for k in range(1, 5):
        assert pk[k] == P_tuple_brute(T, k), (T, k)
    assert pk[1] == ham_brute(T) == ham_py(T)
P("self-test: P_k DP == brute force over permutations (k <= 4, n <= 6, 40 random tournaments); H DP == brute force")

# ---- 1. all classes of order <= 5, all ordered port pairs; sample at order 6
total_cases = 0
fail = 0
bound_viol = 0
attained = {}
three_plus = 0
kernels_checked = 0
census = {}
mats = {}
for m in (3, 4, 5, 6):
    cls = read_classes(m)
    if m == 6:
        pairs_all = [(c, x, y) for (c, h) in cls for x in range(m) for y in range(m) if x != y]
        sample = rng.sample(pairs_all, 150)
    else:
        sample = [(c, x, y) for (c, h) in cls for x in range(m) for y in range(m) if x != y]
    for (code, x, y) in sample:
        K = decode(code, m)
        N = visit_matrix_dfs(K, x, y)
        Npos = {k: v for k, v in N.items() if k[0] >= 1 and k[1] >= 1}
        mats[(m, code, x, y)] = Npos
        maxi = max(k[0] for k in Npos)
        maxj = max(k[1] for k in Npos)
        if maxi > m - 1 or maxj > m - 1:
            bound_viol += 1
        if maxi == m - 1 or maxj == m - 1:
            attained[m] = attained.get(m, 0) + 1
        if maxi >= 3 or maxj >= 3:
            three_plus += 1
        kernels_checked += 1
        # random modules: sizes chosen so that the total order stays <= 15
        ntest = 2 if m <= 5 else 2
        for t in range(ntest):
            na = rng.randint(1, 5 if m <= 4 else 4)
            nb = rng.randint(1, 5 if m <= 4 else 4)
            A = random_tournament(na, rng)
            B = random_tournament(nb, rng)
            pa = P_tuple(A, max(maxi, 1))
            pb = P_tuple(B, max(maxj, 1))
            pred = sum(v * pa[i] * pb[j] for (i, j), v in Npos.items())
            T = substitute(K, [x, y], [A, B])
            mats.setdefault("jobs", []).append((T, pred, (m, code, x, y, na, nb, maxi, maxj)))
jobs = mats.pop("jobs")
Hs = ham_c([j[0] for j in jobs])
tp_cases = 0
for (T, pred, info), h in zip(jobs, Hs):
    total_cases += 1
    if info[6] >= 3 or info[7] >= 3:
        tp_cases += 1
    if h != pred:
        fail += 1
        if fail < 5:
            P("   FAILURE %s: H = %d, predicted %d" % (info, h, pred))
P("1. visit matrix (walks enumerated by DFS): kernels (class, ordered ports) checked %d (all of orders 3, 4, 5; 150 sampled at order 6);"
  % kernels_checked)
P("   kernels with a port visited >= 3 times: %d; module cases %d (of which %d on kernels with >= 3 visits), failures %d"
  % (three_plus, total_cases, tp_cases, fail))
P("2. bound i, j <= m - 1: violations %d; kernels where some port is visited exactly m - 1 times, by order: %s"
  % (bound_viol, sorted(attained.items())))

# N_0j is NOT zero in general (wording of Theorem 2.1)
K = decode(5, 3)
N = visit_matrix_dfs(K, 0, 1)
P("   remark: for C_3 (code 5), ports (0, 1), the walk counts with i = 0 or j = 0 are %s (not zero; they are simply not used)"
  % sorted((k, v) for k, v in N.items() if k[0] == 0 or k[1] == 0))

# ---- 3. factorisation for two-visit kernels
fcases = 0
ffail = 0
form_fail = 0
jobs = []
for key, Npos in mats.items():
    m, code, x, y = key
    if max(k[0] for k in Npos) > 2 or max(k[1] for k in Npos) > 2:
        continue
    n11, n12, n21, n22 = (Npos.get((1, 1), 0), Npos.get((1, 2), 0), Npos.get((2, 1), 0), Npos.get((2, 2), 0))
    det = n11 * n22 - n12 * n21
    al = n11 - 2 * n12 - 2 * n21 + 4 * n22
    be = n21 - 2 * n22
    ga = n12 - 2 * n22
    de = n22
    assert al * de - be * ga == det
    if rng.random() < 0.5:
        continue
    K = decode(code, m)
    na, nb = rng.randint(1, 4), rng.randint(1, 4)
    A = random_tournament(na, rng)
    B = random_tournament(nb, rng)
    pa = P_tuple(A, 2)
    pb = P_tuple(B, 2)
    hA, hB = pa[1], pb[1]
    piA, piB = pa[2] + 2 * hA, pb[2] + 2 * hB
    jobs.append((substitute(K, [x, y], [A, B]), (n11, n12, n21, n22, det, al, be, ga, de, hA, hB, pa[2], pb[2], piA, piB)))
    # transitive modules
    a, b = rng.randint(1, 4), rng.randint(1, 4)
    jobs.append((substitute(K, [x, y], [transitive(a), transitive(b)]), (n11, n12, n21, n22, det, al, be, ga, de, 1, 1, 2 ** a - 2, 2 ** b - 2, 2 ** a, 2 ** b)))
Hs = ham_c([j[0] for j in jobs])
for (T, (n11, n12, n21, n22, det, al, be, ga, de, hA, hB, p2A, p2B, piA, piB)), h in zip(jobs, Hs):
    fcases += 1
    LA = n22 * p2A + n12 * hA
    LB = n22 * p2B + n21 * hB
    if n22 * h != LA * LB + det * hA * hB:
        ffail += 1
    if h != al * hA * hB + be * piA * hB + ga * hA * piB + de * piA * piB:
        form_fail += 1
P("3. two-visit kernels: N22 H = L(A) L'(B) + det G H(A) H(B): cases %d, failures %d; form al hh' + be pi h' + ga h pi' + de pi pi': failures %d"
  % (fcases, ffail, form_fail))


# ---- 4. Proposition 2.3 from the definitions of THM-5130
def pair_state(B, O):
    """(r, M) = (H(B) - W_B(O), U_B(O)) by brute force over ordered pairs of paths (definitions of THM-5130)."""
    n = len(B)
    h = 0
    W = 0
    U = 0
    for p in itertools.permutations(range(n)):
        arcs = [p[i + 1] in B[p[i]] for i in range(n - 1)]
        if all(arcs):
            h += 1
        # ordered pairs (P1, P2): P1 = p[:c], P2 = p[c:], both paths, P2 nonempty (c <= n - 1), P1 possibly empty
        for c in range(0, n):
            if all(arcs[:max(c - 1, 0)]) and all(arcs[c:]):
                # P1 = p[:c] is a path (arcs 0..c-2), P2 = p[c:] is a path (arcs c..n-2)
                if p[c] in O:
                    U += 1                       # P2 starts in O, P1 possibly empty
                    if c >= 1 and p[c - 1] in O and p[c - 1] != p[c]:
                        W += 1                   # both nonempty, P1 ends at e in O, P2 starts at f in O
    return h - W, U, h


def closing(B, O, X):
    """C(B, O; X): B => X => v, v -> O, B - O -> v.  Returns out-sets with vertex order B, X, v and the offsets."""
    nb, nx = len(B), len(X)
    N = nb + nx + 1
    out = [set() for _ in range(N)]
    v = nb + nx
    for a in range(nb):
        out[a] |= set(B[a])
        out[a] |= set(range(nb, nb + nx))
        if a in O:
            out[v].add(a)
        else:
            out[a].add(v)
    for a in range(nx):
        out[nb + a] |= {nb + u for u in X[a]}
        out[nb + a].add(v)
    return out


cases = 0
bad = 0
emptyO = 0
elliptic_claim_fail = 0
jobs = []
for trial in range(140):
    nb = rng.randint(1, 4)
    B0 = random_tournament(nb, rng)
    if trial < 12:
        O = set()
    else:
        O = set(v for v in range(nb) if rng.random() < 0.6)
    r, M, h = pair_state(B0, O)
    K = closing(B0, O, transitive(2))          # ports x = nb, y = nb + 1 with x -> y
    x, y = nb, nb + 1
    N = visit_matrix_dfs(K, x, y)
    Npos = {k: v for k, v in N.items() if k[0] >= 1 and k[1] >= 1}
    two = max(k[0] for k in Npos) <= 2 and max(k[1] for k in Npos) <= 2
    n11, n12, n21, n22 = (Npos.get((1, 1), 0), Npos.get((1, 2), 0), Npos.get((2, 1), 0), Npos.get((2, 2), 0))
    al = n11 - 2 * n12 - 2 * n21 + 4 * n22
    be = n21 - 2 * n22
    ga = n12 - 2 * n22
    de = n22
    det = n11 * n22 - n12 * n21
    cases += 1
    if not (two and (al, be, ga, de) == (r, 0, 0, M) and det == r * M and r % 2 == 1):
        bad += 1
        P("   FAILURE Prop 2.3: %s" % ((B0, O, r, M, two, al, be, ga, de, det),))
    if not O:
        emptyO += 1
        if not (det == 0 and M == 0 and r > 0):
            P("   unexpected for O empty")
    # the claim 'elliptic iff r > 0, hyperbolic iff r < 0'
    if (r > 0) != (det > 0) or (r < 0) != (det < 0):
        elliptic_claim_fail += 1
    na, nb2 = rng.randint(1, 4), rng.randint(1, 4)
    A = random_tournament(na, rng)
    Bm = random_tournament(nb2, rng)
    pa = P_tuple(A, 2)
    pb = P_tuple(Bm, 2)
    jobs.append((substitute(K, [x, y], [A, Bm]), r * pa[1] * pb[1] + M * (pa[2] + 2 * pa[1]) * (pb[2] + 2 * pb[1])))
Hs = ham_c([j[0] for j in jobs])
hfail = sum(1 for (T, pred), h in zip(jobs, Hs) if h != pred)
P("4. Proposition 2.3 (state from the definitions of THM-5130): pairs %d, failures of [two visits, (al,be,ga,de) = (r,0,0,M), det = rM, r odd]: %d;"
  % (cases, bad))
P("   H(K[A,B]) = r H H' + M Pi Pi' on random modules: failures %d of %d" % (hfail, len(jobs)))
P("   pairs with O EMPTY: %d (all have M = 0, det G = 0, r = H(B0) > 0: parabolic although r > 0);" % emptyO)
P("   the statement 'elliptic iff r > 0, hyperbolic iff r < 0' fails on %d of the %d pairs (exactly the pairs with O empty: %s)"
  % (elliptic_claim_fail, cases, elliptic_claim_fail == emptyO))

# ---- 5. the two named kernels
# Gaussian: B0 = point, O = {point}
r, M, h = pair_state([set()], {0})
K4 = closing([set()], {0}, transitive(2))
N = visit_matrix_dfs(K4, 1, 2)
G4 = [[N.get((1, 1), 0), N.get((1, 2), 0)], [N.get((2, 1), 0), N.get((2, 2), 0)]]
sc = sorted(len(s) for s in K4)
P("5. Gaussian kernel: state (r, M) = (%d, %d), G = %s, det %d, scores %s (strong 4-tournament has scores [1,1,2,2]), max visits %s"
  % (r, M, G4, G4[0][0] * G4[1][1] - G4[0][1] * G4[1][0], sc, max(k[0] for k in N), ))
B2 = transitive(2)
r, M, h = pair_state(B2, {0, 1})
K5 = closing(B2, {0, 1}, transitive(2))
N = visit_matrix_dfs(K5, 2, 3)
G5 = [[N.get((1, 1), 0), N.get((1, 2), 0)], [N.get((2, 1), 0), N.get((2, 2), 0)]]
P("   first hyperbolic kernel: state (r, M) = (%d, %d), G = %s, det %d, H(K) = %d"
  % (r, M, G5, G5[0][0] * G5[1][1] - G5[0][1] * G5[1][0], ham_py(K5)))
jobs = []
for trial in range(150):
    na, nb2 = rng.randint(1, 5), rng.randint(1, 5)
    A = random_tournament(na, rng)
    Bm = random_tournament(nb2, rng)
    pa = P_tuple(A, 2)
    pb = P_tuple(Bm, 2)
    piA, piB = pa[2] + 2 * pa[1], pb[2] + 2 * pb[1]
    jobs.append((substitute(K4, [1, 2], [A, Bm]), piA * piB + pa[1] * pb[1]))
    jobs.append((substitute(K5, [2, 3], [A, Bm]), 3 * piA * piB - pa[1] * pb[1]))
Hs = ham_c([j[0] for j in jobs])
P("   H = Pi Pi' + H H' (order 4) and H = 3 Pi Pi' - H H' (order 5) on 150 random module pairs each (orders up to 5): failures %d"
  % sum(1 for (T, pred), h in zip(jobs, Hs) if h != pred))
with open(os.path.join(HERE, "a1_visit_matrix.out"), "w") as fh:
    fh.write("\n".join(out_lines) + "\n")
