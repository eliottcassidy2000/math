"""aud_a_gadget_law.py - audit tasks 2(a), 2(b).

T0  my subset recursion agrees with plain permutation enumeration (all labelled tournaments of order <= 4,
    every O; random pairs of orders 5..8).
T1  gadget law (Theorem A): for random (B, O, X, P) with |B|, |X| <= 4 the state of (C(B,O;X), P), counted
    literally, equals G(X,P) s(B,O) with G built from gadget constants counted by permutation enumeration.
    The literal state of B+ is counted by permutation enumeration when |B+| <= 8 and by the subset recursion
    otherwise (and by both for a subsample at order 9).
T1e edge cases: O empty, O = V(B), P empty, P = V(X), |B| = 1, |X| = 1, exhaustive over all labelled B, X of
    order <= 3 and all O, P (subset recursion, cross-checked by enumeration at order <= 7).
T2  Theorem B: constants of (TT_a, {i, i+delta}) for delta = 1..5, i = 1..4, m = 0..3 (a <= 9 by enumeration
    up to a = 7, by set-partition counting beyond), including i = 1, m = 0.
T3  Remarks (1)-(3): empty base convention, sink / source rules, closing law.
"""
import itertools, random, sys, time
from aud_core import *

rng = random.Random(987654321)
allok = True


def report(name, ok, extra=""):
    global allok
    allok &= bool(ok)
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + extra if extra else ""))
    sys.stdout.flush()


def rand_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def rand_subset(n, p=0.5):
    return [i for i in range(n) if rng.random() < p]


# ------------------------------------------------------------------ T0
t0 = time.time()
bad = cnt = 0
for m in range(1, 5):
    for code in range(1 << (m * (m - 1) // 2)):
        A = adj_from_code(m, code)
        for Om in range(1 << m):
            O = [i for i in range(m) if Om >> i & 1]
            cnt += 1
            bad += (brute_state(A, O) != dp_state(A, O))
for n in (5, 6, 7, 8):
    for _ in range(40 if n < 8 else 12):
        A = rand_tournament(n)
        O = rand_subset(n)
        cnt += 1
        bad += (brute_state(A, O) != dp_state(A, O))
report("T0 subset recursion = permutation enumeration", bad == 0, "%d pairs, %.0fs" % (cnt, time.time() - t0))

# ------------------------------------------------------------------ T1
t0 = time.time()
bad = cnt = cnt_brute = cnt_both9 = 0
offsets = set()
fails = []
for trial in range(400):
    n = rng.randint(1, 4)
    k = rng.randint(1, 4)
    A = rand_tournament(n)
    O = rand_subset(n, 0.6)
    X = rand_tournament(k)
    P = rand_subset(k, 0.6)
    s = brute_state(A, O)
    offsets.add(s[0])
    G = gadget_matrix(brute_gadget(X, P))
    pred = mat_apply(G, s)
    Bp = close(A, O, X)
    Pp = [n + p for p in P]
    N = n + k + 1
    if N <= 8:
        lit = brute_state(Bp, Pp)
        cnt_brute += 1
    else:
        lit = dp_state(Bp, Pp)
        if cnt_both9 < 12:
            cnt_both9 += 1
            if brute_state(Bp, Pp) != lit:
                bad += 1
                fails.append(("brute/dp disagree", A, O, X, P))
    cnt += 1
    if lit != pred:
        bad += 1
        fails.append((A, O, X, P, s, pred, lit))
report("T1 gadget law on random (B, O, X, P), orders <= 4 each", bad == 0,
       "%d cases (%d literal by permutations, %d at order 9 by both), offsets seen: %s, %.0fs"
       % (cnt, cnt_brute, cnt_both9, sorted(offsets), time.time() - t0))
for f in fails[:5]:
    print("   FAIL", f)

# ------------------------------------------------------------------ T1e exhaustive small cases incl. edge cases
t0 = time.time()
bad = cnt = 0
gcache = {}
for n in range(1, 4):
    for codeB in range(1 << (n * (n - 1) // 2)):
        A = adj_from_code(n, codeB)
        for Om in range(1 << n):
            O = [i for i in range(n) if Om >> i & 1]
            s = brute_state(A, O)
            for k in range(1, 4):
                for codeX in range(1 << (k * (k - 1) // 2)):
                    X = adj_from_code(k, codeX)
                    for Pm in range(1 << k):
                        P = [i for i in range(k) if Pm >> i & 1]
                        key = (k, codeX, Pm)
                        if key not in gcache:
                            gcache[key] = gadget_matrix(brute_gadget(X, P))
                        pred = mat_apply(gcache[key], s)
                        Bp = close(A, O, X)
                        lit = dp_state(Bp, [n + p for p in P])
                        cnt += 1
                        if lit != pred:
                            bad += 1
                            if bad < 5:
                                print("   FAIL", n, codeB, O, k, codeX, P, s, pred, lit)
report("T1e gadget law, exhaustive: all labelled B, X of order <= 3, every O and P (empty and full included)",
       bad == 0, "%d cases, %.0fs" % (cnt, time.time() - t0))


# ------------------------------------------------------------------ T2 transitive pair constants
def tt_constants_setpartition(a, P):
    """constants of (TT_a, P) by literal enumeration of ordered set partitions (blocks = increasing paths).
    vertices 1..a; P a set of 1-based vertices."""
    Ps = set(P)
    V = list(range(1, a + 1))
    HX = 1
    Pi = w = w3 = u = nu = u3 = 0
    sigma = 1 if 1 in Ps else 0
    eps = 1 if a in Ps else 0
    # ordered pairs (D1, D2), possibly empty
    for lab in itertools.product((0, 1), repeat=a):
        D1 = [v for v, l in zip(V, lab) if l == 0]
        D2 = [v for v, l in zip(V, lab) if l == 1]
        Pi += 1
        if D2 and D2[0] in Ps:
            u += 1
        if D1 and D1[-1] in Ps:
            nu += 1
        if D1 and D2 and D1[-1] in Ps and D2[0] in Ps:
            w += 1
    for lab in itertools.product((0, 1, 2), repeat=a):
        D = [[v for v, l in zip(V, lab) if l == c] for c in range(3)]
        if not (D[0] and D[1] and D[2]):
            continue
        if D[0][-1] in Ps and D[1][0] in Ps:
            w3 += 1
        if D[0][0] in Ps:
            u3 += 1
    return dict(H=HX, Pi=Pi, w=w, w3=w3, sigma=sigma, u=u, eps=eps, nu=nu, u3=u3)


t0 = time.time()
ok = True
cnt = 0
for delta in (1, 2, 3, 4, 5):
    for i in (1, 2, 3, 4):
        for m in (0, 1, 2, 3):
            a = i + delta + m
            if a > 11:
                continue
            P1 = [i, i + delta]                     # 1-based
            g = tt_constants_setpartition(a, P1)
            if a <= 7:
                gb = brute_gadget(transitive(a), [i - 1, i - 1 + delta])
                ok &= (gb == g)
            rX = -1 if delta == 1 else 1 - 2 ** (delta - 1)
            kap = 2 ** (a - delta) * (2 ** delta - 1 - 3 ** (delta - 1))
            exp = dict(rX=rX, kap=kap, u=(2 ** delta + 1) * 2 ** m,
                       mid=2 * (3 ** delta + 2 ** delta) * 2 ** (i - 1) * 3 ** m, sigma=1 if i == 1 else 0,
                       eps=1 if m == 0 else 0, nu=(2 ** delta + 1) * 2 ** (i - 1), H=1, Pi=2 ** a)
            got = dict(rX=g["H"] - g["w"], kap=g["Pi"] - 2 * g["w"] - 2 * g["w3"], u=g["u"],
                       mid=4 * g["u"] - 2 * g["sigma"] + 2 * g["u3"], sigma=g["sigma"], eps=g["eps"], nu=g["nu"],
                       H=g["H"], Pi=g["Pi"])
            # intermediate closed forms of the proof
            wexp = 2 if delta == 1 else 2 ** (delta - 1)
            w3exp = 2 ** (a - delta - 1) * (1 + 3 ** (delta - 1)) - (1 if delta == 1 else 0) - 2 ** (delta - 1)
            u3exp = sum(2 ** (j - 1) * 3 ** (a - j) - 2 ** (a - j + 1) + (1 if j == 1 else 0) for j in P1)
            ok &= (exp == got and g["w"] == wexp and g["w3"] == w3exp and g["u3"] == u3exp)
            cnt += 1
            if exp != got:
                print("   FAIL", delta, i, m, exp, got)
report("T2 Theorem B constants (delta <= 5, i <= 4, m <= 3, a <= 11), also w, w3, u3 of the proof", ok,
       "%d gadgets, %.0fs" % (cnt, time.time() - t0))

# sign-flip criterion on literal states: delta = 1, 2 reverse the offset, delta = 3, 4 do not
ok = True
cnt = 0
for _ in range(60):
    n = rng.randint(1, 4)
    A = rand_tournament(n)
    O = rand_subset(n, 0.6)
    s = brute_state(A, O)
    for delta in (1, 2, 3, 4):
        i = rng.randint(1, 2)
        m = rng.randint(0, 1)
        a = i + delta + m
        if n + a + 1 > 11:
            continue
        lit = dp_state(close(A, O, transitive(a)), [n + i - 1, n + i - 1 + delta])
        rX = -1 if delta <= 2 else 1 - 2 ** (delta - 1)
        kap = 2 ** (a - delta) * (2 ** delta - 1 - 3 ** (delta - 1))
        alpha, gamma = 2 ** delta + 1, 3 ** delta + 2 ** delta
        pred = (rX * s[0] + kap * s[1],
                alpha * 2 ** m * (s[2] + s[0]) + 2 * gamma * 2 ** (i - 1) * 3 ** m * s[1] - (s[3] if i == 1 else 0),
                s[0] + 2 ** a * s[1],
                alpha * 2 ** (i - 1) * s[1] + ((s[0] - s[2]) if m == 0 else 0))
        cnt += 1
        ok &= (lit == pred)
        if delta <= 2:
            ok &= (lit[0] == -s[0])
report("T2b level map (1) and Corollary B on literal states (delta = 1..4, random bases)", ok, "%d cases" % cnt)

# ------------------------------------------------------------------ T3 remarks
ok = True
for k in range(1, 5):
    for codeX in range(1 << (k * (k - 1) // 2)):
        X = adj_from_code(k, codeX)
        for Pm in range(1 << k):
            P = [i for i in range(k) if Pm >> i & 1]
            G = gadget_matrix(brute_gadget(X, P))
            # X => pt
            T = [row + [1] for row in X] + [[0] * (k + 1)]
            ok &= (brute_state(T, P) == mat_apply(G, (1, 0, 1, 0)))
report("T3a Remark (1): G(X,P)(1,0,1,0) = state of (X => pt, P), all labelled X of order <= 4, every P", ok)
ok = True
for _ in range(80):
    n = rng.randint(1, 6)
    A = rand_tournament(n)
    O = rand_subset(n)
    r, M, h, E = brute_state(A, O)
    S = [row + [1] for row in A] + [[0] * (n + 1)]
    ok &= (brute_state(S, O) == (r, 2 * M, h, 0))
    R = [[0] + [1] * n] + [[0] + row for row in A]
    ok &= (brute_state(R, [o + 1 for o in O]) == (r, M, h, E))
report("T3b Remark (2): sink -> (r, 2M, h, 0); source -> unchanged", ok)
ok = True
cnt = 0
for _ in range(80):
    n = rng.randint(1, 4)
    A = rand_tournament(n)
    O = rand_subset(n)
    r, M, h, E = brute_state(A, O)
    for t in (1, 2, 3):
        cnt += 1
        ok &= (brute_H(close(A, O, transitive(t))) == r + 2 ** t * M)
report("T3c closing law H(T_t(B,O)) = r + 2^t M by permutation enumeration (|B| <= 4, t <= 3)", ok, "%d cases" % cnt)
# closing law with arbitrary X: H(C(B,O;X)) = H(X) r + Pi_2(X) M
ok = True
for _ in range(80):
    n = rng.randint(1, 4)
    k = rng.randint(1, 3)
    A = rand_tournament(n)
    O = rand_subset(n)
    X = rand_tournament(k)
    r, M, h, E = brute_state(A, O)
    g = brute_gadget(X, [])
    ok &= (brute_H(close(A, O, X)) == g["H"] * r + g["Pi"] * M)
report("T3d row h: H(C(B,O;X)) = H(X) r + Pi_2(X) M by permutation enumeration", ok)
print("ALL PASS" if allok else "SOME FAILED")
