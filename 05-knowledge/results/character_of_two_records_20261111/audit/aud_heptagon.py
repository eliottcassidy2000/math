"""AUDIT (second addendum): Proposition H of section 5.4 (7-cycles of doubly regular tournaments) and the sentence
"Short cycles are determined by the order" of section 5.3.

usage: python3 aud_heptagon.py /path/aud_cyc REPO_ROOT
 1. symbolic: tr(A^7)/7 - sum_u t_u f_u + N_2 + 2 N_3 with the values of the note equals both closed forms;
 2. the general identities on random tournaments, by brute-force enumeration of closed 7-walks:
    tr(A^7) = 7 (c_7 + N_1 + N_2 + N_3), sum_u t_u f_u = N_1 + 2 N_2 + 3 N_3, at most three repeated vertices,
    every repeated vertex visited exactly twice, N_3 = #{3-cycle, arc z -> u of it, y with z -> y -> u},
    N_2 = (type 1,2,1,3) + (type 1,2,2,2) with the descriptions of the note;
 3. exact c_7 and c_8 (audit's cycle counter) against the formulas, orders 7..23;
 4. c_3..c_9 equal across classes of one order (19: P_19, D_19; 23: all classes of drt23_classes.json + the double).
"""
import sys, os, json, random, itertools, subprocess
import sympy as sp
from aud_common import paley, from_upper_code, is_doubly_regular, cd_tournament, double_tournament

CYC, ROOT = sys.argv[1], sys.argv[2]
lane = os.path.join(ROOT, "04-computation/experiments/sierpinski_ocf_20261022")

print("== 1. symbolic check of the closed form")
t = sp.symbols("t")
n = 4 * t + 3
q = t + 1
tr7 = (2 * t + 1) * ((2 * t + 1) ** 6 + 7 * q ** 3 - 14 * q ** 2 + 7 * q - 1)
sum_tf = 2 * n * t * (t + 1) ** 2 * (2 * t + 1) ** 2
N2 = sp.Rational(3, 2) * n * (n - 1) * t ** 2 * (t + 1)
c3 = n * (n ** 2 - 1) / 24
N3 = 3 * t * c3
c7 = sp.expand(tr7 / 7 - sum_tf + N2 + 2 * N3)
form_n = sp.expand(n * (n + 1) * (n - 1) * (n - 3) * (n - 5) * (n ** 2 - 6 * n + 1) / 896)
form_t = sp.expand(4 * t * (t + 1) * (2 * t - 1) * (2 * t + 1) * (4 * t + 3) * (2 * t ** 2 - 1) / 7)
print("derivation = n-form: %s ; derivation = t-form: %s" % (sp.simplify(c7 - form_n) == 0, sp.simplify(c7 - form_t) == 0))
vals = [int(form_t.subs(t, k)) for k in range(0, 9)]
print("c_7 for n = 3, 7, ..., 35: %s ; all divisible by 8: %s ; integral for t <= 200: %s" % (
    vals, all(v % 8 == 0 for v in vals), all(sp.Integer(4 * k * (k + 1) * (2 * k - 1) * (2 * k + 1) * (4 * k + 3) * (2 * k * k - 1)) % 7 == 0 for k in range(201))))
# power sums of the roots of x^2 + x + q
s = [sp.Integer(2), sp.Integer(-1)]
qq = sp.symbols("q")
for k in range(2, 8):
    s.append(sp.expand(-s[k - 1] - qq * s[k - 2]))
print("s_7 = %s" % s[7])

print("== 2. general identities, brute force on random tournaments")
random.seed(7)


def rand_t(nv):
    A = [[0] * nv for _ in range(nv)]
    for i in range(nv):
        for j in range(i + 1, nv):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def closed_walks7(A):
    nv = len(A)
    out = [[j for j in range(nv) if A[i][j]] for i in range(nv)]
    res = []
    def rec(w):
        if len(w) == 7:
            if A[w[-1]][w[0]]:
                res.append(tuple(w))
            return
        for x in out[w[-1]]:
            w.append(x)
            rec(w)
            w.pop()
    for v in range(nv):
        rec([v])
    return res


def check_general(A):
    nv = len(A)
    walks = closed_walks7(A)
    # rotation classes
    classes = set()
    for w in walks:
        classes.add(min(tuple(w[i:] + w[:i]) for i in range(7)))
    assert len(walks) == 7 * len(classes)
    N = [0, 0, 0, 0, 0]
    maxmult = 0
    types = {}
    for w in classes:
        cnt = {}
        for x in w:
            cnt[x] = cnt.get(x, 0) + 1
        rep = [x for x, c in cnt.items() if c >= 2]
        maxmult = max(maxmult, max(cnt.values()))
        N[len(rep)] += 1
        if len(rep) == 2:
            # segment lengths between consecutive repeated visits, as a cyclic sequence
            pos = [i for i in range(7) if cnt[w[i]] == 2]
            segs = [(pos[(k + 1) % 4] - pos[k]) % 7 for k in range(4)]
            key = min(tuple(segs[i:] + segs[:i]) for i in range(4))
            types[key] = types.get(key, 0) + 1
    # cycles through vertices
    def cyc_through(u, L):
        c = 0
        def rec(w):
            nonlocal c
            if len(w) == L:
                if A[w[-1]][u]:
                    c += 1
                return
            for x in range(nv):
                if A[w[-1]][x] and x not in w:
                    w.append(x)
                    rec(w)
                    w.pop()
        rec([u])
        return c
    tu = [cyc_through(u, 3) for u in range(nv)]
    fu = [cyc_through(u, 4) for u in range(nv)]
    c3cycles = [(a, b, c) for a in range(nv) for b in range(nv) for c in range(nv) if a < b and a < c and A[a][b] and A[b][c] and A[c][a]]
    N3_desc = 0
    for (a, b, c) in c3cycles:
        for (z, u) in ((a, b), (b, c), (c, a)):
            N3_desc += sum(1 for y in range(nv) if A[z][y] and A[y][u])
    # N_2 descriptions
    T1 = T2 = 0
    for u in range(nv):
        for w_ in range(nv):
            if not A[u][w_]:
                continue
            back = [x for x in range(nv) if A[w_][x] and A[x][u]]          # w -> x -> u
            fwd = [x for x in range(nv) if A[u][x] and A[x][w_]]            # u -> x -> w
            T2 += len(fwd) * len(back) * (len(back) - 1)
            for a in back:
                for b in range(nv):
                    if not A[w_][b] or b in (u, a):
                        continue
                    for c in range(nv):
                        if c in (u, w_, a, b):
                            continue
                        if A[b][c] and A[c][u]:
                            T1 += 1
    c7 = N[0]
    ok = (len(walks) == 7 * (c7 + N[1] + N[2] + N[3]) and N[4] == 0 and maxmult <= 2
          and sum(a * b for a, b in zip(tu, fu)) == N[1] + 2 * N[2] + 3 * N[3]
          and N3_desc == N[3] and set(types) <= {(1, 2, 1, 3), (1, 2, 2, 2)}
          and types.get((1, 2, 1, 3), 0) == T1 and types.get((1, 2, 2, 2), 0) == T2)
    return ok, c7, N, types


allok = True
for nv in (7, 8, 9, 10):
    for _ in range(6):
        ok, c7, N, types = check_general(rand_t(nv))
        allok &= ok
ok, c7, N, types = check_general(paley(7))
print("random tournaments of orders 7..10 (24 cases): all identities and the classification hold: %s" % allok)
print("P_7: c_7 = %d, (N_1, N_2, N_3) = %s, types of N_2: %s, identities hold: %s; formulas: N_2 = %d, N_3 = %d" % (
    c7, N[1:4], types, ok, 3 * 7 * 6 * 1 * 2 // 2, 3 * 1 * 14))
ok, c7, N, types = check_general(paley(11))
print("P_11: c_7 = %d, (N_1, N_2, N_3) = %s, types of N_2: %s, identities hold: %s; formulas: N_2 = %d, N_3 = %d" % (
    c7, N[1:4], types, ok, 3 * 11 * 10 * 4 * 3 // 2, 3 * 2 * 55))

print("== 3. exact c_7, c_8 against the formulas")


def text(V):
    return str(len(V)) + "\n" + "\n".join("".join(str(x) for x in r) for r in V) + "\n"


def glob(V):
    out = subprocess.run([CYC, "global"], input=text(V), capture_output=True, text=True, check=True).stdout.split("\n")
    return {int(l.split()[1]): int(l.split()[2]) for l in out if l.strip()}


def f7(nn):
    return nn * (nn + 1) * (nn - 1) * (nn - 3) * (nn - 5) * (nn * nn - 6 * nn + 1) // 896


def f8(nn):
    return (nn + 1) * nn * (nn - 1) * (nn - 3) ** 2 * (nn - 7) * (nn * nn - 7 * nn + 1) // 2048


code19 = [l.strip() for l in open(os.path.join(lane, "drt19_codes.txt")) if l.strip() and not l.startswith("#")][0]
classes23 = json.load(open(os.path.join(lane, "drt23_classes.json")))
cases = [("P_7", paley(7)), ("P_11", paley(11)), ("CD_4", cd_tournament(4)), ("P_19", paley(19)), ("D_19", from_upper_code(code19))]
cases += [(cl["label"], from_upper_code(cl["code"])) for cl in classes23]
cases.append(("double of P_11", double_tournament(paley(11))))
short = {}
for name, T in cases:
    nn = len(T)
    assert is_doubly_regular(T)
    c = glob(T)
    short.setdefault(nn, []).append((name, [c[l] for l in range(3, min(nn, 12) + 1)]))
    print("%-15s n = %2d: c_7 = %d (formula %d: %s); c_8 = %d (Savchenko's formula %d: %s)" % (
        name, nn, c[7], f7(nn), c[7] == f7(nn) % 2 ** 32, c.get(8, 0), f8(nn), c.get(8, 0) == f8(nn) % 2 ** 32))

print("== 4. short cycles across classes of one order (counts mod 2^32; exact below 2^32)")
for nn in (19, 23):
    rows = short[nn]
    upto = {}
    for L in range(3, 13):
        vals = set(r[1][L - 3] for r in rows)
        upto[L] = len(vals)
    print("order %d: %d classes; number of distinct values of c_L for L = 3..12: %s" % (nn, len(rows), [upto[L] for L in range(3, 13)]))
    print("   common values c_3..c_9: %s" % rows[0][1][:7])
