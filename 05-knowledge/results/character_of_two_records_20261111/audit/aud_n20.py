"""AUDIT: (bK) on the two skew Hadamard switching classes of order 20 (classes of (source) + P_19 and (source) + D_19;
by the classification there are exactly two doubly regular tournaments of order 19) and on the single-vertex
switchings of the structured classes of order 24.
usage: python3 aud_n20.py /path/aud_cyc REPO_ROOT"""
import sys, os, random, subprocess
from math import comb
from aud_common import paley, from_upper_code, is_doubly_regular, double_tournament

CYC, ROOT = sys.argv[1], sys.argv[2]
lane = os.path.join(ROOT, "04-computation/experiments/sierpinski_ocf_20261022")
code19 = [l.strip() for l in open(os.path.join(lane, "drt19_codes.txt")) if l.strip() and not l.startswith("#")][0]


def with_source(T):
    n = len(T)
    U = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(1, n + 1):
        U[0][j] = 1
    for a in range(n):
        for b in range(n):
            U[a + 1][b + 1] = T[a][b]
    return U


def switch(U, X):
    N = len(U)
    inX = [False] * N
    for x in X:
        inX[x] = True
    return [[(1 - U[a][b]) if (a != b and inX[a] != inX[b]) else U[a][b] for b in range(N)] for a in range(N)]


def vec(V):
    N = len(V)
    text = str(N) + "\n" + "\n".join("".join(str(x) for x in r) for r in V) + "\n"
    out = subprocess.run([CYC, "global"], input=text, capture_output=True, text=True, check=True).stdout.split("\n")
    c = {int(l.split()[1]): int(l.split()[2]) for l in out if l.strip()}
    return "".join(str(c[l] % 2) for l in range(3, N, 2))


def pred(N):
    m = N // 2 - 1
    return "".join(str(comb(m, k) % 2 if N % 8 == 4 else 0) for k in range(1, m + 1))


random.seed(20)
for name, T in (("P_19", paley(19)), ("D_19", from_upper_code(code19))):
    U = with_source(T)
    N = len(U)
    res = {}
    for w in range(1, N):
        v = vec(switch(U, [w]))
        res[v] = res.get(v, 0) + 1
    res2 = {}
    for _ in range(150):
        X = [v for v in range(1, N) if random.random() < 0.5]
        v = vec(switch(U, X))
        res2[v] = res2.get(v, 0) + 1
    print("N = 20, class of (source) + %s: predicted %s; 19 single-vertex switchings: %s; 150 random members: %s" % (name, pred(N), res, res2))
    sys.stdout.flush()
for name, T, ws in (("P_23", paley(23), [1]), ("double of P_11", double_tournament(paley(11)), list(range(1, 24)))):
    U = with_source(T)
    N = len(U)
    res = {}
    for w in ws:
        v = vec(switch(U, [w]))
        res[v] = res.get(v, 0) + 1
    print("N = 24, class of (source) + %s: predicted %s; single-vertex switchings (%d tested%s): %s" % (
        name, pred(N), len(ws), "; all 23 are isomorphic, Aut is vertex-transitive" if name == "P_23" else "", res))
    sys.stdout.flush()
