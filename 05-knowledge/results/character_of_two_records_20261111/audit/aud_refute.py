"""AUDIT: how far do (bK) and (bS) of HYP-9667 fail at N = 24, and do (b4), (b8) survive on the descendants?

usage: python3 aud_refute.py HAM CYC REPO_ROOT MODE [args]
  MODE = sample K            K random members (seeded) of the switching class of (source) + T for every doubly regular
                             T of order 23 on record (annealed codes 0..9 of drt23_codes_all.txt, P_23, double of P_11):
                             parity vector of the odd cycle counts and H mod 4 (= 1 + 2 c_odd), against (bK), (bS).
  MODE = single CODE         the 23 single-vertex switchings of (source) + T, T = annealed code CODE.
  MODE = desc CODE [CODE..]  the 24 descendants (normalise a vertex to a source and delete it) of the class of
                             (source) + T: each is doubly regular of order 23; H mod 8 and odd cycle parities, against
                             (b4), (b8) and the normalised form of (bK). CODE = 0..9, or P for P_23, or D for the double.
"""
import sys, os, random, subprocess
from aud_common import (paley, from_upper_code, is_doubly_regular, check_tournament, double_tournament)

HAM, CYC, ROOT, MODE = sys.argv[1:5]
lane = os.path.join(ROOT, "04-computation/experiments/sierpinski_ocf_20261022")
codes23 = [l.strip() for l in open(os.path.join(lane, "drt23_codes_all.txt")) if l.strip() and not l.startswith("#")]


def base(name):
    if name == "P":
        return paley(23)
    if name == "D":
        return double_tournament(paley(11))
    return from_upper_code(codes23[int(name)])


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


def text(V):
    return str(len(V)) + "\n" + "\n".join("".join(str(x) for x in r) for r in V) + "\n"


def cyc(V):
    out = subprocess.run([CYC, "global"], input=text(V), capture_output=True, text=True, check=True).stdout.split("\n")
    c = {}
    for l in out:
        if l.strip():
            _, a, b = l.split()
            c[int(a)] = int(b)
    return c


def ham32(V):
    return int(subprocess.run([HAM, "0"], input=text(V), capture_output=True, text=True, check=True).stdout)


def vec(c, N):
    return "".join(str(c[l] % 2) for l in range(3, N, 2)) if N % 2 == 0 else "".join(str(c[l] % 2) for l in range(3, N + 1, 2))


if MODE == "sample":
    K = int(sys.argv[5])
    random.seed(20261111)
    tot = badK = badS = 0
    for name in [str(i) for i in range(10)] + ["P", "D"]:
        T = base(name)
        assert is_doubly_regular(T)
        U = with_source(T)
        N = len(U)
        vs = {}
        h4 = {}
        for _ in range(K):
            X = [v for v in range(1, N) if random.random() < 0.5]
            c = cyc(switch(U, X))
            v = vec(c, N)
            codd = sum(c[l] for l in range(3, N, 2))
            h = (1 + 2 * codd) % 4
            vs[v] = vs.get(v, 0) + 1
            h4[h] = h4.get(h, 0) + 1
            tot += 1
            if v != "0" * 11:
                badK += 1
            if h != 1:
                badS += 1
        print("class of (source) + %s: %d random members; odd parity vectors (c_3..c_23): %s; H mod 4: %s" % (
            {"P": "P_23", "D": "double of P_11"}.get(name, "annealed code " + name), K, dict(sorted(vs.items())), dict(sorted(h4.items()))))
        sys.stdout.flush()
    print("total members %d: violating (bK): %d; violating (bS): %d" % (tot, badK, badS))

elif MODE == "single":
    T = base(sys.argv[5])
    U = with_source(T)
    N = len(U)
    c0 = cyc(U)
    print("normalised member: odd parity vector %s" % vec(c0, N))
    bad = 0
    for w in range(1, N):
        c = cyc(switch(U, [w]))
        v = vec(c, N)
        codd = sum(c[l] for l in range(3, N, 2))
        if v != "0" * 11:
            bad += 1
        print("switched at the single vertex %2d: odd parity vector %s; H mod 4 = %d" % (w, v, (1 + 2 * codd) % 4))
        sys.stdout.flush()
    print("single-vertex switchings violating (bK): %d of %d" % (bad, N - 1))

elif MODE == "desc":
    for name in sys.argv[5:]:
        T = base(name)
        U = with_source(T)
        N = len(U)
        res8 = {}
        vecs = {}
        notdr = 0
        hs = []
        for w in range(N):
            X = [x for x in range(N) if x != w and U[x][w]]          # in-neighbours of w: after switching w is a source
            V = switch(U, X)
            assert all(V[w][y] for y in range(N) if y != w)
            idx = [y for y in range(N) if y != w]
            D = [[V[a][b] for b in idx] for a in idx]
            if not is_doubly_regular(D):
                notdr += 1
            h = ham32(D)
            c = cyc(D)
            v = vec(c, len(D))
            res8[h % 8] = res8.get(h % 8, 0) + 1
            vecs[v] = vecs.get(v, 0) + 1
            hs.append(h)
        print("descendants of the class of (source) + %s: 24 tournaments, not doubly regular: %d; distinct H mod 2^32: %d; H mod 8: %s; odd parity vectors: %s" % (
            {"P": "P_23", "D": "double of P_11"}.get(name, "annealed code " + name), notdr, len(set(hs)), dict(sorted(res8.items())), dict(sorted(vecs.items()))))
        sys.stdout.flush()
