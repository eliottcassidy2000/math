"""AUDIT (addendum, Run I): odd-cycle parities by length, per-vertex counts, and the path form of (bK),
recomputed with the audit's own programs (aud_cyc.c: global / vertex / paths / dfs).

usage: python3 aud_cycles.py /path/aud_cyc REPO_ROOT
"""
import sys, os, json, subprocess
from math import comb
from aud_common import (paley, from_upper_code, is_doubly_regular, check_tournament, cd_tournament, double_tournament)

CYC, ROOT = sys.argv[1], sys.argv[2]
lane = os.path.join(ROOT, "04-computation/experiments/sierpinski_ocf_20261022")


def text(V):
    return str(len(V)) + "\n" + "\n".join("".join(str(x) for x in r) for r in V) + "\n"


def run(V, mode):
    out = subprocess.run([CYC, mode], input=text(V), capture_output=True, text=True, check=True).stdout.split("\n")
    return [l.split() for l in out if l.strip()]


def glob(V):
    return {int(a): int(b) for _, a, b in run(V, "global")}


def predicted_drt(n):
    m = (n - 1) // 2
    return "".join(str(comb(m, k) % 2 if n % 8 == 3 else 0) for k in range(1, m + 1))


def path_form(n):
    m = (n - 1) // 2
    res = []
    for j in range(1, n + 1):
        v = comb(n, j)
        if n % 8 == 3 and j % 2 == 1 and j >= 3:
            v += 2 * comb(m, (j - 1) // 2)
        res.append(v % 4)
    return res


print("== CD_4 (order 15, audit's own Cayley-Dickson construction): global counts, DFS enumeration, per-vertex counts")
A = cd_tournament(4)
c = glob(A)
d = {int(a): int(b) for _, a, b in run(A, "dfs")}
print("cycle counts by length (subset programme): %s" % [c[l] for l in range(3, 16)])
print("cycle counts by length (DFS enumeration):  %s ; equal: %s" % ([d[l] for l in range(3, 16)], all(c[l] == d[l] for l in range(3, 16))))
print("odd parity vector (c_3, ..., c_15): %s ; (bK) for a doubly regular tournament of order 15: %s" % (
    "".join(str(c[l] % 2) for l in range(3, 16, 2)), predicted_drt(15)))
pv = {}
for _, v, l, x in run(A, "vertex"):
    pv.setdefault(int(v), {})[int(l)] = int(x)
ok_sum = all(sum(pv[v][l] for v in range(15)) == l * c[l] for l in range(3, 16))
print("sum over vertices of (cycles of length l through v) = l c_l for every l: %s" % ok_sum)
vecs = {}
for v in range(15):
    s = "".join(str(pv[v][l] % 2) for l in range(3, 16, 2))
    vecs.setdefault(s, []).append(v)
print("parity vectors of the odd cycles through one vertex: %s" % {k: len(v) for k, v in vecs.items()})
print("vertices with an odd number of 13-cycles through them: %d of 15" % sum(1 for v in range(15) if pv[v][13] % 2))

print("== order 23: classes of drt23_classes.json (code = upper triangle, row-major) and the double of P_11")
classes = json.load(open(os.path.join(lane, "drt23_classes.json")))
todo = [(cl["label"], from_upper_code(cl["code"]), cl["H"]) for cl in classes if cl["label"] in ("D23_0", "D23_4", "D23_5", "D23_8", "P_23")]
todo.append(("double of P_11", double_tournament(paley(11)), 15759752857737405))
for name, T, H in todo:
    assert is_doubly_regular(T)
    c = glob(T)
    vec = "".join(str(c[l] % 2) for l in range(3, 24, 2))
    mod4 = "".join(str(c[l] % 4) for l in range(3, 24, 2))
    codd = sum(c[l] for l in range(3, 24, 2))
    ev = "".join(str(c[l] % 2) for l in range(4, 24, 2))
    p = {int(a): int(b) for _, a, b in run(T, "paths")}
    pmod = [p[j] % 4 for j in range(1, 24)]
    print("%-15s odd parities %s (predicted %s); c_l mod 4 (l = 3..23 odd): %s; c_odd mod 8 = %d; even parities %s; 1 + 2 c_odd = %d = H mod 4 (%d); "
          "p_23 = H mod 2^32: %s; p_j = C(23, j) (mod 4) for every j: %s" % (
              name, vec, predicted_drt(23), mod4, codd % 8, ev, (1 + 2 * codd) % 4, H % 4, p[23] == H % 2 ** 32, pmod == path_form(23)))
    sys.stdout.flush()

print("== path form of (bK): p_j mod 4 against C(n, j) + 2 [n = 3 mod 8][j odd, j >= 3] C(m, (j-1)/2)")
code19 = [l.strip() for l in open(os.path.join(lane, "drt19_codes.txt")) if l.strip() and not l.startswith("#")][0]
for name, T in (("P_3", paley(3)), ("P_7", paley(7)), ("P_11", paley(11)), ("CD_4", cd_tournament(4)), ("P_19", paley(19)),
                ("D_19", from_upper_code(code19))):
    n = len(T)
    p = {int(a): int(b) for _, a, b in run(T, "paths")}
    pmod = [p[j] % 4 for j in range(1, n + 1)]
    c = glob(T)
    vec = "".join(str(c[l] % 2) for l in range(3, n + 1, 2))
    # the subset identity (p_j - C(n,j))/2 = sum_l C(n-l, j-l) c_l (mod 2), l odd  [exact counts are below 2^32 only for n <= 15]
    print("%-5s n = %2d: p_j mod 4 = %s ; formula: %s ; equal: %s ; odd cycle parities %s (predicted %s)" % (
        name, n, "".join(map(str, pmod)), "".join(map(str, path_form(n))), pmod == path_form(n), vec, predicted_drt(n)))

print("== the subset identity (p_j - C(n,j))/2 = sum over odd l of C(n-l, j-l) c_l (mod 2), exact counts, random tournaments")
import random
random.seed(5132)
bad = tot = 0
for n in range(3, 13):
    for _ in range(40):
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        c = glob(A)
        p = {int(a): int(b) for _, a, b in run(A, "paths")}
        for j in range(1, n + 1):
            lhs = (p[j] - comb(n, j))
            assert lhs % 2 == 0
            rhs = sum(comb(n - l, j - l) * c[l] for l in range(3, j + 1, 2))
            tot += 1
            if (lhs // 2 - rhs) % 2:
                bad += 1
print("identities tested: %d, violations: %d" % (tot, bad))
