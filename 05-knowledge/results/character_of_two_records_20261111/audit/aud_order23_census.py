"""AUDIT: how many of the 37 doubly regular tournaments of order 23 (count from the literature: Spence; McKay's
catalogue) are reached from the repository's data, and do (b4), (b8) and the normalised form of (bK) hold on all of them?

Sources: the 12 annealed classes + P_23 (drt23_classes.json), the double of P_11, all their converses, and the 24
descendants (normalise a vertex of (source) + T to a source, delete it) of each of these skew Hadamard classes.
Isomorphism classes by nauty canonical forms (dreadnaut). For one representative of every class found:
H mod 2^32 (aud_ham.c) and the parities of the odd cycle counts (aud_cyc.c).
usage: python3 aud_order23_census.py /path/aud_ham /path/aud_cyc REPO_ROOT
"""
import sys, os, json, subprocess, hashlib
from aud_common import paley, from_upper_code, is_doubly_regular, double_tournament

HAM, CYC, ROOT = sys.argv[1:4]
lane = os.path.join(ROOT, "04-computation/experiments/sierpinski_ocf_20261022")


def canon(A):
    n = len(A)
    lines = ["n=%d d" % n, "g"]
    for i in range(n):
        lines.append("%d : %s;" % (i, " ".join(str(j) for j in range(n) if A[i][j])))
    lines += ["c x", "b", "q"]
    out = subprocess.run(["dreadnaut"], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout
    idx = out.rfind("grpsize")
    tail = out[idx:]
    grp = tail.split(";")[0].replace("grpsize=", "").strip()
    adj = "\n".join(l.strip() for l in tail.split("\n") if ":" in l and "grpsize" not in l)
    return hashlib.sha256(adj.encode()).hexdigest()[:20], grp


def conv(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def descendants(T):
    n = len(T)
    N = n + 1
    U = [[0] * N for _ in range(N)]
    for j in range(1, N):
        U[0][j] = 1
    for a in range(n):
        for b in range(n):
            U[a + 1][b + 1] = T[a][b]
    res = []
    for w in range(N):
        inX = [x != w and U[x][w] == 1 for x in range(N)]
        V = [[(1 - U[a][b]) if (a != b and inX[a] != inX[b]) else U[a][b] for b in range(N)] for a in range(N)]
        idx = [y for y in range(N) if y != w]
        res.append([[V[a][b] for b in idx] for a in idx])
    return res


def text(V):
    return str(len(V)) + "\n" + "\n".join("".join(str(x) for x in r) for r in V) + "\n"


classes23 = json.load(open(os.path.join(lane, "drt23_classes.json")))
bases = [(cl["label"], from_upper_code(cl["code"])) for cl in classes23]
bases.append(("double of P_11", double_tournament(paley(11))))
found = {}      # canonical hash -> (representative, origin, |Aut|)
recorded = set()
for name, T in bases:
    for tag, X in ((name, T), (name + " (converse)", conv(T))):
        h, g = canon(X)
        recorded.add(h)
        found.setdefault(h, (X, tag, g))
n_rec = len(recorded)
for name, T in bases:
    for i, D in enumerate(descendants(T)):
        assert is_doubly_regular(D)
        for tag, X in (("descendant %d of the class of %s" % (i, name), D), ("converse of descendant %d of the class of %s" % (i, name), conv(D))):
            h, g = canon(X)
            found.setdefault(h, (X, tag, g))
print("isomorphism classes among the recorded tournaments and their converses: %d" % n_rec)
print("isomorphism classes after adding all descendants and their converses: %d (the literature counts 37 doubly regular tournaments of order 23)" % len(found))
# up to converse
pairs = {}
for h, (X, tag, g) in found.items():
    hc, _ = canon(conv(X))
    key = tuple(sorted((h, hc)))
    pairs.setdefault(key, (X, tag, g, h == hc))
print("classes up to converse: %d (self-converse: %d)" % (len(pairs), sum(1 for v in pairs.values() if v[3])))
res8 = {}
vecs = {}
auts = {}
for key, (X, tag, g, sc) in sorted(pairs.items(), key=lambda kv: kv[1][1]):
    hh = int(subprocess.run([HAM, "0"], input=text(X), capture_output=True, text=True, check=True).stdout)
    out = subprocess.run([CYC, "global"], input=text(X), capture_output=True, text=True, check=True).stdout.split("\n")
    c = {int(l.split()[1]): int(l.split()[2]) for l in out if l.strip()}
    v = "".join(str(c[l] % 2) for l in range(3, 24, 2))
    res8[hh % 8] = res8.get(hh % 8, 0) + 1
    vecs[v] = vecs.get(v, 0) + 1
    auts[g] = auts.get(g, 0) + 1
    print("  %-52s |Aut| = %-4s self-converse: %-5s H mod 2^32 = %10d  H mod 8 = %d  odd cycle parities %s  c_9 = %d" % (
        tag[:52], g, sc, hh, hh % 8, v, c[9]))
    sys.stdout.flush()
print("H mod 8 over the classes (up to converse): %s; odd parity vectors: %s; |Aut| distribution: %s" % (res8, vecs, auts))
