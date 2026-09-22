#!/usr/bin/env python3
"""Adversarial audit of lane lean_paste_audit_w6 (2026-09-22).

Independent recomputation, by different routes than the lane script:
  A. self-loops 2m = k^2 (closed form: m = 2 j^2) and the pasted Fin-n adjacency;
  B. components of Q_n, n = 0..32, by union-find (lane used DFS);
  C. degree <= 1 vertices by a numpy adjacency matrix (lane used adjacency lists);
  D. Hamiltonian paths: bitmask DP (Held-Karp reachability) for n <= 20, and a
     second backtracker with different pruning/order for 21 <= n <= 32; witness
     validation of the two session-lead paths; the n = 15 '4' claim;
  E. the Delta = 4 ladder: the mod-3 proof that (3,7,11) is the ONLY prime
     triple (p, p+4, p+8) at all, primality by sympy;
  F. the vacuous certificate inequality and the two n = 3 witnesses by sympy;
  G. Lean: recompile the lane's scratch file with #print axioms appended for
     EVERY declaration (the lane printed 5 of 22); count Catalan certificates
     in verification.json; confirm the cited package names exist in the package
     source files (not the README).
Usage: python3 [-O] this.py > ../../05-knowledge/results/collatz_mod6_20260922_w6_lean_paste_audit_w6_audit.out
"""
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
import sympy

T0 = time.time()
ROOT = "/tmp/math-wt-collatz-mod6-b"
LEAN_FILE = os.path.join(ROOT, "04-computation/lean/standalone/collatz_mod6_20260922_w6_lean_paste_audit.lean")
PKG = os.path.join(ROOT, "04-computation/lean/CollatzBlueprintAudit")
SCRATCH = "/private/tmp/claude-501/-Users-e-Documents-GitHub-math/e197ec98-d8f9-4475-947b-5af87889cf35/scratchpad"


def is_sq(m):
    r = math.isqrt(m)
    return r * r == m


def edges(n):
    return [(x, y) for x in range(1, n + 1) for y in range(x + 1, n + 1) if is_sq(x + y)]


def uf_components(n):
    parent = list(range(n + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for x, y in edges(n):
        parent[find(x)] = find(y)
    groups = {}
    for v in range(1, n + 1):
        groups.setdefault(find(v), []).append(v)
    return sorted(groups.values())


def adjmat(n):
    A = np.zeros((n + 1, n + 1), dtype=np.int64)
    for x, y in edges(n):
        A[x, y] = A[y, x] = 1
    return A


print("# collatz_mod6_20260922_w6_lean_paste_audit_w6_audit")
print("# independent recomputation for lane lean_paste_audit_w6")
print()

# ---------------------------------------------------------------- A
print("## A self-loops of the pasted adjacency (2m = k^2 <=> m = 2 j^2)")
closed = [2 * j * j for j in range(1, 5)]
brute = [m for m in range(1, 33) if is_sq(2 * m)]
print("closed form 2j^2, j=1..4:", closed, "; brute force m<=32:", brute, "; equal:", closed == brute)
if closed != brute:
    raise RuntimeError("loop census mismatch")
print("pasted Fin-14 indices with Adj i i (value i+1):", [i for i in range(14) if is_sq(2 * (i + 1))],
      "-> values", [i + 1 for i in range(14) if is_sq(2 * (i + 1))])
print("boundary: value 1 (index 0): 1+1 = 2 square?", is_sq(2), "; so the minimal loop is value 2 (index 1)")
print()

# ---------------------------------------------------------------- B
print("## B components of Q_n by union-find, n = 0..32")
ncomp = {}
for n in range(0, 33):
    cs = uf_components(n)
    ncomp[n] = len(cs)
    tag = "  comps=%s" % (cs,) if n in (12, 13, 14) else ""
    print("n=%2d  #comp=%d  edges=%d%s" % (n, len(cs), len(edges(n)), tag))
lane_expect = {1: 1, 2: 2, 3: 2, 13: 2}
lane_expect.update({n: 3 for n in range(4, 13)})
lane_expect.update({n: 1 for n in range(14, 33)})
ok = all(ncomp[n] == lane_expect[n] for n in lane_expect)
print("matches lane S2 / session lead (1,2,2; 3 for 4..12; 2 at 13; 1 for 14..32):", ok, "; n=0 has 0 components (empty graph)")
if not ok:
    raise RuntimeError("component census mismatch")
c12 = uf_components(12)
print("Q_12 components and the summand triple {1,4,6}:", [(c, [v for v in (1, 4, 6) if v in c]) for c in c12])
print("13's neighbours in Q_13:", [y for y in range(1, 13) if is_sq(13 + y)], "; 14's in Q_14:", [y for y in range(1, 14) if is_sq(14 + y)])
print()

# ---------------------------------------------------------------- C
print("## C degree <= 1 vertices via numpy adjacency matrix, n = 14..32")
lowdeg = {}
for n in range(14, 33):
    A = adjmat(n)
    deg = A.sum(axis=1)
    lowdeg[n] = [v for v in range(1, n + 1) if deg[v] <= 1]
    print("n=%2d  deg<=1: %s" % (n, lowdeg[n]))
lead = {18: [16, 17, 18], 19: [16, 18], 31: [], 32: []}
lead.update({n: [18] for n in range(20, 31)})
okC = all(lowdeg[n] == lead[n] for n in lead) and lowdeg[14] == [8, 9, 10] and lowdeg[15] == [8, 9]
print("matches session lead / lane S3:", okC)
if not okC:
    raise RuntimeError("degree census mismatch")
print()

# ---------------------------------------------------------------- D
print("## D Hamiltonian paths: bitmask DP (n <= 20) and an independent backtracker (n <= 32)")


def ham_dp(n):
    """dp[mask] = set of possible last vertices (as bitmask) of a simple path covering mask."""
    if n == 0:
        return False
    nb = [0] * n
    for x, y in edges(n):
        nb[x - 1] |= 1 << (y - 1)
        nb[y - 1] |= 1 << (x - 1)
    full = (1 << n) - 1
    dp = [0] * (1 << n)
    for v in range(n):
        dp[1 << v] = 1 << v
    for mask in range(1, 1 << n):
        ends = dp[mask]
        if not ends:
            continue
        e = ends
        while e:
            v = (e & -e).bit_length() - 1
            e &= e - 1
            cand = nb[v] & ~mask
            while cand:
                w = (cand & -cand).bit_length() - 1
                cand &= cand - 1
                dp[mask | (1 << w)] |= 1 << w
    return dp[full] != 0


def ham_bt(n):
    """Second backtracker: every start vertex, neighbours in DESCENDING degree order,
    prune when an unused vertex other than the current end has < 1 unused neighbour
    (or when >= 2 unused vertices have exactly 1, beyond the one that can be the final end)."""
    nb = {v: [y for y in range(1, n + 1) if y != v and is_sq(v + y)] for v in range(1, n + 1)}
    if n == 1:
        return [1]
    sys.setrecursionlimit(10000)

    def prune(used, last):
        ones = 0
        for v in range(1, n + 1):
            if v in used:
                continue
            free = sum(1 for w in nb[v] if w not in used or w == last)
            if free == 0:
                return True
            if free == 1:
                ones += 1
        return ones > 2  # at most: one adjacent to last (next step) ... conservative bound

    def rec(path, used):
        if len(path) == n:
            return list(path)
        last = path[-1]
        if prune(used, last):
            return None
        for w in sorted(nb[last], key=lambda z: -len(nb[z])):
            if w in used:
                continue
            used.add(w)
            path.append(w)
            r = rec(path, used)
            if r:
                return r
            path.pop()
            used.discard(w)
        return None

    for s in range(1, n + 1):
        r = rec([s], {s})
        if r:
            return r
    return None


def valid_path(p, n):
    return sorted(p) == list(range(1, n + 1)) and all(a != b and is_sq(a + b) for a, b in zip(p, p[1:]))


dp_yes = [n for n in range(1, 21) if ham_dp(n)]
print("DP (n<=20) YES:", dp_yes)
bt = {}
for n in range(1, 33):
    p = ham_bt(n)
    bt[n] = p
    if p is not None and not valid_path(p, n):
        raise RuntimeError("backtracker returned an invalid path at n=%d" % n)
bt_yes = [n for n in bt if bt[n] is not None]
print("backtracker (n<=32) YES:", bt_yes)
print("witnesses 21..32:")
for n in range(21, 33):
    print("  n=%2d %s" % (n, ",".join(map(str, bt[n])) if bt[n] else "none"))
expected = [1, 15, 16, 17, 23] + list(range(25, 33))
if dp_yes != [n for n in expected if n <= 20] or bt_yes != expected:
    raise RuntimeError("Hamiltonian census mismatch: dp=%s bt=%s" % (dp_yes, bt_yes))
print("both agree with lane S4 / session lead (1,15,16,17,23,25..32):", True)
print("n=1 is the one-vertex convention (A090461 starts at 15); n=2: 1+2=3 not a square -> no path")
lead15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
lead23 = [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22]
print("lead paths valid:", valid_path(lead15, 15), valid_path(lead23, 23))
print("squares used: n=15", sorted(set(a + b for a, b in zip(lead15, lead15[1:]))), "; n=23",
      sorted(set(a + b for a, b in zip(lead23, lead23[1:]))))
A15 = adjmat(15)
print("Q_15 degrees:", {v: int(A15[v].sum()) for v in range(1, 16)})
print("deg(4) in Q_15:", int(A15[4].sum()), "neighbours", [y for y in range(1, 16) if A15[4, y]],
      "; position of 4 in lead path:", lead15.index(4), "of 0..14 (interior)")
print("forced endpoints of any Q_15 Hamiltonian path = its degree-1 vertices:", [v for v in range(1, 16) if A15[v].sum() == 1])
print("leaf counts (deg<=1) at n=14,18:", len(lowdeg[14]), len(lowdeg[18]), "(>= 3 leaves forbids a Hamiltonian path; n=14 also fails by leaf count, lane note only said 18)")
print()

# ---------------------------------------------------------------- E
print("## E the Delta = 4 ladder")
ladder = [3, 7, 11, 17]
print("differences:", [b - a for a, b in zip(ladder, ladder[1:])], "; sympy primality of p+4:", [(p + 4, sympy.isprime(p + 4)) for p in ladder])
print("factorint 15, 21:", dict(sympy.factorint(15)), dict(sympy.factorint(21)))
print("PROVED (mod 3): among p, p+4, p+8 the residues mod 3 are p, p+1, p+2, so one is divisible by 3;",
      "a prime triple with difference 4 must contain 3, hence is (3,7,11) and nothing else, for ALL sizes (not only below 100)")
print("check: prime triples (p,p+4,p+8) with p < 10^5:", [(p, p + 4, p + 8) for p in sympy.primerange(2, 10 ** 5) if sympy.isprime(p + 4) and sympy.isprime(p + 8)])
print()

# ---------------------------------------------------------------- F
print("## F the pasted certificate inequality and the n=3 witnesses (sympy)")
n = sympy.symbols("n", positive=True, integer=True)
print("3^0*n + 0 < 2^1*n  <=>  n < 2n  <=>  0 < n:", sympy.simplify(2 * n - (n + 0)) == n)


def orbit(f, t, x):
    out = [x]
    for _ in range(t):
        x = f(x)
        out.append(x)
    return out


cp = orbit(lambda m: m // 2 if m % 2 == 0 else 3 * m + 1, 6, 3)
cm = orbit(lambda m: m // 2 if m % 2 == 0 else 3 * m - 1, 3, 3)
print("plus orbit t=6:", cp, "; 2^4*%d = %d = 3^2*3+5 = %d < 2^4*3 = %d" % (cp[-1], 16 * cp[-1], 27 + 5, 48))
print("minus orbit t=3:", cm, "; 2^3*%d = %d = 3^1*3+7 = %d < 2^3*3 = %d" % (cm[-1], 8 * cm[-1], 9 + 7, 24))
if 16 * cp[-1] != 32 or 8 * cm[-1] != 16:
    raise RuntimeError("witness arithmetic mismatch")
print("plus witness is also the compression_and_lean_audit S15 witness (same t, K, L, B); minus witness is new to this lane (scratch only)")
print()

# ---------------------------------------------------------------- G
print("## G Lean: recompile with #print axioms on every declaration")
src = open(LEAN_FILE, encoding="utf-8").read()
names = []
for line in src.splitlines():
    if line.startswith("theorem ") or line.startswith("def "):
        names.append(line.split()[1])
names = [nm for nm in names if nm not in ("isSq", "adjPaste", "adj", "stepR", "closure", "reach", "connected",
                                          "numComponents", "isPrime", "iterate", "collatz", "DescentCertificate",
                                          "collatzMinus", "DescentCertificateMinus")]
print("declarations audited:", len(names))
os.makedirs(SCRATCH, exist_ok=True)
tmp = os.path.join(SCRATCH, "w6_lean_paste_audit_axioms.lean")
body = src.split("#print axioms")[0]
with open(tmp, "w", encoding="utf-8") as fh:
    fh.write(body + "\n" + "\n".join("#print axioms %s" % nm for nm in names) + "\n")
r = subprocess.run(["lake", "env", "lean", tmp], cwd=PKG, capture_output=True, text=True)
print("exit code:", r.returncode, "; stderr lines:", len(r.stderr.splitlines()))
axiom_sets = {}
for line in r.stdout.splitlines():
    print("  " + line)
    if "depends on axioms:" in line:
        nm = line.split("'")[1]
        axiom_sets[nm] = line.split("axioms:")[1].strip()
    elif "does not depend on any axioms" in line:
        nm = line.split("'")[1]
        axiom_sets[nm] = "[]"
if r.returncode != 0:
    raise RuntimeError("lean recompile failed")
allowed = {"[]", "[propext, Quot.sound]", "[propext]", "[Quot.sound]"}
bad = {k: v for k, v in axiom_sets.items() if v not in allowed}
print("all %d declarations within {propext, Quot.sound}:" % len(axiom_sets), not bad, bad if bad else "")
if bad or len(axiom_sets) != len(names):
    raise RuntimeError("axiom audit failed: %s / %d of %d" % (bad, len(axiom_sets), len(names)))
print("axiom-free declarations:", sorted(k for k, v in axiom_sets.items() if v == "[]"))
code = src
while "/-" in code:
    a = code.index("/-"); b = code.index("-/", a) + 2
    code = code[:a] + code[b:]
code = "\n".join(ln.split("--")[0] for ln in code.splitlines())
print("sorry/axiom/native_decide/import in code:", code.count("sorry"), code.count("\naxiom "), code.count("native_decide"), code.count("\nimport "))
# cited package names
basic = open(os.path.join(PKG, "CollatzBlueprintAudit/Basic.lean"), encoding="utf-8").read()
clock = open(os.path.join(PKG, "CollatzBlueprintAudit/ClockAudit.lean"), encoding="utf-8").read()
readme = open(os.path.join(PKG, "README.md"), encoding="utf-8").read()
print("descent_iff_affine_inequality in Basic.lean:", "theorem descent_iff_affine_inequality" in basic,
      "; GlobalTallyMargin in ClockAudit.lean:", "def GlobalTallyMargin" in clock,
      "; either name in README.md:", ("descent_iff_affine_inequality" in readme) or ("GlobalTallyMargin" in readme))
vj = json.load(open(os.path.join(ROOT, "04-computation/lean/CatalanEllipticAudit/verification.json"), encoding="utf-8"))
print("CatalanEllipticAudit verification.json: status", vj["status"], "; theorems_audited =", vj["theorems_audited"],
      "; all axiom lists empty:", all(v == [] for v in vj["axioms"].values()) if isinstance(vj["axioms"], dict) else "n/a")
cat_src = open(os.path.join(ROOT, "04-computation/lean/CatalanEllipticAudit/CatalanEllipticAudit/Certificates.lean"), encoding="utf-8").read()
print("theorem count in CatalanEllipticAudit/Certificates.lean:", cat_src.count("\ntheorem "))
for f in ("collatz_mod6_20260922_descent_certificate.lean", "collatz_mod6_20260922_certificate_iff_descent_audit.lean"):
    print("standalone exists:", f, os.path.exists(os.path.join(ROOT, "04-computation/lean/standalone", f)))
print()
if __debug__:
    print("runtime %.1fs" % (time.time() - T0))
print("DONE")
