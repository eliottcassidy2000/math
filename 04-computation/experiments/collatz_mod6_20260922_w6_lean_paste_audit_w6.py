#!/usr/bin/env python3
"""Lane lean_paste_audit_w6 (2026-09-22, session collatz-mod6, wave 6).

Audits the pasted Lean block of the wave-6 user message:
  (1) square_sum_graph (Adj x y := exists k, (x+1)+(y+1)=k^2) is NOT loopless;
  (2) square_sum_fourteen_connectivity is TRUE; "merger of the 3 historical
      components" is off (3 components for 4<=n<=12, 2 at n=13, 1 at n=14);
      a core-Lean decide formalization is compiled with the repo toolchain;
  (3) delta_four_prime_step_invariance is TRUE but 15, 21 are composite and
      3,7,11,17 is not an arithmetic progression;
  (4) SignSpecificCertificate.is_positive_sheet : True is vacuous.
Also re-verifies the session lead's Q_n probe (components, degree<=1 vertices,
Hamiltonian paths n<=32) by independent code.
Usage: python3 [-O] this.py > ../../05-knowledge/results/collatz_mod6_20260922_w6_lean_paste_audit_w6.out
"""
import os
import subprocess
import sys
import time
import json
import urllib.request

T0 = time.time()
ROOT = "/tmp/math-wt-collatz-mod6-b"
LEAN_FILE = os.path.join(ROOT, "04-computation/lean/standalone/collatz_mod6_20260922_w6_lean_paste_audit.lean")
PKG = os.path.join(ROOT, "04-computation/lean/CollatzBlueprintAudit")


def is_square(m):
    r = int(m ** 0.5)
    while r * r > m:
        r -= 1
    while (r + 1) * (r + 1) <= m:
        r += 1
    return r * r == m


def qn_adj(n):
    """Square-sum graph on values 1..n, x != y, x+y square."""
    adj = {v: [] for v in range(1, n + 1)}
    for x in range(1, n + 1):
        for y in range(x + 1, n + 1):
            if is_square(x + y):
                adj[x].append(y)
                adj[y].append(x)
    return adj


def components(adj):
    seen = set()
    comps = []
    for v in adj:
        if v in seen:
            continue
        stack = [v]
        comp = []
        seen.add(v)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


def ham_path(adj, n):
    """Backtracking Hamiltonian path search; returns a path or None (exhaustive)."""
    # start from lowest-degree vertices first; endpoints must include all degree-1 vertices
    deg1 = [v for v in adj if len(adj[v]) <= 1]
    if len(deg1) > 2:
        return None
    if any(len(adj[v]) == 0 for v in adj) and n > 1:
        return None
    order = sorted(adj, key=lambda v: len(adj[v]))
    starts = deg1 if deg1 else order
    sys.setrecursionlimit(10000)
    best = [None]

    def rec(path, used):
        if len(path) == n:
            best[0] = list(path)
            return True
        last = path[-1]
        cands = [w for w in adj[last] if w not in used]
        # prune: an unused vertex with all neighbours used and not adjacent to last -> dead
        cands.sort(key=lambda w: sum(1 for z in adj[w] if z not in used))
        for w in cands:
            used.add(w)
            path.append(w)
            if rec(path, used):
                return True
            path.pop()
            used.discard(w)
        return False

    for s in starts:
        if rec([s], {s}):
            return best[0]
    return None


def check_path(path, n):
    if sorted(path) != list(range(1, n + 1)):
        raise RuntimeError("path is not a permutation of 1..n: %r" % (path,))
    for a, b in zip(path, path[1:]):
        if a == b or not is_square(a + b):
            raise RuntimeError("bad consecutive pair %d,%d" % (a, b))
    return True


print("# collatz_mod6_20260922_w6_lean_paste_audit_w6")
print("# lane lean_paste_audit_w6: audit of the pasted Lean block (square-sum graph, Delta=4 ladder, sign certificate)")
print()

# ---------------------------------------------------------------- S1: loopless refutation
print("## S1 pasted adjacency Adj x y := exists k, (x.val+1)+(y.val+1)=k^2 on Fin n, x = y allowed")
loops = [(i, i + 1) for i in range(0, 32) if is_square(2 * (i + 1))]
print("self-loops (index, value) with value <= 32:", loops)
print("minimal witness: index 1, value 2, 2+2 =", 2 + 2, "= 2^2 ->", is_square(4))
print("count of self-loops among values 1..14:", sum(1 for i, v in loops if v <= 14),
      "; values:", [v for i, v in loops if v <= 14])
print("verdict S1: REFUTED as a SimpleGraph (loopless := sorry is unprovable); fix = add x != y")
print()

# ---------------------------------------------------------------- S2: Q_n components / probe
print("## S2 square-sum graph Q_n on values 1..n (x != y): components, session lead probe re-verified")
comp_table = {}
for n in range(1, 33):
    adj = qn_adj(n)
    cs = components(adj)
    comp_table[n] = len(cs)
    extra = ""
    if n in (12, 13, 14):
        extra = "  comps=" + str(cs)
    print("n=%2d  #components=%d  edges=%d%s" % (n, len(cs), sum(len(a) for a in adj.values()) // 2, extra))
probe_ok = all(comp_table[n] == 3 for n in range(4, 13)) and comp_table[13] == 2 and all(comp_table[n] == 1 for n in range(14, 33))
print("session lead probe (3 comps for 4<=n<=12, 2 at 13, connected 14..32):", probe_ok)
if not probe_ok:
    raise RuntimeError("probe mismatch")
print("pasted comment 'N=14 merger of the 3 historical components': at n=13 there are already only 2 ->",
      "OFF by one merger (12->13 merges 3->2 since 13 joins via 13+3=16, 13+12=25; 14 joins via 14+2=16, 14+11=25)")
adj13 = qn_adj(13)
print("neighbours of 13 in Q_13:", adj13[13], "; neighbours of 14 in Q_14:", qn_adj(14)[14])
print()

print("## S3 degree <= 1 vertices of Q_n, n = 14..32 (session lead probe)")
for n in range(14, 33):
    adj = qn_adj(n)
    low = [v for v in adj if len(adj[v]) <= 1]
    print("n=%2d  deg<=1: %s   deg(18)=%s" % (n, low, len(adj[18]) if n >= 18 else "-"))
adj30 = qn_adj(30)
print("18's neighbours in Q_30:", adj30[18], "(18+18=36 is a self-loop, excluded; 18+31=49 appears at n=31)")
adj31 = qn_adj(31)
print("18's neighbours in Q_31:", adj31[18])
probe3 = (sorted(v for v in qn_adj(18) if len(qn_adj(18)[v]) <= 1) == [16, 17, 18]
          and sorted(v for v in qn_adj(19) if len(qn_adj(19)[v]) <= 1) == [16, 18]
          and all([v for v in qn_adj(n) if len(qn_adj(n)[v]) <= 1] == [18] for n in range(20, 31))
          and all([v for v in qn_adj(n) if len(qn_adj(n)[v]) <= 1] == [] for n in (31, 32)))
print("session lead degree probe confirmed:", probe3)
if not probe3:
    raise RuntimeError("degree probe mismatch")
print()

# ---------------------------------------------------------------- S4: Hamiltonian paths
print("## S4 Hamiltonian paths in Q_n, n = 1..32 (exhaustive backtracking; witness printed)")
ham = {}
for n in range(1, 33):
    adj = qn_adj(n)
    if n == 1:
        p = [1]
    else:
        p = ham_path(adj, n)
    ham[n] = p
    if p is not None:
        check_path(p, n)
    print("n=%2d  hamiltonian=%s  %s" % (n, p is not None, ",".join(map(str, p)) if p else ""))
yes = [n for n in ham if ham[n] is not None]
no = [n for n in ham if ham[n] is None]
print("YES:", yes)
print("NO :", no)
expected_yes = [1, 15, 16, 17, 23] + list(range(25, 33))
print("matches session lead / A090461 (1,15,16,17,23,25..32):", yes == expected_yes)
if yes != expected_yes:
    raise RuntimeError("Hamiltonian probe mismatch")
lead15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
lead23 = [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22]
print("session lead path n=15 valid:", check_path(lead15, 15), " n=23 valid:", check_path(lead23, 23))
p15 = lead15
sq15 = sorted(set(a + b for a, b in zip(p15, p15[1:])))
print("squares used by the n=15 lead path:", sq15, "(the pasted 'braiding 9,16,25')")
p23 = lead23
sq23 = sorted(set(a + b for a, b in zip(p23, p23[1:])))
print("squares used by the n=23 lead path:", sq23)
# 'at N=15 the number 4 is the edge to avoid': degree of 4 in Q_15 and its use in the path
adj15 = qn_adj(15)
print("Q_15: deg(4) =", len(adj15[4]), "neighbours", adj15[4], "; 4 sits inside the lead path with neighbours",
      [p15[p15.index(4) - 1], p15[p15.index(4) + 1]], "-> '4 is the edge to avoid' is REFUTED (4 is an interior vertex)")
print("Q_15 endpoints of the lead path:", p15[0], p15[-1], "degrees", len(adj15[8]), len(adj15[9]))
print("Q_15 degree-<=1 vertices:", [v for v in adj15 if len(adj15[v]) <= 1], "(exactly the two forced endpoints 8, 9)")
print("n=18 has 3 leaves [16,17,18] -> no Hamiltonian path by the leaf count alone; n=19..22 and 24 are refuted only by the exhaustive search above")
print()

# ---------------------------------------------------------------- S5: OEIS A090461 (best effort)
print("## S5 OEIS A090461 (network; best effort, failure is non-fatal)")
try:
    req = urllib.request.Request("https://oeis.org/search?q=id:A090461&fmt=json",
                                 headers={"User-Agent": "Mozilla/5.0 (research)"})
    with urllib.request.urlopen(req, timeout=15) as r:
        js = json.loads(r.read().decode())
    res = js if isinstance(js, list) else js.get("results", [])
    if res:
        data = [int(x) for x in res[0]["data"].split(",")]
        print("A090461 name:", res[0]["name"])
        print("A090461 data (first 20):", data[:20])
        print("our YES list agrees with A090461 prefix (n<=32):", data[:len(yes) - 1] == yes[1:] or data[:len(yes)] == yes)
    else:
        print("A090461: empty result")
except Exception as e:  # noqa: BLE001
    print("A090461 fetch failed (non-fatal):", type(e).__name__)
print()

# ---------------------------------------------------------------- S6: Delta = 4 ladder
print("## S6 the 'Delta = 4 prime ladder' 3,7,11,17")


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


ladder = [3, 7, 11, 17]
diffs = [b - a for a, b in zip(ladder, ladder[1:])]
print("differences:", diffs, "-> arithmetic progression:", len(set(diffs)) == 1)
targets = [p + 4 for p in ladder]
print("p+4:", targets, "; membership in {7,11,15,21}:", all(t in {7, 11, 15, 21} for t in targets))
print("primality of p+4:", [(t, is_prime(t)) for t in targets])
print("factorizations: 3*5 =", 3 * 5, "; 3*7 =", 3 * 7)
print("the actual prime AP with difference 4 starting at 3 is only 3,7,11 (15 composite); the longest AP-4 prime run below 100:",
      max(((p, p + 4, p + 8) for p in range(2, 92) if is_prime(p) and is_prime(p + 4) and is_prime(p + 8)), default=None))
print("verdict S6: statement TRUE, 'prime progression' reading REFUTED (15, 21 composite; 17-11 = %d != 4)" % (17 - 11))
print()

# ---------------------------------------------------------------- S7: vacuous sign field
print("## S7 SignSpecificCertificate.is_positive_sheet : True")
print("a field of type True has exactly 1 inhabitant (trivial); it constrains nothing.")
print("with L, K_L, B_L free the remaining field 3^L*n+B_L < 2^K_L*n is satisfiable for every n>=1 by L=0,K_L=1,B_L=0:",
      all(3 ** 0 * n + 0 < 2 ** 1 * n for n in range(1, 1001)), "(n=1..1000 checked)")
print("so the pasted structure is inhabited for every n >= 1 on BOTH sheets: it is not a sign guard nor a descent certificate.")


def collatz(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def collatz_minus(n):
    return n // 2 if n % 2 == 0 else 3 * n - 1


def iterate(f, t, n):
    for _ in range(t):
        n = f(n)
    return n


print("correct minimal shape (compression_and_lean_audit sec.7): exists t K L B, 0<t, 2^K*iterate collatz t n = 3^L*n+B, 3^L*n+B < 2^K*n")
y = iterate(collatz, 6, 3)
print("plus-sheet witness n=3, t=6: orbit", [iterate(collatz, k, 3) for k in range(7)], "y =", y,
      "; K=4 L=2 B=5: 2^4*y =", 2 ** 4 * y, "= 3^2*3+5 =", 9 * 3 + 5, "< 2^4*3 =", 2 ** 4 * 3)
ym = iterate(collatz_minus, 3, 3)
print("minus-sheet witness n=3, t=3: orbit", [iterate(collatz_minus, k, 3) for k in range(4)], "y =", ym,
      "; K=3 L=1 B=7: 2^3*y =", 2 ** 3 * ym, "= 3*3+7 =", 3 * 3 + 7, "< 2^3*3 =", 2 ** 3 * 3)
print("the certificate SHAPE is identical on both sheets (sheet-blind by construction; portrait lane): the sign enters only through the map bound by `iterate`.")
print()

# ---------------------------------------------------------------- S8: Lean run
print("## S8 core-Lean formalization, compiled with the repo toolchain (no Mathlib, no build-target change)")
cmd = ["lake", "env", "lean", LEAN_FILE]
print("cwd:", PKG)
print("cmd:", " ".join(cmd))
t1 = time.time()
r = subprocess.run(cmd, cwd=PKG, capture_output=True, text=True)
t2 = time.time()
print("exit code:", r.returncode)
print("stdout:")
for line in r.stdout.splitlines():
    print("  " + line)
print("stderr lines:", len(r.stderr.splitlines()))
for line in r.stderr.splitlines()[:20]:
    print("  " + line)
src = open(LEAN_FILE, encoding="utf-8").read()
code_lines = [ln.split("--")[0] for ln in src.splitlines()]
# drop the header block comment /- ... -/
code = "\n".join(code_lines)
while "/-" in code:
    a = code.index("/-"); b = code.index("-/", a) + 2
    code = code[:a] + code[b:]
sorry_code = code.count("sorry")
print("source: lines=%d  theorems=%d  'sorry' in code=%d (comment mentions of the paste's sorry: %d)  'axiom' declarations=%d  'native_decide'=%d  imports=%d" % (
    src.count("\n"), src.count("\ntheorem "), sorry_code, src.count("sorry") - sorry_code, code.count("\naxiom "),
    code.count("native_decide"), code.count("\nimport ")))
if r.returncode != 0:
    raise RuntimeError("lean compile failed")
if sorry_code != 0:
    raise RuntimeError("sorry present in scratch file code")
print("toolchain:", open(os.path.join(PKG, "lean-toolchain")).read().strip())
lv = subprocess.run(["lake", "env", "lean", "--version"], cwd=PKG, capture_output=True, text=True)
print("lean --version:", lv.stdout.strip())
print("theorem names:")
for line in src.splitlines():
    if line.startswith("theorem ") or line.startswith("def pastedCertificateAlwaysInhabited"):
        print("  " + line.split(":")[0].strip())
print("(compile wall time not part of the reproducible output)" if not __debug__ else "compile wall time: %.1fs" % (t2 - t1))
print()

# ---------------------------------------------------------------- S9: verdict table
print("## S9 verdict table")
rows = [
    ("square_sum_graph as pasted (Adj := exists k, (x+1)+(y+1)=k^2)", "REFUTED as SimpleGraph", "value 2: 2+2=4; adjPaste_not_irreflexive; fix x != y"),
    ("square_sum_fourteen_connectivity", "TRUE (PROVED by decide)", "connected 14 = true; numComponents 14 = 1"),
    ("'merger of the 3 historical components' at N=14", "REFUTED (off by one)", "numComponents 12 = 3, 13 = 2, 14 = 1"),
    ("Q_13, Q_12 disconnected", "PROVED by decide", "connected 13 = false, connected 12 = false"),
    ("delta_four_prime_step_invariance statement", "TRUE (trivial and decide both close it)", "[3,7,11,17]+4 = [7,11,15,21]"),
    ("'Delta = 4 prime ladder / AP' reading", "REFUTED", "15 = 3*5, 21 = 3*7 composite; 17-11 = 6"),
    ("SignSpecificCertificate.is_positive_sheet : True", "REFUTED as sign guard (vacuous)", "inhabited for all n>=1 with L=0,K=1,B=0"),
    ("minimal correct DescentCertificate", "CITED (compression_and_lean_audit S15)", "sheet-blind shape; witnesses n=3 on both sheets"),
    ("session lead Q_n probe (components, leaves, Hamiltonian 1..32)", "FINITE-EXACT confirmed", "S2-S4"),
    ("'at N=15 the number 4 is the edge to avoid'", "REFUTED", "4 is interior in the lead path (12,4,5); deg(4)=2"),
]
for a, b, c in rows:
    print("| %s | %s | %s |" % (a, b, c))
print()
if __debug__:
    print("runtime %.1fs" % (time.time() - T0))
print("DONE")
