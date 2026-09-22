#!/usr/bin/env python3
"""collatz_mod6_20260922_compression_and_lean_audit.py

Lane compression_and_lean_audit (wave 2026-09-22, session collatz-mod6-20260917).

Part A: the pasted "T_{n-2} XOR edge-flip compression with suffix lookahead reaches
        H(8/pi^2) = 0.704 bits per node, absolute zero entropy" claim.
  A1  8/pi^2 = prod_{p odd} (1 - p^-2) is the coprime density among pairs of odd
      integers: Euler-product numerics, an exact Moebius count and an independent
      brute-force gcd count over the first 10^4 x 10^4 odd pairs.
  A2  H(8/pi^2) to 10+ digits (0.70028..., not 0.704); which p would give 0.704.
  A3  Counting refutation: pigeonhole for lossless codes; a suffix lookahead is a
      bijection; the tournament layout costs C(n,2) bits for T(n-2) payload bits;
      re-check of sum_T h(T) = n! 2^{T(n-2)} for n <= 6 (inherited identity);
      the Bernoulli(8/pi^2) entropy is per EDGE and the "per node" figure grows.
  A4  Consecutive odd Syracuse iterates are always coprime (density 1, not 8/pi^2).
  A5  Paley tournament on F_7: 21 arcs, 14 cyclic triples (two Fano planes), 21
      transitive triples, 2 cyclic triples per arc, Hamiltonian path count.
Part B: the pasted Lean 4 "consolidated blueprint".
  B1  Toolchain / package state; lake build of CollatzBlueprintAudit if cached.
  B2  lean on the verbatim paste, the Unicode rendering, a Mathlib-stubbed
      rendering (autoImplicit true and false), and a correct minimal descent
      certificate scratch file. All outputs verbatim.
  B3  Line-by-line inspection table of the paste.

Runs in well under 8 minutes and under 1 GB. Explicit raise, no bare assert.
"""
import math
import os
import subprocess
import sys
import time
from fractions import Fraction

import numpy as np
import mpmath

T0 = time.time()
WT = "/tmp/math-wt-collatz-mod6-b"
LEAN_STANDALONE = os.path.join(WT, "04-computation/lean/standalone")
PKG = os.path.join(WT, "04-computation/lean/CollatzBlueprintAudit")


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# --------------------------------------------------------------------------
hdr("A1. 8/pi^2 as the coprime density among pairs of ODD integers")

mpmath.mp.dps = 30
eight_over_pi2 = mpmath.mpf(8) / mpmath.pi ** 2
print("8/pi^2 =", mpmath.nstr(eight_over_pi2, 20))
print("6/pi^2 =", mpmath.nstr(mpmath.mpf(6) / mpmath.pi ** 2, 20), "(all-integer coprime density, for contrast)")
print("ratio (8/pi^2)/(6/pi^2) = 4/3 = 1/(1-1/4): removing the prime 2 from the Euler product")


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


P = primes_upto(2_000_000)
prod = mpmath.mpf(1)
for p in P[1:]:  # odd primes only
    prod *= (1 - mpmath.mpf(1) / (int(p) * int(p)))
print("prod_{3<=p<=2e6 odd prime} (1-1/p^2) =", mpmath.nstr(prod, 15))
print("  |prod - 8/pi^2| =", mpmath.nstr(abs(prod - eight_over_pi2), 5),
      " (tail of the product over p > 2e6 is O(1/(2e6 log 2e6)) ~ 3e-8)")

# Exact Moebius count. Odd integers 1,3,...,2N-1 (N of them). For odd d, the number of
# odd multiples of d in [1, 2N-1] is ceil(N/d) ... precisely floor((2N-1)/d + 1)/2.
N = 10_000
M = 2 * N - 1


def odd_multiples(d):
    # count of odd m <= M with d | m  (d odd): m = d*k, k odd, k <= M//d
    return (M // d + 1) // 2


# Moebius via sieve up to M
mu = np.ones(M + 1, dtype=np.int64)
is_p = np.ones(M + 1, dtype=bool)
is_p[:2] = False
for i in range(2, M + 1):
    if is_p[i]:
        is_p[i * i::i] = False
        mu[i::i] *= -1
        mu[i * i::i * i] = 0
coprime_moebius = 0
for d in range(1, M + 1, 2):
    if mu[d]:
        c = odd_multiples(d)
        coprime_moebius += int(mu[d]) * c * c
print(f"N = {N} = 10^4 odd integers 1..{M}; ordered pairs = {N*N} = 10^8")
print("CITED: Hardy-Wright, Theorem 332 (coprime density 6/pi^2); odd restriction removes the factor 1-1/4")
print("Moebius exact count of coprime ordered odd pairs =", coprime_moebius)

# Independent brute force with numpy gcd, chunked (memory << 1 GB).
odds = np.arange(1, M + 1, 2, dtype=np.int32)
check(len(odds) == N, "odd count")
brute = 0
CH = 400
for r0 in range(0, N, CH):
    block = np.gcd(odds[r0:r0 + CH, None], odds[None, :])
    brute += int(np.count_nonzero(block == 1))
print("brute-force gcd count of coprime ordered odd pairs   =", brute)
check(brute == coprime_moebius, "Moebius and brute force disagree")
frac = Fraction(coprime_moebius, N * N)
print("density =", coprime_moebius, "/", N * N, "=", float(frac))
print("8/pi^2  =", float(eight_over_pi2))
print("difference density - 8/pi^2 =", float(frac) - float(eight_over_pi2),
      " (Mertens-type error O(log N / N) ~", round(math.log(N) / N, 6), ")")
# also show the all-integer coprime density for contrast at the same size
allints = np.arange(1, N + 1, dtype=np.int32)
brute_all = 0
for r0 in range(0, N, CH):
    block = np.gcd(allints[r0:r0 + CH, None], allints[None, :])
    brute_all += int(np.count_nonzero(block == 1))
print("contrast: coprime pairs among ALL 1..N x 1..N =", brute_all, "density", brute_all / (N * N),
      "vs 6/pi^2 =", float(6 / mpmath.pi ** 2))

# --------------------------------------------------------------------------
hdr("A2. Binary entropy H(8/pi^2) exactly; the paste's 0.704 is wrong")


def H(p):
    p = mpmath.mpf(p)
    return -(p * mpmath.log(p, 2) + (1 - p) * mpmath.log(1 - p, 2))


Hval = H(eight_over_pi2)
print("H(8/pi^2)  =", mpmath.nstr(Hval, 15), "bits  (10 digits: %s)" % mpmath.nstr(Hval, 10))
print("H(6/pi^2)  =", mpmath.nstr(H(mpmath.mpf(6) / mpmath.pi ** 2), 12), "bits (contrast)")
lo, hi = mpmath.mpf("0.75"), mpmath.mpf("0.9")  # H is strictly decreasing on (1/2, 1): bisection
for _ in range(120):
    mid = (lo + hi) / 2
    if H(mid) > mpmath.mpf("0.704"):
        lo = mid
    else:
        hi = mid
p704 = (lo + hi) / 2
print("the p > 1/2 with H(p) = 0.704 is p =", mpmath.nstr(p704, 10), " (not 8/pi^2 =", mpmath.nstr(eight_over_pi2, 10), ")")
print("H(8/pi^2) - 0.704 =", mpmath.nstr(Hval - mpmath.mpf("0.704"), 6))
print("zero entropy would need p in {0,1}; H(p) > 0 for every 0 < p < 1 (strict concavity).")
print("H(8/pi^2) rounded to 3 decimals = %.3f, to 4 = %.4f" % (float(Hval), float(Hval)))

# --------------------------------------------------------------------------
hdr("A3. Counting refutation of the compression claim")

print("Pigeonhole: number of binary strings of length < L is 2^L - 1 < 2^L = number of L-bit streams.")
print(" L | 2^L streams | strings of length < L | deficit")
for L in (1, 2, 3, 4, 8, 16, 32, 64):
    print(f"{L:3d} | {2**L:12d} | {2**L-1:20d} | {2**L-(2**L-1)}")
print("2^64 =", 2**64)
print("=> every injective (lossless) code of all L-bit streams assigns length >= L to at least one stream.")
print("   Prefix-free codes on the uniform source: Kraft + Shannon give mean length >= L (CITED Shannon 1948).")

# Suffix lookahead is a bijection: exhaustively for L = 12, m = 4.
L, m = 12, 4
seen = set()
for x in range(2 ** L):
    bits = format(x, f"0{L}b")
    reordered = bits[L - m:] + bits[:L - m]  # suffix first, then the rest
    seen.add(reordered)
print(f"suffix lookahead (last m={m} bits first) on all 2^{L} = {2**L} streams: {len(seen)} distinct images")
check(len(seen) == 2 ** L, "reordering is not a bijection")
print("=> a bijection composed with any code is still injective; pigeonhole bound unchanged.")

print()
print("Tournament layout: n vertices, C(n,2) arcs (bits); a fixed Hamiltonian spine forces n-1 arcs;")
print("free arcs T(n-2) = C(n,2) - (n-1) = C(n-1,2).  Storage per T(n-2) payload bits = C(n,2) bits.")
print("  n | C(n,2) | T(n-2) | overhead C(n,2)-T(n-2) | payload/storage")
for n in range(3, 13):
    c2 = n * (n - 1) // 2
    t = (n - 2) * (n - 1) // 2
    check(t == c2 - (n - 1), "T(n-2) identity")
    print(f"{n:3d} | {c2:6d} | {t:6d} | {c2 - t:22d} | {t / c2:.4f}")
print("=> the layout EXPANDS the payload by n-1 bits per tournament; ratio -> 1 but never < 1.")

# Re-check the inherited identity sum_T h(T) = n! 2^{T(n-2)} for n <= 6 (scaffolding audit, section 5).
print("inherited: scaffolding audit Theorem 5.1, sum_T h(T) = n! 2^{T(n-2)}; re-checked below for n <= 6")
import itertools


def ham_paths(n, adj):
    # count directed Hamiltonian paths in tournament adj (adj[i][j]=1 iff i->j)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c:
                continue
            for w in range(n):
                if not (S >> w) & 1 and adj[v][w]:
                    dp[S | (1 << w)][w] += c
    return sum(dp[full])


print()
print("  n | 2^C(n,2) tournaments | sum_T h(T) | n! 2^{T(n-2)} | equal")
for n in range(2, 7):
    pairs = list(itertools.combinations(range(n), 2))
    tot = 0
    for mask in range(1 << len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (mask >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        tot += ham_paths(n, adj)
    rhs = math.factorial(n) * 2 ** ((n - 2) * (n - 1) // 2)
    check(tot == rhs, f"double counting identity n={n}")
    print(f"{n:3d} | {2**len(pairs):20d} | {tot:10d} | {rhs:13d} | {tot == rhs}")

print()
print("Bernoulli(8/pi^2) model: if the T(n-2) free arcs were i.i.d. defects with P(flip)=8/pi^2,")
print("the entropy would be H(8/pi^2) = %.5f bits per ARC; total T(n-2) H bits; per node (n-1)(n-2)/(2n) H:" % float(Hval))
print("  n | T(n-2) | total bits | per node")
for n in (3, 4, 7, 10, 20, 100):
    t = (n - 2) * (n - 1) // 2
    print(f"{n:3d} | {t:6d} | {t*float(Hval):10.3f} | {t*float(Hval)/n:8.3f}")
print("=> 'per node' is not a constant; the only correct per-symbol figure is per arc, and it is a model")
print("   assumption (i.i.d. Bernoulli defects) with no defined link to Collatz orbits or to coprimality.")

# --------------------------------------------------------------------------
hdr("A4. Consecutive odd Syracuse iterates are always coprime (density 1, not 8/pi^2)")


def syracuse(n):
    m = 3 * n + 1
    while m % 2 == 0:
        m //= 2
    return m


LIM = 100_000
bad = 0
cnt = 0
for n in range(1, LIM + 1, 2):
    m = syracuse(n)
    cnt += 1
    if math.gcd(n, m) != 1:
        bad += 1
print(f"odd n <= {LIM}: {cnt} pairs (n, Syr(n)); gcd != 1 in {bad} cases")
check(bad == 0, "coprimality of consecutive odd iterates")
print("Proof: gcd(n, 3n+1) = gcd(n, 1) = 1 and Syr(n) | 3n+1, so gcd(n, Syr(n)) = 1. Same for the")
print("unaccelerated map: gcd(n, 3n+1)=1 and gcd(2k, k)=k>1 only on even steps, which are not odd pairs.")
print("=> any 'coprime mask' evaluated on consecutive odd orbit values has density 1, not 8/pi^2 = %.5f." % float(eight_over_pi2))

# --------------------------------------------------------------------------
hdr("A5. Paley tournament on F_7 (session lead probe, verified)")

QR = {1, 2, 4}
arcs = {(i, j) for i in range(7) for j in range(7) if i != j and (j - i) % 7 in QR}
print("arcs:", len(arcs), " (expected 21)")
check(len(arcs) == 21, "21 arcs")
check(all(((i, j) in arcs) != ((j, i) in arcs) for i in range(7) for j in range(7) if i != j), "tournament")
cyc, trans = [], []
for tri in itertools.combinations(range(7), 3):
    a, b, c = tri
    out = {v: sum(1 for w in tri if (v, w) in arcs) for v in tri}
    if sorted(out.values()) == [1, 1, 1]:
        cyc.append(tri)
    else:
        trans.append(tri)
print("cyclic triples:", len(cyc), " transitive triples:", len(trans), " (p^3-p)/24 = (7^3-7)/24 =", (343 - 7) // 24)
check(len(cyc) == 14 and len(trans) == 21, "14/21 split")
dev13 = {tuple(sorted(((s + d) % 7) for d in (0, 1, 3))) for s in range(7)}
dev15 = {tuple(sorted(((s + d) % 7) for d in (0, 1, 5))) for s in range(7)}
print("dev{0,1,3} lines:", sorted(dev13))
print("dev{0,1,5} lines:", sorted(dev15))
print("shorthand: dev{0,1,3} =", " ".join("".join(map(str, t)) for t in sorted(dev13)), "; dev{0,1,5} =", " ".join("".join(map(str, t)) for t in sorted(dev15)))
print("cyclic triples == dev{0,1,3} U dev{0,1,5}:", set(cyc) == dev13 | dev15, "; disjoint:", not (dev13 & dev15))
check(set(cyc) == dev13 | dev15 and not (dev13 & dev15), "two Fano planes")
# orientation on lines {s, s+1, s+3}: s -> s+1 -> s+3 -> s
ok = all(((s, (s + 1) % 7) in arcs and ((s + 1) % 7, (s + 3) % 7) in arcs and ((s + 3) % 7, s) in arcs) for s in range(7))
print("on every line {s,s+1,s+3}: s->s+1->s+3->s (octonion rule e_i e_{i+1} = e_{i+3}):", ok)
check(ok, "line orientation")
per_arc = {a: sum(1 for t in cyc if a[0] in t and a[1] in t) for a in arcs}
print("cyclic triples through each arc: set =", sorted(set(per_arc.values())), " (2-(7,3,2) design)")
check(set(per_arc.values()) == {2}, "lambda = 2")
adj7 = [[1 if (i, j) in arcs else 0 for j in range(7)] for i in range(7)]
h7 = ham_paths(7, adj7)
print("Hamiltonian paths in the Paley tournament T_7: h =", h7)
print("directed 3-cycles as CYCLES = 14; as (line, rotation) pairs 14 x 3 = 42; the paste's 7 x 3 = 21")
print("counts three rotations of each of only 7 of the 14 cycles. Also each Fano plane has 7 lines and 7 points;")
print("the cyclic-triple design has 14 blocks, so 'lines <-> 3-cycle triads' is 2-to-1 at best.")
print("THM-1370 (canon): the Hamiltonian-path spectrum omits 7 and 21; here h(T_7) =", h7, "is neither.")

# --------------------------------------------------------------------------
hdr("B1. Lean toolchain and package state")


def run(cmd, cwd=None, timeout=300):
    t = time.time()
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        out = (r.stdout + r.stderr)
        code = r.returncode
    except subprocess.TimeoutExpired:
        out, code = "<TIMEOUT>", -1
    except FileNotFoundError as e:
        out, code = f"<NOT FOUND: {e}>", -2
    return code, out, time.time() - t


code, out, dt = run(["lean", "--version"])
print("lean --version:", out.strip(), "(exit", code, ")")
code, out, dt = run(["lake", "--version"])
print("lake --version:", out.strip(), "(exit", code, ")")
print("package:", PKG)
print("  lean-toolchain:", open(os.path.join(PKG, "lean-toolchain")).read().strip())
print("  lakefile.toml:", open(os.path.join(PKG, "lakefile.toml")).read().strip().replace("\n", " | "))
mathlib_hits = []
for root in (os.path.join(WT, "04-computation/lean"), os.path.expanduser("~/.elan")):
    for dp, dn, fn in os.walk(root):
        if "Mathlib.olean" in fn:
            mathlib_hits.append(dp)
        if dp.count(os.sep) - root.count(os.sep) > 7:
            dn[:] = []
print("  Mathlib.olean found anywhere under the lean dir or ~/.elan:", mathlib_hits if mathlib_hits else "NONE")
cached = os.path.isdir(os.path.join(PKG, ".lake", "build"))
print("  .lake/build cached:", cached)
if cached:
    code, out, dt = run(["lake", "build"], cwd=PKG, timeout=300)
    print(f"  lake build: exit {code} in {dt:.1f}s")
    print("  " + out.strip().replace("\n", "\n  "))
    check(code == 0, "lake build failed")
else:
    print("  lake build skipped (no cache; the README's verify.py rebuilds from scratch in ~10 s on this machine)")
import json
vj = json.load(open(os.path.join(PKG, "verification.json")))
names = sorted(vj.get("axioms", {}).keys())
print(f"  verification.json lists {len(names)} audited theorems; axiom sets used:",
      sorted({tuple(v) for v in vj["axioms"].values()}))
for nm in names:
    print("   ", nm, vj["axioms"][nm])
cat = os.path.join(WT, "04-computation/lean/CatalanEllipticAudit")
vj2 = json.load(open(os.path.join(cat, "verification.json")))
n2 = sorted(vj2.get("axioms", {}).keys())
print(f"  CatalanEllipticAudit/verification.json lists {len(n2)} theorems; axiom sets:",
      sorted({tuple(v) for v in vj2["axioms"].values()}))
for nm in n2:
    print("   ", nm, vj2["axioms"][nm])

# README controls of both packages, recomputed here so the note can cite them.
def collatz_std(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1
v, steps = 7, 0
while v != 1:
    v = collatz_std(v); steps += 1
print("  README control: 7 reaches 1 in", steps, "standard steps")
check(steps == 16, "seven in 16")
def signed_std(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1
for start, per in ((-1, 2), (-5, 5), (-17, 18)):
    v = start
    for _ in range(per):
        v = signed_std(v)
    print(f"  README control: {start} returns to itself after {per} signed steps:", v == start)
    check(v == start, "signed return")
print("  Catalan package controls: 2^11 - 3^7 =", 2**11 - 3**7, "; 17*139 =", 17 * 139,
      "; point (10,49) on y^2 = 2x^3 + 4x^2 + 1:", 49**2 == 2 * 10**3 + 4 * 10**2 + 1)
check(2**11 - 3**7 == -139 and 17 * 139 == 2363 and 49**2 == 2 * 10**3 + 4 * 10**2 + 1, "catalan controls")

# --------------------------------------------------------------------------
hdr("B2. lean on the pasted blueprint (three renderings) and on the descent-certificate scratch")

files = [
    ("verbatim ASCII paste", "collatz_mod6_20260922_paste_verbatim.lean", []),
    ("Unicode rendering", "collatz_mod6_20260922_paste_unicode.lean", []),
    ("Mathlib-stubbed rendering, autoImplicit default (true)", "collatz_mod6_20260922_paste_stubbed.lean", []),
    ("Mathlib-stubbed rendering, -DautoImplicit=false (Mathlib project setting)", "collatz_mod6_20260922_paste_stubbed.lean", ["-DautoImplicit=false"]),
    ("correct minimal descent certificate (core Lean, no sorry)", "collatz_mod6_20260922_descent_certificate.lean", []),
]
summary = []
for label, fn, extra in files:
    path = os.path.join(LEAN_STANDALONE, fn)
    code, out, dt = run(["lean"] + extra + [path], cwd=LEAN_STANDALONE, timeout=240)
    lines = [ln for ln in out.splitlines() if ln.strip() and not ln.lstrip().startswith("Hint:")]
    nerr = sum(1 for ln in lines if ": error" in ln)
    nsorry = sum(1 for ln in lines if "declaration uses `sorry`" in ln)
    summary.append((label, fn, code, nerr, nsorry, dt))
    print(f"--- {label}: {fn} {' '.join(extra)}")
    print(f"    exit {code}, {nerr} error lines, {nsorry} sorry warnings, {dt:.1f}s")
    for ln in lines:
        print("    " + ln.replace(LEAN_STANDALONE + "/", ""))
print()
print("summary: label | exit | errors | sorry warnings")
for label, fn, code, nerr, nsorry, dt in summary:
    print(f"  {label} | {code} | {nerr} | {nsorry}")
check(summary[0][2] != 0 and summary[1][2] != 0, "paste should fail")
check(summary[4][2] == 0 and summary[4][3] == 0 and summary[4][4] == 0, "certificate scratch must compile clean")
# Type* is Mathlib notation (Mathlib.Tactic.TypeStar), not core Lean: verify on a one-liner.
tstar = os.path.join(LEAN_STANDALONE, "collatz_mod6_20260922_typestar_probe.lean")
with open(tstar, "w") as fh:
    fh.write("def probe (V : Type*) : Type := V\n")
code, out, dt = run(["lean", tstar], cwd=LEAN_STANDALONE, timeout=60)
print()
print("probe `def probe (V : Type*) : Type := V` in core Lean: exit", code)
for ln in out.splitlines():
    if ln.strip():
        print("    " + ln.replace(LEAN_STANDALONE + "/", ""))
check(code != 0, "Type* should not parse in core Lean")
os.remove(tstar)
print()
print("descent-certificate scratch file (verbatim):")
print(open(os.path.join(LEAN_STANDALONE, "collatz_mod6_20260922_descent_certificate.lean")).read())

# --------------------------------------------------------------------------
hdr("B3. Line-by-line inspection of the paste (Unicode rendering line numbers)")

src = open(os.path.join(LEAN_STANDALONE, "collatz_mod6_20260922_paste_unicode.lean")).read().splitlines()
for i, ln in enumerate(src, 1):
    print(f"{i:3d}: {ln}")
print()
inspection = [
    (1, "import Mathlib.Data.Matrix.Basic", "no Mathlib build available: 'unknown module prefix Mathlib' (E1: missing dependency; aborts the file)"),
    (2, "import Mathlib.Combinatorics.SimpleGraph.Kuratowski", "E2: no such Mathlib module: CITED external check 2026-09-22, https://leanprover-community.github.io/mathlib4_docs/Mathlib/Combinatorics/SimpleGraph/Kuratowski.html returns HTTP 404 (Mathlib has no Kuratowski theorem and no SimpleGraph.Planar)"),
    (3, "import Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse", "exists in Mathlib but is unrelated (arcsin/arccos); logb is Real.logb in Mathlib.Analysis.SpecialFunctions.Log.Base (CITED external check 2026-09-22: 'noncomputable def Real.logb (b x : R) : R'), which is not imported"),
    (6, "Type* / [Fintype ...]", "Type* is Mathlib syntax (parse error in core Lean); Fintype is Mathlib; matrices over Z carry no tournament axioms (no antisymmetry, no 0/1 entries): E3 the structure does not say 'tournament'"),
    (15, "simple_graph_of_tournament_union B", "E4: undefined identifier (with autoImplicit=false: 'Unknown identifier'; with autoImplicit=true: 'Function expected')"),
    (16, "G.Planar", "E5: SimpleGraph has no field Planar in Mathlib; 'Invalid field notation'"),
    (16, "by sorry", "E6: sorry (theorem 1 of 4)"),
    (19, "points : Fin 7", "E7 (semantic): a single point of Fin 7, not a point set; the type FanoPlane then has 7 x |lines-data| inhabitants"),
    (21, "forall l in lines, Set.Card l = 3", "E8: 'Unknown constant Set.Card' (Mathlib: Set.ncard / Set.encard)"),
    (23, "paley_3_cycles ... := sorry", "E9: a definition by sorry (its value is unspecified, so nothing about it can be proved)"),
    (26, "exists psi : FanoPlane equiv paley_3_cycles M7, True", "E10 (semantic): an Equiv between the structure TYPE FanoPlane and the subtype of a set of triples, followed by 'True'; the body is provable (with sorry-defined set) or vacuous, and never mentions incidence"),
    (28, "def shannon_entropy_limit : R := ...", "E11: real arithmetic is noncomputable in Mathlib: 'failed to compile definition, consider marking it as noncomputable'"),
    (32, "Matrix (Fin n) (Fin n) Z", "E12: n is unbound ('Unknown identifier n' under Mathlib's autoImplicit=false; silently auto-bound as an implicit otherwise); the statement 'exists layout, true' is trivially true and encodes no compression claim"),
    (32, "by sorry", "E6: sorry (theorem 3 of 4)"),
    (35, "collatz_functor_path step = none", "E13: undefined identifier; the hypotheses n, h, B are unused; the conclusion mentions neither n nor a Collatz map"),
    (35, "by sorry", "E6: sorry (theorem 4 of 4)"),
]
print("line | fragment | error class")
for ln, frag, cls in inspection:
    print(f"{ln:4d} | {frag} | {cls}")
print()
print("the paste section-4 tokens K_5, K_{3,3}, S_2 x S_3, B^3 = -I, omega do not occur anywhere in the Lean file.")
print("error classes: E1 missing dependency; E2 nonexistent module; E3 no tournament axioms; E4/E13 undefined identifiers;")
print("E5 nonexistent field; E6 four sorries; E7/E10 mis-typed Fano statement; E8 unknown constant; E9 sorry-defined data;")
print("E11 noncomputable real def; E12 unbound n.  Count: 4 theorems, 4 sorry, 0 proved statements, 3 undefined names,")
print("1 nonexistent Mathlib module, 1 unknown constant, 1 unbound variable.")

print()
print(f"total runtime {time.time() - T0:.1f}s")
