#!/usr/bin/env python3
"""collatz_mod6_20260922_compression_and_lean_audit_audit.py

Adversarial audit of lane compression_and_lean_audit (wave 2026-09-22).
Every key number of the lane note is recomputed by an INDEPENDENT route:
  A1 8/pi^2: PARI/GP prodeuler + mpmath; coprime odd-pair count by per-row inclusion-exclusion
     (not Moebius over d, not numpy gcd); all-integer count by 2*sum phi - 1.
  A2 H(8/pi^2), H(6/pi^2), the p with H(p)=0.704: double precision + PARI/GP solve.
  A3 pigeonhole; entropy invariance under the suffix reordering (chain rule, numeric check);
     sum_T h(T) = n! 2^{T(n-2)} by PERMUTATION enumeration (n<=5) and a second DP (n=6).
  A4 gcd census along Syracuse and standard orbits.
  A5 Paley F_7: recount 21/14/21/2, both Fano planes, both line-orientation rules, h(T_7)=189 by
     brute permutation enumeration; the lines->cycles direction (identity map exists).
  A6 sidecar: 8/pi^2 is ALSO the squarefree density among odd integers (census).
  A7 descent certificate numerics: actual counters at the first strict-descent time, n<=2000.
  B  Lean: rerun the standalone files, plus a new scratch proving certificate <-> strict descent;
     Mathlib.olean search without depth cap; TournamentH7 state; timed cold rebuild of the package;
     verification.json counts; decide count in CatalanEllipticAudit; external HTTP checks.
  C  citation checks: THM-1370 file collision, prior occurrences of h(T_7)=189.
Explicit raise, no bare assert; < 1 GB, < 2 minutes.
"""
import itertools
import json
import math
import os
import subprocess
import sys
import time
from fractions import Fraction

import mpmath

T0 = time.time()
WT = "/tmp/math-wt-collatz-mod6-b"
STAND = os.path.join(WT, "04-computation/lean/standalone")
PKG = os.path.join(WT, "04-computation/lean/CollatzBlueprintAudit")
RES = os.path.join(WT, "05-knowledge/results")


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


def run(cmd, cwd=None, timeout=300, inp=None):
    t = time.time()
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, input=inp)
        return r.returncode, r.stdout + r.stderr, time.time() - t
    except subprocess.TimeoutExpired:
        return -1, "<TIMEOUT>", time.time() - t
    except FileNotFoundError as e:
        return -2, f"<NOT FOUND: {e}>", time.time() - t


# --------------------------------------------------------------------------
hdr("A1. 8/pi^2 and the coprime odd-pair count, by independent routes")
mpmath.mp.dps = 30
e8 = mpmath.mpf(8) / mpmath.pi ** 2
print("mpmath 8/pi^2 =", mpmath.nstr(e8, 20), "; 6/pi^2 =", mpmath.nstr(mpmath.mpf(6) / mpmath.pi ** 2, 20))
gp_script = r"""
default(realprecision,30);
print("gp 8/Pi^2 = ", 8/Pi^2);
print("gp prodeuler(p=3..2e6, 1-1/p^2) = ", prodeuler(p=3,2000000,1-1/p^2));
print("gp prodeuler - 8/Pi^2 = ", prodeuler(p=3,2000000,1-1/p^2) - 8/Pi^2);
h(p) = -(p*log(p)+(1-p)*log(1-p))/log(2);
print("gp H(8/Pi^2) = ", h(8/Pi^2));
print("gp H(6/Pi^2) = ", h(6/Pi^2));
print("gp H(8/Pi^2) - 0.704 = ", h(8/Pi^2) - 0.704);
print("gp p in (0.75,0.9) with H(p) = 0.704: ", solve(p=0.75,0.9,h(p)-0.704));
print("gp binary entropy of 0.704 bits at p = 8/Pi^2? ", h(8/Pi^2) == 0.704);
"""
code, out, dt = run(["gp", "-q"], inp=gp_script, timeout=120)
print(out.strip())
check(code == 0 and "gp H(8/Pi^2) = 0.7002786534" in out, "gp entropy route")

N = 10_000
M = 2 * N - 1


def odd_mult(d):
    return (M // d + 1) // 2


# smallest-prime-factor sieve to M
spf = list(range(M + 1))
for i in range(2, int(M ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, M + 1, i):
            if spf[j] == j:
                spf[j] = i


def distinct_primes(a):
    ps = []
    while a > 1:
        p = spf[a]
        ps.append(p)
        while a % p == 0:
            a //= p
    return ps


# Row-wise inclusion-exclusion: for each odd a, #odd b <= M coprime to a = sum_{d | rad(a)} mu(d) oddmult(d)
count_ie = 0
for a in range(1, M + 1, 2):
    ps = distinct_primes(a)
    s = 0
    for k in range(len(ps) + 1):
        for sub in itertools.combinations(ps, k):
            d = 1
            for p in sub:
                d *= p
            s += (-1) ** k * odd_mult(d)
    count_ie += s
print(f"row-wise inclusion-exclusion: coprime ordered pairs among {N} x {N} odd integers 1..{M} = {count_ie}")
check(count_ie == 81058757, "odd coprime pair count 81058757")
dens = Fraction(count_ie, N * N)
print("density =", float(dens), "; minus 8/pi^2 =", float(dens) - float(e8), "; 8/pi^2 =", float(e8))
# all-integer count = 2 * sum_{k<=N} phi(k) - 1
phi = list(range(N + 1))
for p in range(2, N + 1):
    if spf[p] == p:
        for j in range(p, N + 1, p):
            phi[j] -= phi[j] // p
sphi = sum(phi[1:N + 1])
print(f"all-integer contrast: 2*sum_{{k<={N}}} phi(k) - 1 = 2*{sphi} - 1 = {2*sphi-1}")
check(2 * sphi - 1 == 60794971, "all-integer coprime count 60794971")
print("lane numbers 81058757 and 60794971: CONFIRMED by independent routes.")

# --------------------------------------------------------------------------
hdr("A2. Binary entropy, double precision cross-check")


def Hd(p):
    return -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


p8 = 8 / math.pi ** 2
print("double: H(8/pi^2) = %.12f ; H(6/pi^2) = %.12f ; H(8/pi^2)-0.704 = %.8f" % (Hd(p8), Hd(6 / math.pi ** 2), Hd(p8) - 0.704))
check(abs(Hd(p8) - 0.700278653411826) < 1e-12, "H(8/pi^2)")
check(abs(Hd(6 / math.pi ** 2) - 0.966124241906) < 1e-11, "H(6/pi^2)")
lo, hi = 0.75, 0.9
for _ in range(200):
    mid = (lo + hi) / 2
    if Hd(mid) > 0.704:
        lo = mid
    else:
        hi = mid
print("double bisection: p with H(p)=0.704 is %.10f" % lo)
check(abs(lo - 0.8087879981) < 2e-10, "p704")
print("lane numbers 0.7002786534 / 0.966124241906 / -0.00372135 / 0.8087879981: CONFIRMED.")

# --------------------------------------------------------------------------
hdr("A3. Coding claim: pigeonhole, entropy invariance under reordering, double counting")
for L in (1, 2, 3, 4, 8, 16, 32, 64):
    check(sum(2 ** k for k in range(L)) == 2 ** L - 1, "strings shorter than L")
print("strings of length < L (including the empty string) = sum_{k<L} 2^k = 2^L - 1: checked L in {1,2,3,4,8,16,32,64}")

# Entropy invariance: for ANY distribution on {0,1}^L, H(sigma(X)) = H(X) for the suffix reordering sigma,
# and the chain rule H(X) = H(suffix) + H(prefix | suffix). Numeric check on a non-uniform source (L=12, m=4).
import random
random.seed(20260922)
L, m = 12, 4
w = [random.random() ** 3 for _ in range(2 ** L)]
Z = sum(w)
P = [x / Z for x in w]


def ent(ps):
    return -sum(p * math.log2(p) for p in ps if p > 0)


HX = ent(P)
Q = {}
for x in range(2 ** L):
    b = format(x, f"0{L}b")
    y = b[L - m:] + b[:L - m]
    Q[y] = Q.get(y, 0.0) + P[x]
HY = ent(Q.values())
# chain rule: H(X) = H(S) + H(prefix | S), S = suffix
PS = {}
for x in range(2 ** L):
    s = format(x, f"0{L}b")[L - m:]
    PS[s] = PS.get(s, 0.0) + P[x]
HS = ent(PS.values())
Hcond = 0.0
for x in range(2 ** L):
    b = format(x, f"0{L}b")
    s = b[L - m:]
    if P[x] > 0:
        Hcond -= P[x] * math.log2(P[x] / PS[s])
print(f"non-uniform source on 2^{L} streams: H(X) = {HX:.9f}, H(reordered) = {HY:.9f}, |diff| = {abs(HX-HY):.2e}")
print(f"chain rule: H(suffix) + H(prefix|suffix) = {HS:.9f} + {Hcond:.9f} = {HS+Hcond:.9f}")
check(abs(HX - HY) < 1e-9 and abs(HS + Hcond - HX) < 1e-9, "entropy invariance / chain rule")
print("=> a suffix lookahead permutes the coordinates; the source entropy is unchanged (PROVED: bijection invariance).")

# sum_T h(T) = n! 2^{T(n-2)}: permutation enumeration for n <= 5 (independent of the lane's DP)
print()
print("  n | sum_T h(T) by permutations | n! 2^{T(n-2)}")
for n in range(2, 6):
    pairs = list(itertools.combinations(range(n), 2))
    perms = list(itertools.permutations(range(n)))
    tot = 0
    for mask in range(1 << len(pairs)):
        arc = set()
        for k, (i, j) in enumerate(pairs):
            arc.add((i, j) if (mask >> k) & 1 else (j, i))
        tot += sum(1 for pi in perms if all((pi[q], pi[q + 1]) in arc for q in range(n - 1)))
    rhs = math.factorial(n) * 2 ** ((n - 1) * (n - 2) // 2)
    check(tot == rhs, f"double counting n={n}")
    print(f"{n:3d} | {tot:26d} | {rhs}")


def hpaths(n, arc):
    # second DP: dp over (last vertex, visited-set) built from paths of increasing length
    cur = {(v, 1 << v): 1 for v in range(n)}
    for _ in range(n - 1):
        nxt = {}
        for (v, S), c in cur.items():
            for w in range(n):
                if not (S >> w) & 1 and (v, w) in arc:
                    key = (w, S | (1 << w))
                    nxt[key] = nxt.get(key, 0) + c
        cur = nxt
    return sum(cur.values())


n = 6
pairs = list(itertools.combinations(range(n), 2))
tot = 0
for mask in range(1 << len(pairs)):
    arc = set()
    for k, (i, j) in enumerate(pairs):
        arc.add((i, j) if (mask >> k) & 1 else (j, i))
    tot += hpaths(n, arc)
print(f"  6 | {tot:26d} | {math.factorial(6) * 2**10}  (second DP)")
check(tot == 737280, "n=6 double counting")
print("T(n-2) = C(n,2) - (n-1) = C(n-1,2) for 3 <= n <= 12:",
      all((n - 1) * (n - 2) // 2 == n * (n - 1) // 2 - (n - 1) == math.comb(n - 1, 2) for n in range(3, 13)))
Hv = Hd(p8)
print("per-node figures T(n-2) H / n at n = 3, 7, 10, 20, 100:",
      ", ".join("%.3f" % ((n - 1) * (n - 2) / 2 * Hv / n) for n in (3, 7, 10, 20, 100)))
check(abs((99 * 98 / 2) * Hv / 100 - 33.971) < 5e-4, "per node at 100")

# --------------------------------------------------------------------------
hdr("A4. gcd census on orbits")


def syr(n):
    m = 3 * n + 1
    while m % 2 == 0:
        m //= 2
    return m


bad = sum(1 for n in range(1, 100_001, 2) if math.gcd(n, syr(n)) != 1)
print("odd n <= 100000: gcd(n, Syr(n)) != 1 in", bad, "cases of 50000")
check(bad == 0, "Syracuse coprimality")
# standard map: non-coprime consecutive pairs are exactly the even steps with n >= 4
nc = [(n, n // 2) for n in range(2, 2001, 2) if math.gcd(n, n // 2) != 1]
print("standard map, n <= 2000: non-coprime consecutive pairs (n, T(n)) are exactly the even n >= 4:",
      nc == [(n, n // 2) for n in range(4, 2001, 2)], "; count", len(nc))
odd_bad = sum(1 for n in range(1, 2001, 2) if math.gcd(n, 3 * n + 1) != 1)
print("standard map, odd n <= 2000: gcd(n, 3n+1) != 1 in", odd_bad, "cases")
check(odd_bad == 0 and nc == [(n, n // 2) for n in range(4, 2001, 2)], "standard map gcd pattern")

# --------------------------------------------------------------------------
hdr("A5. Paley tournament on F_7, independent recount")
QR = {1, 2, 4}
check(QR == {(x * x) % 7 for x in range(1, 7)}, "QR mod 7")
arc = {(i, j) for i in range(7) for j in range(7) if i != j and (j - i) % 7 in QR}
outdeg = [sum(1 for j in range(7) if (i, j) in arc) for i in range(7)]
print("arcs:", len(arc), "; out-degrees:", outdeg)
check(len(arc) == 21 and outdeg == [3] * 7, "regular tournament")
# cyclic triples: exactly those with a directed 3-cycle in some rotation
cyc = set()
for a, b, c in itertools.combinations(range(7), 3):
    if ((a, b) in arc and (b, c) in arc and (c, a) in arc) or ((a, c) in arc and (c, b) in arc and (b, a) in arc):
        cyc.add((a, b, c))
print("cyclic triples:", len(cyc), "; transitive:", 35 - len(cyc), "; regular-tournament formula C(7,3) - 7*C(3,2) =", 35 - 21)
check(len(cyc) == 14, "14 cyclic triples")
d13 = {tuple(sorted((s + d) % 7 for d in (0, 1, 3))) for s in range(7)}
d15 = {tuple(sorted((s + d) % 7 for d in (0, 1, 5))) for s in range(7)}
print("cyc == dev{0,1,3} | dev{0,1,5}:", cyc == d13 | d15, "; disjoint:", not (d13 & d15), "; sizes", len(d13), len(d15))
check(cyc == d13 | d15 and not (d13 & d15), "two Fano planes")
# each is a 2-(7,3,1) design; union a 2-(7,3,2) design
for name, D in (("dev{0,1,3}", d13), ("dev{0,1,5}", d15), ("union", d13 | d15)):
    lam = {sum(1 for t in D if i in t and j in t) for i, j in itertools.combinations(range(7), 2)}
    print(f"  {name}: every pair of points in exactly {sorted(lam)} block(s)")
    check(lam == ({1} if name != "union" else {2}), "design lambda")
r13 = all(((s, (s + 1) % 7) in arc and ((s + 1) % 7, (s + 3) % 7) in arc and ((s + 3) % 7, s) in arc) for s in range(7))
r15 = all(((s, (s + 1) % 7) in arc and ((s + 1) % 7, (s + 5) % 7) in arc and ((s + 5) % 7, s) in arc) for s in range(7))
print("line rule dev{0,1,3}: s->s+1->s+3->s :", r13, " (e_i e_{i+1} = e_{i+3})")
print("line rule dev{0,1,5}: s->s+1->s+5->s :", r15, " (differences 1, 4, 2 all in QR)")
check(r13 and r15, "line orientations")
# direction of the paste's psi: 7 lines of dev{0,1,3} -> 7 of the 14 cyclic triples is the IDENTITY (injective);
# cycles -> lines of one plane is 2-to-1 at best.
print("lines of dev{0,1,3} are themselves cyclic triples (identity map, injective):", d13 <= cyc)
print("so 'the 7 lines map to 7 directed 3-cycle triads' is TRUE for one plane and INCOMPLETE (7 more cycles in dev{0,1,5});")
print("'exactly 21 directed 3-cycles' stays false: 14 cycles, 42 rooted cycles, 21 = arcs = transitive triples.")
h7 = sum(1 for pi in itertools.permutations(range(7)) if all((pi[q], pi[q + 1]) in arc for q in range(6)))
print("h(T_7) by brute force over 5040 permutations =", h7)
check(h7 == 189, "h(T_7) = 189")
print("189 is odd (Redei) and is neither 7 nor 21 (THM-1370-h-spectrum-omits-7-21-all-n.md).")

# --------------------------------------------------------------------------
hdr("A6. Sidecar: 8/pi^2 is also the squarefree density among ODD integers")
LIM = 2_000_000
sf = bytearray([1]) * (LIM + 1)
for q in range(2, int(LIM ** 0.5) + 1):
    for j in range(q * q, LIM + 1, q * q):
        sf[j] = 0
odd_sf = sum(sf[n] for n in range(1, LIM + 1, 2))
all_sf = sum(sf[n] for n in range(1, LIM + 1))
print(f"squarefree odd n <= {LIM}: {odd_sf} of {LIM//2} odd integers, density {odd_sf/(LIM//2):.6f} vs 8/pi^2 = {p8:.6f}")
print(f"squarefree n <= {LIM}: {all_sf}, density {all_sf/LIM:.6f} vs 6/pi^2 = {6/math.pi**2:.6f}; odd squarefree over ALL integers {odd_sf/LIM:.6f} vs 4/pi^2 = {4/math.pi**2:.6f}")
check(abs(odd_sf / (LIM // 2) - p8) < 2e-4, "odd squarefree density")
print("PROVED (Euler product): density of squarefree among odd integers = prod_{p odd}(1-1/p^2) = 8/pi^2 (same product as A1).")

# --------------------------------------------------------------------------
hdr("A7. Descent certificate: actual counters at the first strict-descent time")


def collatz(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def counters(n, t):
    K = L = B = 0
    v = n
    for _ in range(t):
        if v % 2 == 0:
            K += 1
            v //= 2
        else:
            L += 1
            B = 3 * B + 2 ** K
            v = 3 * v + 1
    return K, L, B, v


maxt = 0
worst = None
for n in range(2, 2001):
    v, t = n, 0
    while v >= n:
        v = collatz(v)
        t += 1
    K, L, B, y = counters(n, t)
    check(y == v and 2 ** K * y == 3 ** L * n + B, f"tally identity n={n}")
    check(3 ** L * n + B < 2 ** K * n, f"margin at first descent n={n}")
    for s in range(1, t):
        K2, L2, B2, y2 = counters(n, s)
        check(2 ** K2 * y2 == 3 ** L2 * n + B2 and not (3 ** L2 * n + B2 < 2 ** K2 * n), f"no margin before descent n={n} s={s}")
    if t > maxt:
        maxt, worst = t, (n, K, L, B, y)
print("n = 2..2000: identity 2^K y = 3^L n + B holds at every time; margin holds at the first strict-descent time")
print("  and fails at every earlier positive time (equivalence with y < n, as descent_iff_affine_inequality states).")
print("  largest first-descent time:", maxt, "at (n, K, L, B, y) =", worst)
print("  witnesses: n=3:", counters(3, 6), " n=7:", counters(7, 11))
check(counters(3, 6) == (4, 2, 5, 2) and counters(7, 11) == (7, 4, 73, 5), "lane witnesses")
check(16 * 2 == 32 == 9 * 3 + 5 and 32 < 48 and 128 * 5 == 640 == 81 * 7 + 73 and 640 < 896, "witness arithmetic")

# --------------------------------------------------------------------------
hdr("B. Lean reruns and package state")
code, out, dt = run(["lean", "--version"])
print(out.strip())
files = [
    ("collatz_mod6_20260922_paste_verbatim.lean", [], 1),
    ("collatz_mod6_20260922_paste_unicode.lean", [], 1),
    ("collatz_mod6_20260922_paste_stubbed.lean", [], 5),
    ("collatz_mod6_20260922_paste_stubbed.lean", ["-DautoImplicit=false"], 7),
    ("collatz_mod6_20260922_descent_certificate.lean", [], 0),
    ("collatz_mod6_20260922_certificate_iff_descent_audit.lean", [], 0),
]
for fn, extra, expect_err in files:
    code, out, dt = run(["lean"] + extra + [os.path.join(STAND, fn)], cwd=STAND, timeout=240)
    lines = [ln for ln in out.splitlines() if ln.strip()]
    nerr = sum(1 for ln in lines if ": error" in ln)
    nsorry = sum(1 for ln in lines if "declaration uses `sorry`" in ln)
    print(f"  {fn} {' '.join(extra)}: exit {code}, {nerr} error lines, {nsorry} sorry warnings")
    check(nerr == expect_err and (code == 0) == (expect_err == 0), f"lean rerun {fn} {extra}")
print("  the new scratch certificate_iff_descent_audit.lean proves, in core Lean with no sorry:")
print("    DescentCertificate n <-> exists t, 0 < t and iterate collatz t n < n   (for 0 < n)")
print("  i.e. the S15 honesty clause is machine-checked; the file is a scratch, NOT added to any package.")
src = open(os.path.join(STAND, "collatz_mod6_20260922_certificate_iff_descent_audit.lean")).read()
check("sorry" not in src and "axiom" not in src, "scratch is sorry/axiom free")
print("  scratch contains no 'sorry' and no 'axiom':", "sorry" not in src and "axiom" not in src)

hits = []
for root in (os.path.join(WT, "04-computation/lean"), os.path.expanduser("~/.elan")):
    for dp, dn, fnames in os.walk(root):
        if "Mathlib.olean" in fnames:
            hits.append(dp)
print("  Mathlib.olean anywhere under the lean dir or ~/.elan (no depth cap):", hits if hits else "NONE")
t7 = os.path.join(WT, "04-computation/lean/TournamentH7")
lk = open(os.path.join(t7, "lakefile.toml")).read()
print("  TournamentH7: requires mathlib:", 'name = "mathlib"' in lk, "; has .lake:", os.path.isdir(os.path.join(t7, ".lake")))
check('name = "mathlib"' in lk and not os.path.isdir(os.path.join(t7, ".lake")), "TournamentH7 state")

# timed cold rebuild of the package (lake clean; lake build), restoring the cache for the lane script
code, out, dt = run(["lake", "clean"], cwd=PKG, timeout=120)
check(code == 0, "lake clean")
code, out, dt = run(["lake", "build"], cwd=PKG, timeout=300)
print(f"  CollatzBlueprintAudit cold rebuild: lake clean + lake build exit {code} in {dt:.1f}s; last line: {out.strip().splitlines()[-1] if out.strip() else ''}")
check(code == 0, "cold lake build")
vj = json.load(open(os.path.join(PKG, "verification.json")))
ax = vj["axioms"]
print("  verification.json: theorems =", len(ax), "; distinct axiom sets =", sorted({tuple(v) for v in ax.values()}))
check(len(ax) == 35 and all(set(v) <= {"propext", "Quot.sound"} for v in ax.values()), "35 theorems, axioms")
vj2 = json.load(open(os.path.join(WT, "04-computation/lean/CatalanEllipticAudit/verification.json")))
cert = open(os.path.join(WT, "04-computation/lean/CatalanEllipticAudit/CatalanEllipticAudit/Certificates.lean")).read()
print("  CatalanEllipticAudit: theorems =", len(vj2["axioms"]), "; all axiom-free:", all(v == [] for v in vj2["axioms"].values()),
      "; occurrences of 'decide' in Certificates.lean =", cert.count("decide"))
check(len(vj2["axioms"]) == 12 and cert.count("decide") == 12, "12 decide certificates")
# statement shapes cited in S14
basic = open(os.path.join(PKG, "CollatzBlueprintAudit/Basic.lean")).read()
clock = open(os.path.join(PKG, "CollatzBlueprintAudit/ClockAudit.lean")).read()
for frag in ("def EventuallyReachesOne (f : Nat → Nat) : Prop :=\n  ∀ n, 0 < n → ∃ k, iterate f k n = 1",
             "def EventuallyStrictlyDescends (f : Nat → Nat) : Prop :=\n  ∀ n, 1 < n → ∃ k, 0 < k ∧ iterate f k n < n",
             "EventuallyReachesOne f ↔ EventuallyStrictlyDescends f",
             "y < n ↔ 3 ^ L * n + B < 2 ^ K * n"):
    check(frag in basic, "Basic.lean fragment: " + frag[:40])
for frag in ("def GlobalTallyMargin : Prop :=\n  ∀ n, 1 < n → ∃ t, 0 < t ∧\n    3 ^ oddSteps n t * n + standardCarry n t < 2 ^ evenSteps n t * n",
             "EventuallyReachesOne collatz ↔ GlobalTallyMargin"):
    check(frag in clock, "ClockAudit.lean fragment: " + frag[:40])
print("  S14 paraphrases of reachesOne_iff_strictDescent / descent_iff_affine_inequality / GlobalTallyMargin: match the source.")

# external checks (network may be unavailable; report status either way)
for label, url in (("Kuratowski module page", "https://leanprover-community.github.io/mathlib4_docs/Mathlib/Combinatorics/SimpleGraph/Kuratowski.html"),
                   ("Log.Base module page", "https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/SpecialFunctions/Log/Base.html")):
    code, out, dt = run(["curl", "-s", "-o", "/dev/null", "-m", "20", "-w", "%{http_code}", url], timeout=40)
    print(f"  {label}: HTTP {out.strip() if code == 0 else 'unreachable (' + str(code) + ')'}")
code, out, dt = run(["curl", "-s", "-m", "20", "https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/SpecialFunctions/Log/Base.html"], timeout=40)
print("  Log.Base page mentions 'Real.logb':", "Real.logb" in out, "; mentions 'noncomputable':", "noncomputable" in out)

# --------------------------------------------------------------------------
hdr("C. Citation checks")
thm1370 = sorted(f for f in os.listdir(os.path.join(WT, "01-canon/theorems")) if f.startswith("THM-1370"))
print("files named THM-1370*:", thm1370)
check(len(thm1370) == 2 and any("h-spectrum" in f for f in thm1370), "THM-1370 collision")
hs = open(os.path.join(WT, "01-canon/theorems", [f for f in thm1370 if "h-spectrum" in f][0])).read()
print("h-spectrum file status line:", [ln for ln in hs.splitlines() if ln.startswith("**Status:**")][0][:120])
print("h-spectrum statement contains 'No tournament, on any number of vertices, has exactly 7 or exactly 21':",
      "No tournament, on any number of vertices, has exactly 7 or exactly 21" in hs)
for fn, needle in (("arithmetic_braids2_20260917_squarefree_symmetry.md", "189"),
                   ("arithmetic_braids2_20260917_fano_code.md", "189"),
                   ("collatz_mod6_20260922_paley_fano_octonion_design.out", "h(T_7) = number of directed Hamiltonian paths = 189"),
                   ("collatz_mod6_20260922_block_spectrum_audit.md", "14 cyclic triples")):
    txt = open(os.path.join(RES, fn)).read()
    print(f"  {fn}: contains '{needle}':", needle in txt)
    check(needle in txt, "prior occurrence " + fn)
sq = open(os.path.join(RES, "arithmetic_braids2_20260917_squarefree_symmetry.md")).read()
print("  squarefree_symmetry line with the 189 count:", [ln.strip() for ln in sq.splitlines() if "OCF count is" in ln][0][:160])
print("=> h(T_7) = 189 was already recorded (as H = 189 for the 16 sign-gauge orientations, Paley among them) on 2026-09-17;")
print("   the lane's 'not previously recorded' sidecar claim is a novelty overclaim and is corrected in the note.")
print()
print(f"total runtime {time.time() - T0:.1f}s")
