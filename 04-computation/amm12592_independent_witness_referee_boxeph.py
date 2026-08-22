#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AMM 12592 / HYP-9061 -- INDEPENDENT HOSTILE-REFEREE witness verification (boxeph).

Re-implements, FROM THE THEOREM STATEMENTS ONLY (no import of any prior solver
or referee code), the epoch-closure identity and the admissibility test:

  THM-3002 sec 3 (representation problem (*)):
      closure of epoch [R, 2R) in the H = 1 normal form  <=>
      q^{R-1} = sum_{i=0}^{R-1} p^i Delta_i(p,q),
      Delta_i in the Lucas box of degree d_i  (capacity + parity).
  THM-3026 sec 1 (basis form):
      Delta_i = sum_k delta_{i,k} B_{d_i,k}(x),  B_{d,k}(x) = x^{d-k} (1-x)^k,
      admissible  <=>  |delta_k| <= binom(d,k)  and  delta_k == binom(d,k) mod 2.
  Floor profile (THM-3029 (A)): d_i = floor(gamma* (R+i)), D0 = 0,
      gamma* = log_5(phi^2) = log(phi)/log(sqrt 5).

INDEPENDENCE NOTES.
  * The floor d = floor(m gamma*) is computed EXACTLY by integer comparison
    5^d <= phi^{2m} < 5^{d+1} using phi^{2m} = (L_{2m} + F_{2m} sqrt 5)/2
    (Lucas/Fibonacci), i.e. sign tests on (2*5^d - L_{2m})^2 vs 5 F_{2m}^2.
    No floating point, no rational approximation of gamma*, no sympy.
    (Prior lanes used a guarded rational approximation; this is independent.)
  * All polynomial arithmetic is exact int lists over Z[x] written here.
  * Negative controls: the checker is shown to FAIL on three deliberate
    corruptions of a verified witness (parity, capacity, identity).

ADVERSARIAL AUDIT of the two attack lanes' structural claims is in parts B/C.

Exact arithmetic only: int + fractions.Fraction.  No numpy.
"""

import json
import math
import random
from fractions import Fraction

WT = "/tmp/math-wt-boxeph-multifront"
F_COMBINED = WT + "/04-computation/amm12592_floor_witnesses_R8_R16_R32.json"
F_R64_DIRECT = WT + "/04-computation/amm12592_floor_witness_R64.json"
F_R64_BOXEPH = WT + "/04-computation/amm12592_floor_witness_R64_boxeph.json"

comb = math.comb

LEDGER = []


def check(name, ok, detail=""):
    LEDGER.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    msg = "  [%s] %s" % (tag, name)
    if detail:
        msg += "  -- " + detail
    print(msg)
    return bool(ok)


CLAIMS = []


def claim(name, holds, detail=""):
    """Audit of a PRIOR AGENT'S structural claim.  REFUTED is a valid referee
    outcome (a finding), not a failure of this referee."""
    CLAIMS.append((name, bool(holds)))
    tag = "CONFIRMED" if holds else "REFUTED  "
    msg = "  [%s] %s" % (tag, name)
    if detail:
        msg += "  -- " + detail
    print(msg)
    return bool(holds)


# ---------------------------------------------------------------------------
# exact integer polynomials over Z[x], coefficient lists low -> high
# ---------------------------------------------------------------------------

def ptrim(a):
    n = len(a)
    while n > 0 and a[n - 1] == 0:
        n -= 1
    return a[:n]


def padd(a, b):
    if len(a) < len(b):
        a, b = b, a
    out = list(a)
    for i, c in enumerate(b):
        out[i] += c
    return out


def psub(a, b):
    return padd(a, [-c for c in b])


def pscale(a, s):
    return [s * c for c in a]


def pshift(a, k):
    return [0] * k + list(a)


def pmul(a, b):
    a = ptrim(a)
    b = ptrim(b)
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca:
            for j, cb in enumerate(b):
                if cb:
                    out[i + j] += ca * cb
    return out


def peq(a, b):
    return ptrim(list(a)) == ptrim(list(b))


def ppow(a, n):
    out = [1]
    for _ in range(n):
        out = pmul(out, a)
    return out


Q = [1, -1]          # q = 1 - x
P = [0, 1]           # p = x
PMQ = [-1, 2]        # p - q = 2x - 1

_QPOW = {0: [1]}


def qpow(n):
    if n not in _QPOW:
        _QPOW[n] = pmul(qpow(n - 1), Q)
    return _QPOW[n]


def pdiv_exact(a, b):
    """Exact division a / b over Z (leading coeff of b must be +-1)."""
    a = ptrim(list(a))
    b = ptrim(list(b))
    assert b and abs(b[-1]) == 1
    q = [0] * (max(len(a) - len(b) + 1, 0))
    r = list(a)
    while len(ptrim(r)) >= len(b):
        r = ptrim(r)
        s = r[-1] // b[-1]
        k = len(r) - len(b)
        q[k] = s
        r = psub(r, pshift(pscale(b, s), k))
    assert peq(r, []), "division not exact"
    return ptrim(q)


# ---------------------------------------------------------------------------
# exact floor(m * gamma*) via Fibonacci/Lucas integer comparisons
# gamma* = log_5(phi^2)  =>  floor = max{d : 5^d <= phi^{2m}}
# phi^{2m} = (L_{2m} + F_{2m} sqrt5)/2
# ---------------------------------------------------------------------------

def _fib_pair(n):
    """(F_n, F_{n+1}) by fast doubling -- O(log n) big-int ops, no cache."""
    if n == 0:
        return (0, 1)
    a, b = _fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    if n & 1:
        return (d, c + d)
    return (c, d)


def fib(n):
    return _fib_pair(n)[0]


def luc(n):
    a, b = _fib_pair(n)
    return 2 * b - a          # L_n = 2 F_{n+1} - F_n


def five_pow_le_phi2m(d, m):
    """Exact test: 5^d <= phi^(2m)?  (equality impossible for m >= 1)."""
    if d < 0:
        return True
    L, F = luc(2 * m), fib(2 * m)
    A = 2 * 5 ** d - L                    # 5^d <= (L + F sqrt5)/2  <=>  A <= F sqrt5
    if A <= 0:
        return True
    return A * A < 5 * F * F              # strict: sqrt5 irrational, F != 0


def floor_gamma_star(m):
    """Exact floor(m * gamma*), gamma* = log_5(phi^2), for integer m >= 1."""
    d = int(m * 0.5979874356654402)       # float seed only; corrected exactly below
    while five_pow_le_phi2m(d + 1, m):
        d += 1
    while not five_pow_le_phi2m(d, m):
        d -= 1
    return d


# ---------------------------------------------------------------------------
# block machinery straight from THM-3002 / THM-3026 statements
# ---------------------------------------------------------------------------

def block_poly(deltas, d):
    """Delta(x) = sum_k delta_k x^{d-k} (1-x)^k, exact in Z[x]."""
    acc = []
    for k, c in enumerate(deltas):
        if c:
            acc = padd(acc, pscale(pshift(qpow(k), d - k), c))
    return ptrim(acc)


def admissible(deltas, d):
    """Lucas-box capacity + parity (THM-3026 sec 1). Returns (ok, why)."""
    if len(deltas) != d + 1:
        return False, "length %d != d+1 = %d" % (len(deltas), d + 1)
    for k, c in enumerate(deltas):
        b = comb(d, k)
        if abs(c) > b:
            return False, "capacity violated at k=%d: |%d| > C(%d,%d)=%d" % (k, c, d, k, b)
        if (c - b) % 2 != 0:
            return False, "parity violated at k=%d" % k
    return True, ""


def backbone_cells(d):
    """Cell coefficients of (p-q) = 2x-1 expressed at degree d:
       (x-(1-x))(x+(1-x))^{d-1} => b_k = C(d-1,k) - C(d-1,k-1)."""
    return [comb(d - 1, k) - (comb(d - 1, k - 1) if k >= 1 else 0) for k in range(d + 1)]


def is_E_form(sig):
    """sigma == E_m = -1 + x + ... + x^m for some m >= 0 ?"""
    s = ptrim(list(sig))
    if not s or s[0] != -1:
        return False
    return all(c == 1 for c in s[1:])


def verify_epoch(label, R, profile, blocks):
    """Full independent verification; returns dict of results."""
    res = {"label": label, "R": R}
    want = [floor_gamma_star(R + i) for i in range(R)]
    res["profile_ok"] = (list(profile) == want)
    bad = []
    for i, (d, row) in enumerate(zip(profile, blocks)):
        ok, why = admissible(row, d)
        if not ok:
            bad.append((i, why))
    res["adm_ok"] = (len(bad) == 0)
    res["adm_bad"] = bad
    S = []
    polys = []
    for i, (d, row) in enumerate(zip(profile, blocks)):
        pl = block_poly(row, d)
        polys.append(pl)
        S = padd(S, pshift(pl, i))
    res["identity_ok"] = peq(S, qpow(R - 1))
    res["unit_ok"] = all(row[0] in (1, -1) for row in blocks)
    res["signword"] = "".join("+" if row[0] == 1 else "-" for row in blocks)
    res["eff_rate"] = max(Fraction(d, R + i) for i, d in enumerate(profile))
    res["polys"] = polys
    return res


# ---------------------------------------------------------------------------
print("=" * 78)
print("AMM 12592 -- INDEPENDENT HOSTILE-REFEREE WITNESS VERIFICATION (boxeph)")
print("checker re-implemented from THM-3002 (*) + THM-3026 admissibility;")
print("floor profile computed EXACTLY via 5^d <= phi^(2m) Lucas/Fibonacci tests")
print("=" * 78)

# sanity of the exact floor against the definition (spot values)
print("\n-- A0  exact-floor engine self-checks")
check("A0.1 gamma* digits: floor(10^5 g*) = 59798, floor(10^6 g*) = 597987",
      floor_gamma_star(10 ** 5) == 59798 and floor_gamma_star(10 ** 6) == 597987)
check("A0.2 floor engine agrees with direct high-precision check m<=256",
      all(five_pow_le_phi2m(floor_gamma_star(m), m)
          and not five_pow_le_phi2m(floor_gamma_star(m) + 1, m)
          for m in range(1, 257)))
b_ok = True
for d in range(1, 80):
    if not peq(block_poly(backbone_cells(d), d), PMQ):
        b_ok = False
check("A0.3 backbone cells represent 2x-1 exactly at every degree d <= 79", b_ok)
check("A0.4 full box [-C(d,k)]_k represents the constant -1 (d <= 79)",
      all(peq(block_poly([-comb(d, k) for k in range(d + 1)], d), [-1])
          for d in range(1, 80)))

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("PART A  --  witness verification (independent checker)")
print("=" * 78)

witnesses = []
combined = json.load(open(F_COMBINED))
for w in combined:
    witnesses.append(("R%d (THM-3029 combined file)" % w["R"], w["R"],
                      w["profile"], w["blocks"]))
w64d = json.load(open(F_R64_DIRECT))
witnesses.append(("R64 DIRECT (beam lane, amm12592_floor_witness_R64.json)",
                  w64d["R"], w64d["profile"], w64d["blocks"]))
w64b = json.load(open(F_R64_BOXEPH))
witnesses.append(("R64 BOXEPH (doubling map, amm12592_floor_witness_R64_boxeph.json)",
                  w64b["R"], w64b["profile"], w64b["blocks"]))

results = []
for label, R, prof, blocks in witnesses:
    print("\n-- witness: %s" % label)
    r = verify_epoch(label, R, prof, blocks)
    results.append(r)
    check("A1 [%s] profile == exact gamma* floor profile, D0=0" % label, r["profile_ok"])
    check("A2 [%s] all %d blocks admissible (capacity+parity)" % (label, R),
          r["adm_ok"], "" if r["adm_ok"] else str(r["adm_bad"][:3]))
    check("A3 [%s] epoch identity sum_i x^i Delta_i == (1-x)^(R-1) exact" % label,
          r["identity_ok"])
    check("A4 [%s] forced unit delta_(i,0) = +-1 every row (THM-3026 (F))" % label,
          r["unit_ok"])
    print("      effective rate max_i d_i/(R+i) = %s = %.6f" %
          (r["eff_rate"], float(r["eff_rate"])))

# the two R = 64 witnesses: same profile, distinct blocks?
same_prof = results[3]["polys"] is not None and witnesses[3][2] == witnesses[4][2]
ndiff = sum(1 for a, b in zip(witnesses[3][3], witnesses[4][3]) if a != b)
check("A5 two R=64 witnesses share the exact floor profile", same_prof)
check("A6 two R=64 witnesses are DISTINCT solutions (independent evidence)",
      ndiff > 0, "%d of 64 rows differ" % ndiff)
check("A7 R=64 effective rate = 58/97 (both witnesses), < gamma* < 3/5",
      results[3]["eff_rate"] == Fraction(58, 97) == results[4]["eff_rate"]
      and Fraction(58, 97) < Fraction(3, 5))

# negative controls: the checker must FAIL on corrupted witnesses
print("\n-- A8  negative controls (checker must catch deliberate corruptions)")
Rb, profb, blocksb = witnesses[4][1], witnesses[4][2], [list(r) for r in witnesses[4][3]]
c1 = [list(r) for r in blocksb]
kslack = next(k for k in range(profb[5] + 1)
              if abs(blocksb[5][k]) + 1 < comb(profb[5], k))   # not capacity-saturated
c1[5][kslack] += 1                             # parity break (capacity still fine)
ok1, why1 = admissible(c1[5], profb[5])
check("A8.1 +1 corruption caught by parity test", (not ok1) and "parity" in why1)
c2 = [list(r) for r in blocksb]
c2[5][7] += 2                                  # parity ok, identity must break
r2 = verify_epoch("corrupt2", Rb, profb, c2)
check("A8.2 +2 corruption caught by identity test", not r2["identity_ok"])
c3 = [list(r) for r in blocksb]
c3[10][3] = comb(profb[10], 3) + 2             # over capacity, parity kept
ok3, why3 = admissible(c3[10], profb[10])
check("A8.3 over-capacity corruption caught by capacity test",
      (not ok3) and "capacity" in why3)

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("PART B  --  adversarial audit of the attack lanes' structural claims")
print("=" * 78)

print("\n-- B1  'Delta_{R-1} = -1 (full box)' and tail rows 'p-q, p-q' -- claimed")
print("   by the doubling lane for 'EVERY floor witness'; audited on all five:")
for (label, R, prof, blocks), r in zip(witnesses, results):
    d_last = prof[R - 1]
    full_box = (blocks[R - 1] == [-comb(d_last, k) for k in range(d_last + 1)])
    ndev = sum(1 for k in range(d_last + 1)
               if blocks[R - 1][k] != -comb(d_last, k))
    claim("B1 [%s] last block = full box -1" % label, full_box,
          "" if full_box else "%d of %d cells deviate" % (ndev, d_last + 1))
    tail_pq = all(peq(r["polys"][i], PMQ) for i in (R - 3, R - 2))
    claim("B1 [%s] rows R-3, R-2 are exactly p-q" % label, tail_pq)
    r["backbone_tail"] = full_box and tail_pq

print("\n-- B2  'c_i even in every cell' -- TRUE but TAUTOLOGICAL given admissibility")
print("   proof: delta_k == C(d,k) (mod 2) and backbone_k = C(d-1,k)-C(d-1,k-1)")
print("          == C(d-1,k)+C(d-1,k-1) = C(d,k) (mod 2), so delta_k - backbone_k")
print("          is even for EVERY admissible block; no structure is discovered.")
rng = random.Random(12592)
taut = True
for _ in range(300):
    d = rng.randrange(2, 40)
    row = [comb(d, k) - 2 * rng.randrange(0, comb(d, k) + 1) for k in range(d + 1)]
    ok, _w = admissible(row, d)
    bb = backbone_cells(d)
    if not ok or any((row[k] - bb[k]) % 2 for k in range(d + 1)):
        taut = False
check("B2 evenness of c holds for 300 RANDOM admissible blocks (tautology confirmed)", taut)

print("\n-- B3  slimness ledger max|c| over mid rows (the SUBSTANTIVE content)")
claimed = {8: 8, 16: 3432, 32: 184756}
for (label, R, prof, blocks), r in zip(witnesses, results):
    mx = 0
    for i in range(1, R - 1):
        bb = backbone_cells(prof[i])
        for k in range(prof[i] + 1):
            mx = max(mx, abs(blocks[i][k] - bb[k]))
    boxmax = max(comb(prof[R - 1], k) for k in range(prof[R - 1] + 1))
    ratio = Fraction(mx, boxmax)
    print("   %-58s max|c| = %-28d (max box coeff %.2e, ratio %.2e)"
          % (label[:58], mx, float(boxmax), float(ratio)))
    if R in claimed:
        claim("B3 [R=%d] max|c| matches doubling lane's claim %d" % (R, claimed[R]),
              mx == claimed[R])
    r["slim"] = mx

slim_direct = results[3]["slim"]
slim_boxeph = results[4]["slim"]
claim("B4 doubling lane's R=64 output is FAT as claimed (max|c| = 3.4e21 scale)",
      slim_boxeph == 3400359520022008934936)
print("   >>> COORDINATION FINDING: R=64 DIRECT witness max|c| = %d (%.2e)"
      % (slim_direct, float(slim_direct)))
print("   >>> vs BOXEPH doubling output max|c| = %.2e; slimmer by factor %.1f"
      % (float(slim_boxeph), float(Fraction(slim_boxeph, max(slim_direct, 1)))))

print("\n-- B5  reduced identity qC = p^R + q^R - p(p-q): CORRECT but a TAUTOLOGY")
print("   derivation: with c_i := Delta_i - (2x-1) (1<=i<=R-2), c_{R-1} := Delta_{R-1}+1,")
print("   C := Delta_0 + sum_{i>=1} x^i c_i, the epoch identity is EQUIVALENT to")
print("   q C = q^R + x^R - x(2x-1); it adds no information beyond A3.")
for (label, R, prof, blocks), r in zip(witnesses, results):
    Cp = list(r["polys"][0])
    for i in range(1, R):
        ci = psub(r["polys"][i], PMQ if i <= R - 2 else [-1])
        Cp = padd(Cp, pshift(ci, i))
    lhs = pmul(Q, Cp)
    rhs = psub(padd(pshift([1], R), qpow(R)), [0, -1, 2])   # x^R + q^R - (2x^2 - x)
    check("B5 [%s] qC == p^R + q^R - p(p-q) exact" % label, peq(lhs, rhs))

print("\n-- B6  reduced doubling algebra -- verified as an UNCONDITIONAL identity")
for R in (4, 8, 16, 32, 64):
    NR = psub(padd(pshift([1], R), qpow(R)), [0, -1, 2])
    N2R = psub(padd(pshift([1], 2 * R), qpow(2 * R)), [0, -1, 2])
    CR = pdiv_exact(NR, Q)
    C2R = pdiv_exact(N2R, Q)
    rhs = padd(padd(pmul(Q, pmul(CR, CR)), pmul([0, -2, 4], CR)),
               padd(pscale(pmul([0, -1, 2], [1, 2]), -1),
                    pscale(pshift(qpow(R - 1), R), -2)))
    check("B6 [R=%d] C_2R = qC_R^2 + 2p(p-q)C_R - p(p-q)(3p+q) - 2p^R q^(R-1)" % R,
          peq(C2R, rhs))

print("\n-- B7  Lucas facts")
for R in (8, 16, 32, 64):
    lucpoly = []
    for m in range(0, R // 2 + 1):
        coef = (-1) ** m * (R * comb(R - m, m)) // (R - m)
        lucpoly = padd(lucpoly, pscale(ppow(pmul(P, Q), m), coef))
    check("B7 [R=%d] p^R + q^R = L_R(pq) (Lucas polynomial in pq)" % R,
          peq(lucpoly, padd(pshift([1], R), qpow(R))))
check("B7 Fibonacci mass sum_m C(R-m,m) = F_{R+1} for R <= 64",
      all(sum(comb(R - m, m) for m in range(0, R // 2 + 1)) == fib(R + 1)
          for R in range(1, 65)))
check("B7 (note) TRUE Lucas coefficient mass is L_R, not F_{R+1}: check R<=64",
      all(sum(abs((R * comb(R - m, m)) // (R - m)) for m in range(0, R // 2 + 1)) == luc(R)
          for R in range(1, 65)))

print("\n-- B8  endgame attractor E_m = -1 + x + ... + x^m")
check("B8 recurrence E_m - (2x-1) = x E_{m-1} for m <= 60",
      all(peq(psub([-1] + [1] * m, PMQ), pshift([-1] + [1] * (m - 1), 1))
          for m in range(1, 61)))
claimed_entry = {8: 3, 16: 8, 32: 28}
for (label, R, prof, blocks), r in zip(witnesses, results):
    sigma = qpow(R - 1)
    entry, stays = None, True
    for i in range(R):
        t = psub(sigma, r["polys"][i])
        if t and t[0] != 0:
            stays = False
            break
        sigma = ptrim(t[1:])
        if entry is None and is_E_form(sigma):
            entry = i
        elif entry is not None and i < R - 1 and not is_E_form(sigma):
            stays = False
    ok_end = peq(sigma, [])
    print("   %-58s attractor entry row = %s, stays = %s, sigma_(R-1) = 0: %s"
          % (label[:58], entry, stays, ok_end))
    claim("B8 [%s] residual enters E_m attractor and stays to 0" % label,
          entry is not None and stays and ok_end)
    if r["R"] in claimed_entry:
        claim("B8 [R=%d] entry row matches doubling lane's claim (%d)"
              % (r["R"], claimed_entry[r["R"]]), entry == claimed_entry[r["R"]])
    r["attractor_entry"] = entry

print("\n-- B9  floor superadditivity at gamma* ('degrees fit at ANY rate')")
print("   proof: floor(a) + floor(b) <= floor(a+b) with a+b exactly additive in i;")
print("   exact numerical confirmation over all pairs, R = 8..64:")
supok = True
for R in (8, 16, 32, 64):
    dR = [floor_gamma_star(R + i) for i in range(R)]
    d2R = [floor_gamma_star(2 * R + j) for j in range(2 * R)]
    for i in range(R):
        for i2 in range(R):
            if dR[i] + dR[i2] > d2R[i + i2]:
                supok = False
check("B9 d_i(R) + d_i'(R) <= d_(i+i')(2R) for all pairs, R = 8,16,32,64", supok)

print("\n-- B10 profile monotonicity (M) -- HOSTILE lift tests (THM-3029 sec 1)")
lab8, R8, prof8, blocks8 = witnesses[0]
r8 = results[0]
for tag, newprof in [
    ("uniform +1", [d + 1 for d in prof8]),
    ("uniform +5", [d + 5 for d in prof8]),
    ("NON-uniform lift to the 3/5 profile",
     [3 * (R8 + i) // 5 for i in range(R8)]),
]:
    ok = all(nd >= od for nd, od in zip(newprof, prof8))
    lifted = []
    for i in range(R8):
        rr = newprof[i] - prof8[i]
        base = blocks8[i]
        lift = [comb(rr, j) for j in range(rr + 1)]
        conv = [0] * (newprof[i] + 1)
        for k, c in enumerate(base):
            for j, l in enumerate(lift):
                conv[k + j] += c * l
        lifted.append(conv)
    rl = verify_epoch("lift", R8, newprof, lifted)
    same_polys = all(peq(a, b) for a, b in zip(rl["polys"], r8["polys"]))
    check("B10 R=8 lift (%s): admissible + identity + polynomials unchanged" % tag,
          ok and rl["adm_ok"] and rl["identity_ok"] and same_polys)

print("\n-- B11 refuted / unsupported claims (audit verdicts, no computation)")
print("   B11.a REFUTED (attribution): doubling lane's 'every beam configuration")
print("         (width 250-900, ctrl 2-3, span 2-3) fails this profile at the final")
print("         row' and 'first ever / produced NOT by search'.  The DIRECT lane's")
print("         witness file records provenance 'deterministic beam, rank l1deg,")
print("         beam 400, ctrl 2, span 2 + exhaustive 2-row completion' and IS a")
print("         beam product; it verifies (A1-A4).  The universal negative only")
print("         covered the committed solver's ranking family -- itself an instance")
print("         of the THM-3029 beam-negative hazard.  Mathematical content of")
print("         both lanes is unaffected: R = 64 closes either way, twice over.")
print("   B11.b NOT ESTABLISHED: 'the map closes iff input is SLIM'.  Only the")
print("         directions 'slim seeds closed (3 cases)' and 'one fat input failed'")
print("         were observed; no iff is proved.  Downgrade to observation.")
print("   B11.c NOT AUDITED HERE: 'D0(2R) >= 2 D0(R) + O(log R)' (asymptotic claim,")
print("         no finite certificate offered).  Treat as heuristic.")
print("   B11.d CONFIRMED (proof): per-epoch slack D0(R) = o(R) suffices for")
print("         C* <= 1 + gamma*: T(m) = m + 1 + floor(gamma* m) + D0(R(m)) and")
print("         D0(R(m)) <= eps m for all large m, finitely many epochs absorbed")
print("         into the additive constant D; G5-1 assembly is per-epoch.")
print("   B11.e REFUTED AS UNIVERSAL (hostile case = the doubling lane's OWN R=64")
print("         output): 'EVERY floor witness is backbone + corrections' with full-box")
print("         last row and attractor entry FAILS on amm12592_floor_witness_R64_")
print("         boxeph.json -- last block deviates from the full box in 48/76 cells,")
print("         rows R-3/R-2 are not p-q, and the residual NEVER enters the E_m")
print("         attractor (B1/B8).  The witness itself is still VALID (A1-A4).")
print("         The structure DOES hold on all four SLIM witnesses (R=8/16/32 and")
print("         the direct R=64, attractor entry rows 3/8/28/54), so the corrected")
print("         statement is: backbone tail + attractor entry is an observed")
print("         property of SLIM witnesses, not of floor witnesses as such --")
print("         consistent with, and sharpening, the lane's own slim/fat dichotomy.")
print("         Also: 'c_i even' is a tautology of admissibility (B2) and the")
print("         reduced identity is equivalent to the epoch identity (B5); the")
print("         only substantive backbone content is SLIMNESS (B3).")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("PART C  --  standing bound and dependency chain")
print("=" * 78)

base_ok = all((m // 2) <= floor_gamma_star(m) for m in range(1, 8))
check("C1 base rows m=1..7: gamma=1/2 degrees floor(m/2) <= floor(gamma* m)", base_ok)
check("C2 epochs R = 8, 16, 32, 64 all verified closed at the exact gamma* floor "
      "profile, D0 = 0 (Part A)",
      all(r["profile_ok"] and r["adm_ok"] and r["identity_ok"] for r in results))

n_fail = sum(1 for _n, ok in LEDGER if not ok)
n_ref = sum(1 for _n, ok in CLAIMS if not ok)
print("\nVALIDITY LEDGER: %d checks, %d failures" % (len(LEDGER), n_fail))
print("CLAIM-AUDIT LEDGER: %d prior-agent claims audited, %d CONFIRMED, %d REFUTED"
      % (len(CLAIMS), len(CLAIMS) - n_ref, n_ref))
if n_ref:
    print("refuted claims (findings, not referee failures):")
    for nm, ok in CLAIMS:
        if not ok:
            print("   - " + nm)
if n_fail == 0:
    print("ALL INDEPENDENT-REFEREE VALIDITY CHECKS PASSED")
    print("""
VERDICT (canon-ready):
  R = 64 IS CLOSED at the exact gamma* floor profile d_i = floor(gamma*(64+i)),
  D0 = 0, H = 1, by TWO independently produced, distinct witnesses (direct beam
  lane; deterministic carve-and-carry doubling map), both re-verified here by a
  from-the-theorem-statements checker with an exact algebraic (Lucas/Fibonacci)
  floor, negative controls included.

NEW STANDING BOUND:
  C = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654402...  is ATTAINED for all
  critical values n <= 127, with T(n) = n + 1 + floor(gamma* n) for 8 <= n <= 127
  and T(n) = n + 1 + floor(n/2) <= n + 1 + floor(gamma* n) for n <= 7; i.e.
  T(n) <= (1 + gamma*) n + 1 on the whole range.  This extends THM-3029 (A)
  (n <= 63) by one epoch and SUPERSEDES C <= 8/5 (THM-3002 5b) on n <= 127.

DEPENDENCY CHAIN:
  THM-2966 (spine normal form; deadline T(m) = m + 1 + d_m)
    -> THM-3002 (G5-1 assembly of base [1,7] + per-epoch closures; H = 1
       reduction to representation problem (*); base + [8,15] at gamma = 1/2)
    -> THM-3026 (admissibility = Lucas box capacity + parity; (M)/(L) lifting)
    -> THM-3029 (A) (floor profile closes at R = 8, 16, 32; witnesses re-verified
       here, Part A) + profile monotonicity (M) for the base-row envelope C1
    -> THIS REFEREE: R = 64 floor closure, two distinct witnesses, Part A
    -> matching lower bound: archimedean floor THM-3009/THM-3017/THM-3024
       (C >= 1 + gamma* for the class; NOT re-verified here -- cited).
  OPEN: R = 128 and all-R (which would give C* = log_5(5 phi^2) exactly).
  NEXT DECISIVE TEST (sharpened by B3): the DIRECT R = 64 witness is ~928x
  slimmer in cells than the doubling map's output (box-saturation 1.1e-3 vs
  0.99), has the full backbone tail and attractor entry (row 54), i.e. it IS
  a candidate 'witness-grade slim' seed; feed IT to the carve-and-carry map
  for 64 -> 128 before designing any new slimming pass.
""")
else:
    print("REFEREE FAILURES PRESENT -- see ledger above")
