#!/usr/bin/env python3
"""
amm12592_r128_independent_referee_boxeph.py          (independent referee, 2026-08-03)

FRESH independent referee for AMM 12592 epoch-closure witnesses at the gamma* floor
profile.  Written directly from the problem statement (*); shares NO code with any
solver/verifier in this repo (stdlib only: json, fractions, math.comb).

Statement (*), epoch [R, 2R):
    sum_{i=0}^{R-1} x^i Delta_i(x) = (1-x)^{R-1},
    Delta_i = sum_{k=0}^{d_i} delta_{i,k} B_{d_i,k},   B_{d,k}(x) = x^{d-k} (1-x)^k,
    |delta_{i,k}| <= C(d_i,k),   delta_{i,k} == C(d_i,k)  (mod 2),
    d_i = floor(gamma* (R+i)),   gamma* = log(phi)/log(sqrt(5)).

Witness JSON format (inferred from data, confirmed by the identity holding):
    {"R": R, "profile": [d_0..d_{R-1}], "blocks": [[delta_{i,0}..delta_{i,d_i}]]_i}

Per-witness checks:
  STRUCT  len(profile)==len(blocks)==R, len(blocks[i])==d_i+1
  BOX     |delta_{i,k}| <= C(d_i,k) for every cell
  PARITY  delta_{i,k} == C(d_i,k) (mod 2) for every cell
  IDENT   sum_i x^i Delta_i == (1-x)^{R-1} exactly in Z[x] (all higher coeffs 0)
  PROFILE profile[i] == true floor(gamma*(R+i)) by exact integer Fibonacci/Lucas
          sign tests  5^d < phi^(2m) < 5^(d+1), phi^(2m) = (L_{2m}+F_{2m} sqrt5)/2
          (no floats, no rational proxy); slack profile d_i - floor_i recorded (D0)
  EFF     effective rate max_i d_i/(R+i) as exact Fraction, proved < gamma* by the
          same integer test (a/b < gamma*  <=>  5^a < phi^(2b))
  UNIT    forced-cell sanity: Delta_i(0) = delta_{i,d_i} = +-1 every row
Negative controls (per witness): single-cell corruptions must FAIL:
  +2 on an interior in-box cell  -> identity FAILS (box+parity still pass)
  +1 on a cell                   -> parity FAILS
  cell := C(d,k)+2               -> box FAILS

AUDIT (task 2): hostile audit of Angle 4's structural claims, p = x, q = 1-x,
E_m := -1 + x + ... + x^m:
  A1 closed form  C = q^{R-1} - x*E_{R-2}  solves  q*C = p^R + q^R - p(p-q),
     exhaustively R = 2..1024 (incl. odd, prime, 2^t+-1 -- beyond Angle 4's 3..256)
  A2 helper identity (p-q)*(1+...+x^{R-2}) = E_{R-2} + 2x^{R-1}, R = 2..1024;
     master-equation equivalence as an identity in ARBITRARY c_i: for hostile R,
     random integer c_i:  [sum_{i<=R-2} x^i((p-q)+c_i) - x^{R-1}] - q^{R-1}
                       == [sum x^i c_i] - [q^{R-1} - E_{R-1}]
  A3 dyadic parity theorem: q^{R-1} - E_{R-1} == 0 mod 2 coefficientwise IFF R=2^t.
     Exhaustive R = 2..1024 via mod-2 Pascal rows (no Lucas bit-trick, so the Lucas
     step itself is being tested); big-integer binomial cross-check R <= 300.
  A4 instantiation sweep over ALL available witnesses: Delta_{R-1} == -1 (minus
     full box), ballot backbone c_i = Delta_i - (p-q) satisfies
     sum_{i<=R-2} x^i c_i == q^{R-1} - E_{R-1}, and every c-cell is even.
     (Note: c-cell evenness is automatic given witness parity, since the exact
     Bernstein-d decode of (p-q) has cell k = C(d-1,k)-C(d-1,k-1) == C(d,k) mod 2;
     the contentful checks are the last row and the master equation.)
"""

import json
import sys
from fractions import Fraction
from math import comb

BASE = "/tmp/math-wt-boxeph-multifront/04-computation"
OUTDIR = "/tmp/math-wt-boxeph-multifront/05-knowledge/results"

WITNESS_FILES = [
    "amm12592_floor_witness_R128_direct.json",
    "amm12592_floor_witness_R128_lp.json",
    "amm12592_floor_witness_R128_rule_boxeph.json",
    "amm12592_floor_witness_R64_slim.json",
]

LOG_LINES = []


def log(*args):
    s = " ".join(str(a) for a in args)
    LOG_LINES.append(s)
    print(s, flush=True)


# ----------------------------------------------------------------------------
# exact gamma* floor machinery (Fibonacci/Lucas fast doubling, integer signs)
# ----------------------------------------------------------------------------

def fib_lucas(n):
    """(F_n, L_n) by fast doubling, exact integers."""
    def fd(n):
        if n == 0:
            return (0, 1)
        a, b = fd(n >> 1)            # F_k, F_{k+1}
        c = a * (2 * b - a)          # F_{2k}
        d = a * a + b * b            # F_{2k+1}
        if n & 1:
            return (d, c + d)
        return (c, d)
    F, F1 = fd(n)
    return F, 2 * F1 - F


def cmp_5a_phi2b(a, b):
    """sign of 5^a - phi^(2b) (never 0 for b >= 1): phi^(2b) = (L+F*sqrt5)/2."""
    F, L = fib_lucas(2 * b)
    A = 2 * 5 ** a - L
    if A <= 0:
        return -1
    return -1 if A * A < 5 * F * F else 1


def gamma_floor(m):
    """exact floor(gamma* * m), gamma* = log(phi)/log(sqrt5): the unique d with
    5^d < phi^(2m) < 5^(d+1)  (equality impossible: phi^(2m) irrational)."""
    d = (3 * m) // 5                       # crude start, corrected by exact tests
    while cmp_5a_phi2b(d + 1, m) < 0:
        d += 1
    while cmp_5a_phi2b(d, m) > 0:
        d -= 1
    assert cmp_5a_phi2b(d, m) < 0 and cmp_5a_phi2b(d + 1, m) > 0
    return d


# ----------------------------------------------------------------------------
# polynomial helpers: dense int coefficient lists, index = power of x
# ----------------------------------------------------------------------------

def qpow(n):
    """(1-x)^n coefficients."""
    return [(-comb(n, j) if (j & 1) else comb(n, j)) for j in range(n + 1)]


def delta_poly(d, cells):
    """Delta = sum_k cells[k] * x^(d-k) (1-x)^k  as dense coeff list, degree <= d."""
    p = [0] * (d + 1)
    for k, v in enumerate(cells):
        if v:
            base = d - k
            for j in range(k + 1):
                c = comb(k, j)
                p[base + j] += (-v * c) if (j & 1) else (v * c)
    return p


def poly_eq(a, b):
    la, lb = len(a), len(b)
    if la < lb:
        a, b, la, lb = b, a, lb, la
    return a[:lb] == b and all(v == 0 for v in a[lb:])


# ----------------------------------------------------------------------------
# the referee proper
# ----------------------------------------------------------------------------

def verify_witness(R, profile, blocks, check_profile=True, name=""):
    """Full independent verification; returns dict of results."""
    rep = {"name": name, "R": R}
    ok = True

    # STRUCT
    struct_ok = (len(profile) == R and len(blocks) == R and
                 all(len(blocks[i]) == profile[i] + 1 for i in range(R)))
    rep["struct_ok"] = struct_ok
    ok &= struct_ok
    if not struct_ok:
        rep["fail"] = "structure"
        return rep

    # BOX + PARITY + UNIT
    box_ok, parity_ok, unit_ok = True, True, True
    first_box = first_par = None
    for i in range(R):
        d = profile[i]
        row = blocks[i]
        for k in range(d + 1):
            c = comb(d, k)
            v = row[k]
            if abs(v) > c:
                box_ok = False
                if first_box is None:
                    first_box = (i, k, v, c)
            if (v - c) & 1:
                parity_ok = False
                if first_par is None:
                    first_par = (i, k, v, c)
        if row[d] not in (1, -1):
            unit_ok = False
    rep["box_ok"], rep["parity_ok"], rep["unit_ok"] = box_ok, parity_ok, unit_ok
    if first_box:
        rep["first_box_violation"] = first_box
    if first_par:
        rep["first_parity_violation"] = first_par
    ok &= box_ok and parity_ok and unit_ok

    # IDENT
    S = [0] * (R + profile[-1] + 1)
    for i in range(R):
        dp = delta_poly(profile[i], blocks[i])
        for j, cv in enumerate(dp):
            if cv:
                S[i + j] += cv
    identity_ok = poly_eq(S, qpow(R - 1))
    rep["identity_ok"] = identity_ok
    ok &= identity_ok

    # PROFILE (exact gamma* floors) + slack
    if check_profile:
        floors = [gamma_floor(R + i) for i in range(R)]
        slacks = [profile[i] - floors[i] for i in range(R)]
        rep["profile_is_exact_floor"] = all(s == 0 for s in slacks)
        rep["D0_max_slack"] = max(slacks)
        rep["min_slack"] = min(slacks)
        if not rep["profile_is_exact_floor"]:
            rep["first_slack_row"] = next(i for i, s in enumerate(slacks) if s != 0)
        ok &= (rep["min_slack"] >= 0)   # profile must dominate the floor pointwise

    # EFF rate
    eff = max(Fraction(profile[i], R + i) for i in range(R))
    rep["eff_rate"] = eff
    rep["eff_lt_gamma"] = (cmp_5a_phi2b(eff.numerator, eff.denominator) < 0)
    ok &= rep["eff_lt_gamma"]

    rep["VERDICT"] = "PASS" if ok else "FAIL"
    return rep


def corrupt(blocks, i, k, newval):
    b = [list(row) for row in blocks]
    b[i][k] = newval
    return b


def run_negative_controls(R, profile, blocks, name):
    """Three single-cell corruptions; each must FAIL for the stated reason."""
    i0 = R // 2
    d0 = profile[i0]
    k0 = 2                                   # interior cell, box C(d0,2) huge
    v0 = blocks[i0][k0]
    results = {}

    # +2: parity & box preserved, identity must break
    nb = corrupt(blocks, i0, k0, v0 + (2 if v0 + 2 <= comb(d0, k0) else -2))
    r = verify_witness(R, profile, nb, check_profile=False, name=name + "+2")
    results["ctrl_identity"] = (r["VERDICT"] == "FAIL" and not r["identity_ok"]
                                and r["box_ok"] and r["parity_ok"])

    # +1: parity must break
    nb = corrupt(blocks, i0, k0, v0 + 1)
    r = verify_witness(R, profile, nb, check_profile=False, name=name + "+1")
    results["ctrl_parity"] = (r["VERDICT"] == "FAIL" and not r["parity_ok"])

    # over box
    nb = corrupt(blocks, i0, k0, comb(d0, k0) + 2)
    r = verify_witness(R, profile, nb, check_profile=False, name=name + "box")
    results["ctrl_box"] = (r["VERDICT"] == "FAIL" and not r["box_ok"])

    return results


# ----------------------------------------------------------------------------
# AUDIT of Angle 4 structural claims (hostile sweep)
# ----------------------------------------------------------------------------

def E_poly(m):
    """E_m = -1 + x + ... + x^m."""
    return [-1] + [1] * m


def mul_q(c):
    """(1-x)*c."""
    out = list(c) + [0]
    for j in range(len(c)):
        out[j + 1] -= c[j]
    return out


def audit_A1_closed_form(Rmax=1024):
    """C = q^{R-1} - x*E_{R-2} solves q*C = p^R + q^R - p(p-q), R = 2..Rmax."""
    bad = []
    for R in range(2, Rmax + 1):
        C = qpow(R - 1) + [0] * 0
        C = C + [0] * (R - len(C))            # ensure length >= R
        # subtract x*E_{R-2}: coeff x^1 gets -(-1)=+1; x^2..x^{R-1} get -1
        C[1] += 1
        for j in range(2, R):
            C[j] -= 1
        lhs = mul_q(C)
        rhs = [0] * (R + 1)
        for j, c in enumerate(qpow(R)):
            rhs[j] += c
        rhs[R] += 1                            # p^R
        rhs[1] += 1                            # - p(p-q) = -2x^2 + x
        rhs[2] -= 2
        if not poly_eq(lhs, rhs):
            bad.append(R)
    return bad


def audit_A2_master_equation(Rmax=1024, hostile_Rs=(5, 6, 7, 96, 100, 129, 250, 255),
                             trials=3):
    """helper identity sweep + equivalence-as-identity with random c_i."""
    bad_helper = []
    for R in range(2, Rmax + 1):
        lhs = [0] * R
        for i in range(R - 1):                 # (2x-1)*(1+..+x^{R-2})
            lhs[i] -= 1
            lhs[i + 1] += 2
        rhs = E_poly(R - 2) + [0]
        rhs[R - 1] += 2
        if not poly_eq(lhs, rhs):
            bad_helper.append(R)

    import random
    rng = random.Random(12592)
    bad_equiv = []
    for R in hostile_Rs:
        for _ in range(trials):
            cs = [[rng.randint(-9, 9) for _ in range(rng.randint(1, 11))]
                  for _ in range(R - 1)]
            # LHS(*) - RHS(*) with backbone Delta_i = (p-q)+c_i, Delta_{R-1} = -1
            L = [0] * (R + 12)
            for i in range(R - 1):
                L[i] -= 1
                L[i + 1] += 2
                for j, cv in enumerate(cs[i]):
                    L[i + j] += cv
            L[R - 1] -= 1                       # Delta_{R-1} = -1
            for j, cv in enumerate(qpow(R - 1)):
                L[j] -= cv                      # minus RHS(*)
            # claimed equal: sum x^i c_i - (q^{R-1} - E_{R-1})
            M = [0] * (R + 12)
            for i in range(R - 1):
                for j, cv in enumerate(cs[i]):
                    M[i + j] += cv
            for j, cv in enumerate(qpow(R - 1)):
                M[j] -= cv
            for j, cv in enumerate(E_poly(R - 1)):
                M[j] += cv
            if not poly_eq(L, M):
                bad_equiv.append(R)
                break
    return bad_helper, bad_equiv


def audit_A3_dyadic_parity(Rmax=1024, cross_max=300):
    """q^{R-1}-E_{R-1} all-even coefficientwise IFF R = 2^t; Pascal rows mod 2."""
    all_even_set = []
    row = [1]                                   # binomials of n=0 mod 2
    for n in range(0, Rmax):                    # n = R-1
        R = n + 1
        # coeff_j(q^{R-1}-E_{R-1}) mod 2 = C(R-1,j) + 1  (E has all odd coeffs)
        if all(b == 1 for b in row):
            all_even_set.append(R)
        row = [1] + [(row[j - 1] + row[j]) % 2 for j in range(1, n + 1)] + [1]
    powers = [2 ** t for t in range(0, 11) if 2 ** t <= Rmax]
    iff_ok = (all_even_set == powers)

    cross_ok = True
    for R in range(2, cross_max + 1):           # big-int binomial cross-check
        diff_even = all((comb(R - 1, j) + 1) % 2 == 0 for j in range(R))
        # rebuild full coefficient of q^{R-1} - E_{R-1} exactly for extra hostility
        full = qpow(R - 1)
        for j, cv in enumerate(E_poly(R - 1)):
            full[j] -= cv
        really_even = all(c % 2 == 0 for c in full)
        if diff_even != really_even or really_even != ((R & (R - 1)) == 0):
            cross_ok = False
            break
    return all_even_set, iff_ok, cross_ok


def audit_A4_ballot_form(R, profile, blocks):
    """Delta_{R-1} == -1 (minus full box); ballot backbone sums to q^{R-1}-E_{R-1};
    all c-cells even (the latter is automatic given parity)."""
    rep = {}
    last = delta_poly(profile[-1], blocks[-1])
    rep["last_row_is_minus_one"] = (last[0] == -1 and all(v == 0 for v in last[1:]))
    rep["last_row_cells_full_box"] = all(
        blocks[-1][k] == -comb(profile[-1], k) for k in range(profile[-1] + 1))

    # c_i = Delta_i - (p-q); accumulate sum x^i c_i
    Ssum = [0] * (2 * R + max(profile) + 2)
    cells_even = True
    for i in range(R - 1):
        d = profile[i]
        dp = delta_poly(d, blocks[i])
        dp[0] += 1                              # subtract (2x-1)
        dp[1] -= 2
        for j, cv in enumerate(dp):
            Ssum[i + j] += cv
        # cell-level: decode of (p-q) at degree d is C(d-1,k)-C(d-1,k-1)
        for k in range(d + 1):
            ballot = comb(d - 1, k) - (comb(d - 1, k - 1) if k >= 1 else 0)
            if (blocks[i][k] - ballot) & 1:
                cells_even = False
    rep["c_cells_all_even"] = cells_even
    target = qpow(R - 1)
    for j, cv in enumerate(E_poly(R - 1)):
        target[j] -= cv
    rep["master_equation_holds"] = poly_eq(Ssum, target)
    return rep


A4_EXTRA_FILES = [
    "amm12592_floor_witness_R64.json",            # original slim direct-beam
    "amm12592_floor_witness_R64_boxeph.json",     # fat (doubling)
    "amm12592_floor_witness_R64_cellslim.json",   # cellmax-rank direct beam
    "amm12592_floor_witness_R64_rule_boxeph.json",
    "amm12592_floor_witness_R64_slim.json",       # chip-fired slimming of fat
    "amm12592_floor_witness_R128_direct.json",
    "amm12592_floor_witness_R128_lp.json",
    "amm12592_floor_witness_R128_rule_boxeph.json",
]


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    log("=" * 78)
    log("AMM 12592 -- INDEPENDENT REFEREE (fresh implementation from statement (*))")
    log("gamma* floors certified by exact Fib/Lucas integer sign tests; no floats,")
    log("no rational proxy, no imports from any solver in this repo.")
    log("=" * 78)

    verdicts = {}
    for fn in WITNESS_FILES:
        path = f"{BASE}/{fn}"
        with open(path) as fh:
            w = json.load(fh)
        R, profile, blocks = w["R"], w["profile"], w["blocks"]
        log(f"\n--- {fn}  (R={R}) ---")
        log(f"    claimed verified flag: {w.get('verified')}   "
            f"source: {w.get('source', w.get('source_label', '?'))!r}")
        rep = verify_witness(R, profile, blocks, check_profile=True, name=fn)
        log(f"    struct={rep['struct_ok']} box={rep['box_ok']} "
            f"parity={rep['parity_ok']} unit={rep['unit_ok']} "
            f"identity={rep['identity_ok']}")
        log(f"    profile==exact gamma* floor: {rep['profile_is_exact_floor']} "
            f"(D0 max slack {rep['D0_max_slack']}, min slack {rep['min_slack']})")
        eff = rep["eff_rate"]
        log(f"    effective rate max d_i/(R+i) = {eff} "
            f"~ {eff.numerator / eff.denominator:.6f}  (< gamma*: {rep['eff_lt_gamma']})")
        log(f"    VERDICT: {rep['VERDICT']}")

        ctr = run_negative_controls(R, profile, blocks, fn)
        log(f"    negative controls (must all be True): "
            f"identity-break={ctr['ctrl_identity']} parity-break={ctr['ctrl_parity']} "
            f"box-break={ctr['ctrl_box']}")
        rep["controls_ok"] = all(ctr.values())
        rep["eff_rate"] = [eff.numerator, eff.denominator]
        verdicts[fn] = rep

        # checkpoint after every witness
        with open(f"{OUTDIR}/amm12592_r128_independent_referee_verdicts_boxeph.json",
                  "w") as fh:
            json.dump(verdicts, fh, indent=1, default=str)

    # ---------------- Angle 4 audit ----------------
    log("\n" + "=" * 78)
    log("ANGLE 4 STRUCTURAL AUDIT (hostile sweeps, exact)")
    log("=" * 78)

    badA1 = audit_A1_closed_form(1024)
    log(f"A1 closed form C = q^(R-1) - x*E_(R-2) solves qC = p^R+q^R-p(p-q):")
    log(f"   R = 2..1024 exhaustive: {'CONFIRMED (no counterexample)' if not badA1 else f'FAILS at R={badA1[:5]}'}")

    badH, badE = audit_A2_master_equation(1024)
    log(f"A2 helper identity (p-q)(1+..+x^(R-2)) = E_(R-2)+2x^(R-1), R=2..1024: "
        f"{'CONFIRMED' if not badH else f'FAILS at {badH[:5]}'}")
    log(f"   master-equation equivalence (identity in arbitrary c_i), hostile "
        f"R in (5,6,7,96,100,129,250,255) x3 random trials: "
        f"{'CONFIRMED' if not badE else f'FAILS at {badE}'}")

    aeset, iff_ok, cross_ok = audit_A3_dyadic_parity(1024, 300)
    log(f"A3 dyadic parity theorem (q^(R-1)-E_(R-1) even IFF R=2^t):")
    log(f"   all-even set in [1,1024] = {aeset}")
    log(f"   equals powers of two: {iff_ok}; big-int binomial + full-coefficient "
        f"cross-check R<=300: {cross_ok}")

    log("A4 ballot-normal-form instantiation sweep (last row == -1 as poly, "
        "master eq sum x^i c_i == q^(R-1)-E_(R-1), c-cells even):")
    a4_all = {}
    with open(f"{BASE}/amm12592_floor_witnesses_R8_R16_R32.json") as fh:
        for w in json.load(fh):
            a4 = audit_A4_ballot_form(w["R"], w["profile"], w["blocks"])
            a4_all[f"combined_R{w['R']}"] = a4
            log(f"   R{w['R']:>3} (combined slim):        last=-1 "
                f"{a4['last_row_is_minus_one']!s:5}  master-eq "
                f"{a4['master_equation_holds']!s:5}  c-even {a4['c_cells_all_even']}")
    for fn in A4_EXTRA_FILES:
        with open(f"{BASE}/{fn}") as fh:
            w = json.load(fh)
        a4 = audit_A4_ballot_form(w["R"], w["profile"], w["blocks"])
        a4_all[fn] = a4
        log(f"   {fn:<44} last=-1 {a4['last_row_is_minus_one']!s:5}  master-eq "
            f"{a4['master_equation_holds']!s:5}  c-even {a4['c_cells_all_even']}")

    audit = {"A1_bad": badA1, "A2_helper_bad": badH, "A2_equiv_bad": badE,
             "A3_all_even_set": aeset, "A3_iff_ok": iff_ok, "A3_cross_ok": cross_ok,
             "A4": a4_all}

    # ---------------- summary ----------------
    log("\n" + "=" * 78)
    log("SUMMARY")
    log("=" * 78)
    all_pass = True
    for fn, rep in verdicts.items():
        line = (f"{fn}: {rep['VERDICT']}  eff={rep['eff_rate'][0]}/{rep['eff_rate'][1]}"
                f"  floor-exact={rep['profile_is_exact_floor']}"
                f"  controls={rep['controls_ok']}")
        log(line)
        all_pass &= (rep["VERDICT"] == "PASS" and rep["controls_ok"])
    r128_closed = all(verdicts[f]["VERDICT"] == "PASS"
                      for f in WITNESS_FILES if "R128" in f)
    log(f"\nR=128 exact-gamma*-floor epoch closure: "
        f"{'CONFIRMED (three independent witnesses)' if r128_closed else 'NOT confirmed'}")
    log(f"All witnesses + controls pass: {all_pass}")

    with open(f"{OUTDIR}/amm12592_r128_independent_referee_verdicts_boxeph.json",
              "w") as fh:
        json.dump({"witnesses": verdicts, "audit": audit,
                   "r128_closed": r128_closed, "all_pass": all_pass},
                  fh, indent=1, default=str)

    out_text = "\n".join(LOG_LINES) + "\n"
    for p in (f"{BASE}/amm12592_r128_independent_referee_boxeph.out",
              f"{OUTDIR}/amm12592_r128_independent_referee_boxeph.out"):
        with open(p, "w") as fh:
            fh.write(out_text)
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
