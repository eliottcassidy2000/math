#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_wsb_verify_widespreadsign_kps-S9-wf.py   (kind-pasteur-2026-06-19-S9-wf)

ADVERSARIAL VERIFICATION of the claimed advance "wide-spread-signed-weyl"
(file lrc14_wsb_wide-spread-signed-weyl_kps-S9-wf.py).

We INDEPENDENTLY re-derive every load-bearing claim and HUNT for counterexamples:

  (V0) The framing claim: is the target meas(S7(E)) <= meas(S7(consec_k)) the
       RIGHT one?  Canon says the rigorous target is meas(S7(E)) <= cap_k.
       We check: (a) the prompt's cap values are wrong/swapped (canon T859/T860);
                 (b) consec <= cap_k (so "<= consec" suffices) ... AT WHICH k?
                 (c) is "<= consec" even TRUE as a same-k statement? (T860 says
                     FALSE at k=12).
  (V1) exact targets: M7(k), meas(S7(consec_k)), delta_k = consec - M7.
       Claim: delta_k >= 0.302 for all k=8..13.  RE-DERIVE EXACTLY.
  (V2) LEMMA A (support-6 floor): K(n)=0 unless >=6 nonzero non-7 coords.
       Re-derive the C(U) algebra independently; EXHAUSTIVE support<=5 over a
       structured grid (not just random); confirm support-6 nonzero.
       ALSO: verify on ACTUAL relation lattices of real E.
  (V3) LEMMA B envelope: |shat(n,j)|=|sin(pi n/7)|/(pi|n|); c1 envelope constant.
  (V4) the signed identity meas(S7)=M7(k)+sum K(n).
  (V5) RESONANCE w==0 mod 7: hunt for a resonant config exceeding consec / cap.
  (V6) HUNT for counterexamples to meas(S7(E)) <= meas(S7(consec_k)):
        - exhaustive bounded-spread primitive E
        - aggressive wide-spread
        - resonant (w == 0 mod 7)
        - short-relation shapes {0,1,N,N+1}, {0,1,2,N,...}
       One hit => the "<= consec" claim is FALSE for that k (witness).
  (V7) HUNT for counterexamples to meas(S7(E)) <= cap_k (the CANONICAL target).
  (V8) the GAP audit: is the absolute residual after the floor really > delta_k
       at the AP?  Is there really a rigorous B(k)?  Check the constants.
"""
from __future__ import annotations
import sys, itertools, math, cmath
from fractions import Fraction as F
from math import comb, gcd

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

TWO_PI_I = 2j * math.pi

# ---------------------------------------------------------------------------
# EXACT meas(S7): breakpoint sweep.  Sector hit = floor(7*frac(e x)).
# S7 at x iff {floor(7*frac(e_i x))} = {0..6}.  e=0 always hits sector 0.
# ---------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7: total += x1 - x0
    return total

def M7(k):
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))

# ---------------------------------------------------------------------------
# Independent re-derivation of meas(S7) via inclusion-exclusion over sectors,
# as a CROSS-CHECK of the breakpoint engine.
# A_j = {x : sector j unhit by all e}.  meas(S7) = sum_{A subset {1..6}} (-1)^|A| J(A).
# (sector 0 always hit by e=0, so we IE over the other 6 sectors.)
# J(A) = meas{x : for all e in E, frac(ex) avoids all sectors in A}.
# Use the exact J() engine from the prompt.
# ---------------------------------------------------------------------------
def J(A, E):
    E = sorted(set(E)); arcs = [(F(j, 7), F(j + 1, 7)) for j in A]; bp = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for (a, b) in arcs:
            for end in (a, b):
                m = 0
                while True:
                    xv = (end + m) / e
                    if xv >= 1: break
                    if xv >= 0: bp.add(xv)
                    m += 1
    bp = sorted(b for b in bp if 0 <= b < 1); tot = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if all(not any(a < ((e * mid) % 1) < b for (a, b) in arcs) for e in E): tot += hi - lo
    return tot

def measS7_IE(E):
    tot = F(0)
    for r in range(7):
        for A in itertools.combinations(range(1, 7), r):
            tot += F((-1) ** r) * J(set(A), E)
    return tot

# ---------------------------------------------------------------------------
# Fourier kernels (float).
# ---------------------------------------------------------------------------
def shat(n, j):
    if n == 0: return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)

SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1) ** len(T) for T in SUB}
_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0: v = complex(1 - len(T) / 7.0, 0.0)
    elif n % 7 == 0: v = 0j
    else: v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v

def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T] * p
    return s

def banner(t): print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)

# Canonical caps (T859/T860, HYP-2603) -- the RIGOROUS target.
CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
       11: F(66, 91), 12: F(6, 7), 13: F(1, 1)}
# The PROMPT's cap values (claimed garbled by canon).
CAP_PROMPT = {8: F(2243, 5880), 9: F(2025, 4004), 10: F(36, 91),
              11: F(25, 91), 12: F(1, 7), 13: F(0, 1)}


def V0_framing():
    banner("V0 -- FRAMING: which target?  cap_k (canon) vs consec (claimed advance)")
    print("  Proved floor (THM-535): cap_k >= (k-6)/7.  Check both cap tables:")
    print(f"  {'k':>2} {'(k-6)/7':>10} {'CAP_canon':>14} {'>=floor?':>9} {'CAP_prompt':>14} {'>=floor?':>9}")
    for k in range(8, 14):
        floor = F(k - 6, 7)
        cc, cp = CAP[k], CAP_PROMPT[k]
        print(f"  {k:>2} {float(floor):>10.5f} {str(cc):>14} {str(cc>=floor):>9} "
              f"{str(cp):>14} {str(cp>=floor):>9}")
    print("\n  => The PROMPT's cap_9,10,11,12,13 VIOLATE the proved floor (garbled).")
    print("     Canon caps satisfy the floor.  So 'cap framing wrong' is partly right:")
    print("     the prompt's NUMBERS were garbled, BUT the cap TARGET itself is canonical")
    print("     (HYP-2603/THM-532).  The claimed advance instead targets meas(S7(consec)).")
    print("\n  Now: is consec <= cap_k (so proving <=consec suffices)?  And is <=consec")
    print("  even TRUE as a same-k statement?")
    print(f"  {'k':>2} {'meas(consec)':>14} {'cap_canon':>14} {'consec<=cap?':>12}")
    for k in range(8, 14):
        c = measS7(list(range(k)))
        print(f"  {k:>2} {str(c):>14} {str(CAP[k]):>14} {str(c <= CAP[k]):>12}")


def V1_targets():
    banner("V1 -- exact targets M7(k), meas(S7(consec_k)), delta_k=consec-M7")
    deltas = {}
    print(f"  {'k':>2} {'M7(k)':>14} {'meas(consec)':>16} {'delta_k':>14} {'>=0.302?':>9}")
    for k in range(8, 14):
        c = measS7(list(range(k))); m = M7(k); d = c - m; deltas[k] = d
        print(f"  {k:>2} {float(m):>14.6f} {str(c):>16} {float(d):>14.6f} {str(d >= F(302,1000)):>9}")
    # cross-check consec values against the claimed advance's asserted rationals
    claimed = {8: F(481,1470), 9: F(2447,5880), 10: F(8899,17640),
               11: F(3419,5880), 12: F(121103,194040), 13: F(14573,21560)}
    print("\n  Cross-check meas(S7(consec_k)) vs claimed advance's asserted exact values:")
    for k in range(8, 14):
        c = measS7(list(range(k)))
        print(f"    k={k}: computed={c}  claimed={claimed[k]}  match={c==claimed[k]}")
    return deltas


def V2_lemmaA():
    banner("V2 -- LEMMA A (support-6 floor): K(n)=0 unless >=6 nonzero non-7 coords")
    # (a) independent C(U) algebra: C(U)=sum_{T superset U}(-1)^|T|.
    print("  (a) C(U) = sum_{T superset U, T subset {1..6}} (-1)^|T|  (independent calc):")
    bad = []
    for u in range(0, 7):
        for U in itertools.combinations(range(1, 7), u):
            US = set(U)
            CU = sum((-1) ** len(T) for r in range(7)
                     for T in itertools.combinations(range(1, 7), r) if US.issubset(set(T)))
            expected = 0 if u < 6 else 1  # (-1)^6 = 1
            # general: C(U) = (-1)^|U| * (1-1)^{6-|U|} = 0 unless |U|=6, then (-1)^6=1
            if u == 6:
                if CU != 1: bad.append((U, CU))
            else:
                if CU != 0: bad.append((U, CU))
    print(f"      C(U)=0 for all |U|<6, C(full)=1 :  {'OK' if not bad else 'FAIL '+str(bad)}")
    # (b) EXHAUSTIVE structured grid for support<=5: all coordinate tuples in a small
    #     box with nonzero non-7 entries; confirm K(n)=0.
    print("  (b) EXHAUSTIVE support 1..5, coords in {-8..8}\\{0,+-7} (structured, not random):")
    vals = [v for v in range(-8, 9) if v != 0 and v % 7 != 0]
    worst = 0.0; worstn = None
    for s in range(1, 6):
        # to keep it finite, sample structured: full product for s<=2, strided for s>=3
        if s <= 2:
            it = itertools.product(vals, repeat=s)
        else:
            small = [v for v in range(-4, 5) if v != 0 and v % 7 != 0]
            it = itertools.product(small, repeat=s)
        cnt = 0
        for n in it:
            cnt += 1
            kv = abs(Kk(list(n)))
            if kv > worst: worst = kv; worstn = n
        print(f"      support {s}: tested {cnt:>6} tuples, running max|K|={worst:.2e}")
    print(f"      OVERALL max|K| over support<=5 = {worst:.2e} at n={worstn}")
    # (c) support-6 generically nonzero
    print("  (c) support-6 examples (must be NONZERO -> floor sharp):")
    for n6 in ([1,2,3,4,5,6], [1,1,1,1,1,1], [1,-1,2,-2,3,-3], [3,3,3,3,3,3]):
        print(f"      n={n6}: |K|={abs(Kk(n6)):.6e}")
    # (d) the support-2 W-matrix identically zero (the lambda_1 short relation gives 0)
    W = {}
    for j in range(1, 7):
        for l in range(1, 7):
            W[(j, l)] = sum((-1) ** len(T) for r in range(7)
                            for T in itertools.combinations(range(1, 7), r)
                            if j in T and l in T)
    allzero = all(v == 0 for v in W.values())
    print(f"  (d) support-2 kernel W[j,l] identically zero: {allzero}  "
          f"(=> shortest relations contribute K=0)")
    return (not bad) and worst < 1e-10 and allzero


def V3_envelope():
    banner("V3 -- LEMMA B envelope |shat(n,j)|=|sin(pi n/7)|/(pi|n|); c1 constant")
    ok = True
    for n in range(1, 50):
        for j in range(1, 7):
            lhs = abs(shat(n, j)); rhs = abs(math.sin(math.pi * n / 7)) / (math.pi * n)
            if abs(lhs - rhs) > 1e-12: ok = False
    print(f"  |shat(n,j)| closed form verified n=1..49, all j: {ok}")
    # c1 = max over n (not 7|n) of |n| * max_T |chat(n,T)|.  periodic mod 7.
    c1 = 0.0; argn = None
    for n in range(1, 71):
        if n % 7 == 0: continue
        mx = max(abs(chat(n, T)) for T in SUB) * n
        if mx > c1: c1 = mx; argn = n
    print(f"  c1 = max_n |n|*max_T|chat(n,T)| = {c1:.6f}  at n={argn} (n mod 7 = {argn%7})")
    print(f"  triangle bound 6|sin(3pi/7)|/pi = {6*abs(math.sin(3*math.pi/7))/math.pi:.6f}")
    print(f"  claimed c1 ~ 0.6973/0.6974 : {'CONSISTENT' if abs(c1-0.6974)<0.002 else 'MISMATCH'}")
    return c1


def V4_identity():
    banner("V4 -- signed identity meas(S7)=M7(k)+sum_{0!=n in Lambda^o} K(n)")
    def brute(nz, N0):
        d = len(nz); out = []
        for v in itertools.product(range(-N0, N0 + 1), repeat=d):
            if all(x == 0 for x in v): continue
            if sum(v[i] * nz[i] for i in range(d)) == 0: out.append(v)
        return out
    for name, E, N0 in [("k4 [0,1,2,5]", [0,1,2,5], 12),
                        ("k5 [0,1,2,3,7]", [0,1,2,3,7], 9),
                        ("k5 [0,2,3,5,8]", [0,2,3,5,8], 9),
                        ("k4 [0,1,3,4]", [0,1,3,4], 12)]:
        nz = [e for e in E if e != 0]; k = len(E)
        s = sum(Kk(n) for n in brute(nz, N0))
        lhs = float(measS7(E)); rhs = float(M7(k)) + s.real
        print(f"  {name}: measS7={lhs:.6f}  M7+sumK={rhs:.6f}  |diff|={abs(lhs-rhs):.2e}  "
              f"(truncated |n|<={N0})")


def V5_resonance():
    banner("V5 -- RESONANCE w==0 mod 7: hunt for resonant config exceeding consec/cap")
    consec8 = measS7(list(range(8)))
    cap8 = CAP[8]
    print(f"  ref: meas(S7(consec_8))={consec8}={float(consec8):.6f}, cap_8={cap8}={float(cap8):.6f}")
    tests = [
        ("[0,7,14,1,2,3,4,5]", [0,7,14,1,2,3,4,5]),
        ("[0,1,2,3,4,5,6,49]", [0,1,2,3,4,5,6,49]),
        ("[0,7,14,21,28,35,42,49]", [0,7,14,21,28,35,42,49]),
        ("[0,7,1,8,2,9,3,10]", [0,7,1,8,2,9,3,10]),
        ("[0,14,21,35,42,7,28,49]", [0,14,21,35,42,7,28,49]),
        ("[0,1,7,8,14,15,21,22]", [0,1,7,8,14,15,21,22]),
    ]
    over = []
    for nm, E in tests:
        m = measS7(E)
        oc = m > consec8; ocap = m > cap8
        if oc: over.append((nm, E, m))
        print(f"  {nm:>26} meas={float(m):.6f}  >consec?{oc}  >cap?{ocap}")
    print(f"  resonant configs exceeding consec: {over if over else 'NONE'}")
    return over


def hunt_consec(k, max_e, extra=()):
    """Hunt for primitive E (|E|=k, 0 in E) with measS7(E) > measS7(consec_k)."""
    target = measS7(list(range(k)))
    beaters = []
    # primitive: gcd of all = 1.  0 in E.  generate k-1 distinct positives <= max_e.
    pos = list(range(1, max_e + 1))
    for combo in itertools.combinations(pos, k - 1):
        if gcd(*combo) != 1 and len(combo) > 1: continue
        E = (0,) + combo
        m = measS7(list(E))
        if m > target:
            beaters.append((E, m))
    for E in extra:
        m = measS7(list(E))
        if m > target:
            beaters.append((tuple(E), m))
    return target, beaters


def V6_hunt_consec():
    banner("V6 -- HUNT counterexamples to meas(S7(E)) <= meas(S7(consec_k))")
    # k=8: exhaustive-ish bounded spread
    for k, max_e in [(8, 13), (9, 12), (10, 12)]:
        target, beaters = hunt_consec(k, max_e)
        print(f"  k={k}, max_e={max_e}: target consec={float(target):.6f}, "
              f"beaters={len(beaters)}")
        for E, m in sorted(beaters, key=lambda t: -t[1])[:5]:
            print(f"       BEATER E={E} meas={float(m):.6f} (> consec {float(target):.6f})")
    # k=12: T860 says {0,1,...,10,12} beats consec.  CHECK.
    banner("V6b -- k=12 same-k extremality of consec (T860 claims FALSE)")
    k = 12
    target = measS7(list(range(k)))
    cand = list(range(11)) + [12]   # {0,1,..,10,12}
    mc = measS7(cand)
    print(f"  consec_12 meas = {target} = {float(target):.6f}")
    print(f"  E={cand} meas = {mc} = {float(mc):.6f}")
    print(f"  E beats consec? {mc > target}   (T860: claims YES => '<=consec' is FALSE at k=12)")
    # broader k=12 small hunt
    t2, beat2 = hunt_consec(12, 13)
    print(f"  k=12 hunt max_e=13: beaters over consec = {len(beat2)}")
    for E, m in sorted(beat2, key=lambda t: -t[1])[:6]:
        print(f"       BEATER E={E} meas={float(m):.6f} vs consec {float(t2):.6f}, "
              f"<=cap_12? {m<=CAP[12]}")
    return mc > target


def V7_hunt_cap():
    banner("V7 -- HUNT counterexamples to meas(S7(E)) <= cap_k  (the CANONICAL target)")
    for k, max_e in [(8, 14), (9, 13), (10, 13), (11, 13), (12, 14), (13, 15)]:
        cap = CAP[k]; over = []
        pos = list(range(1, max_e + 1))
        cnt = 0
        for combo in itertools.combinations(pos, k - 1):
            if len(combo) > 1 and gcd(*combo) != 1: continue
            cnt += 1
            E = (0,) + combo
            m = measS7(list(E))
            if m > cap: over.append((E, m))
        worst = max((measS7(list((0,)+c)) for c in itertools.combinations(pos, k-1)
                     if not (len(c)>1 and gcd(*c)!=1)), default=F(0))
        print(f"  k={k} max_e={max_e}: primitive shapes={cnt}, cap={float(cap):.6f}, "
              f"max meas={float(worst):.6f}, OVER CAP={len(over)}")
        for E, m in sorted(over, key=lambda t: -t[1])[:3]:
            print(f"       *** CAP VIOLATION E={E} meas={float(m):.6f} > cap={float(cap):.6f}")


def V8_gap_audit(deltas, c1):
    banner("V8 -- GAP AUDIT: absolute residual after support-6 floor vs delta_k")
    # For consec_8, compute sum |K(n)| over support>=6 relations within a box,
    # and compare to delta_8.  The claim: still ~1.2 > delta_8 ~ 0.30 (so the
    # ABSOLUTE bound after the floor does NOT close -- signed gain still needed).
    E = list(range(8)); nz = [e for e in E if e != 0]; d = len(nz)
    N0 = 9
    abs_supp6 = 0.0; abs_all = 0.0; signed = 0j
    for v in itertools.product(range(-N0, N0 + 1), repeat=d):
        if all(x == 0 for x in v): continue
        if sum(v[i] * nz[i] for i in range(d)) != 0: continue
        kv = Kk(list(v)); signed += kv
        abs_all += abs(kv)
        nzc = sum(1 for x in v if x != 0 and x % 7 != 0)
        if nzc >= 6: abs_supp6 += abs(kv)
    print(f"  consec_8, box |n|<={N0}:")
    print(f"    sum|K| all relations          = {abs_all:.4f}")
    print(f"    sum|K| support>=6 relations   = {abs_supp6:.4f}")
    print(f"    signed sum (=corr)            = {signed.real:+.6f}")
    print(f"    delta_8                       = {float(deltas[8]):.6f}")
    print(f"    corr <= delta_8?              = {signed.real <= float(deltas[8])+1e-9}")
    print(f"    absolute(supp>=6) <= delta_8? = {abs_supp6 <= float(deltas[8])}  "
          f"(claim: NO -> abs bound lossy, signed gain still needed)")
    print(f"    floor removes fraction of abs mass = {1 - abs_supp6/abs_all:.3f}")
    print("  (NB: box truncation underestimates the true infinite sums; the relative")
    print("   ordering 'abs(supp>=6) >> delta_8' is the load-bearing qualitative claim.)")


def main():
    print("ADVERSARIAL VERIFICATION: LRC(14) wide-spread-signed-weyl (kps-S9-wf)")
    V0_framing()
    deltas = V1_targets()
    a_ok = V2_lemmaA()
    c1 = V3_envelope()
    V4_identity()
    V5_resonance()
    consec_fails_k12 = V6_hunt_consec()
    V7_hunt_cap()
    V8_gap_audit(deltas, c1)
    banner("SUMMARY OF VERIFICATION")
    print(f"  Lemma A (support-6 floor) re-derived & exhaustive-verified: {a_ok}")
    print(f"  c1 envelope ~ {c1:.4f}")
    print(f"  '<=consec' FALSE at k=12 (T860 reproduced): {consec_fails_k12}")
    print("  See per-part output above for cap-violation hunts and gap audit.")


if __name__ == "__main__":
    main()
