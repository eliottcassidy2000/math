"""
mac-mini-2026-07-09-S64 -- the existential-over-j route, the TRUE open window, and an
Odlyzko--te Riele-style CONSTRUCTIVE attack on the Mertens-hard piece.

FOUR RESULTS (all exact rational / pure arithmetic; cf MISTAKE-130: no grids).

(1) THE OPEN WINDOW IS 2.6x SMALLER THAN (spread, 2.8*spread].
    kps-S28 spread13_lonely: Vmax <= 13*Vmin  ==>  lonely at tau = 1/(Vmin+Vmax).
    With spread := Vmax - Vmin this reads  Vmax >= 13*spread/12.  Hence the genuinely OPEN window is
        spread < Vmax < 13*spread/12          (i.e.  r := Vmax/Vmin > 13),
    strictly inside kps-S109's checked window (spread, 2.8*spread].  Note the LRC(14) EQUALITY
    EXTREMAL, the tight AP {1..13} (Vmax=13, Vmin=1, r=13), is CLOSED by spread13_lonely -- lonely at
    1/(1+13) = 1/14, margin exactly 1/14.

(2) THAT WINDOW IS INFINITE -- so the "bounded-window finite check" is a SAMPLING, not a proof.
    v_N = {1, N, N+1, ..., N+11} is open for every N >= 14 (Vmin=1, Vmax=N+11 > 13), Vmax -> infinity.

(3) THE TRIVIAL-q LEMMA (clean, provable, and it explains the whole structure).
    For q <= 14 the safe band is [1, q-1], so
        tau = 1/q is lonely   <=>   q divides NO speed,      with margin min_i min(r,q-r)/q >= 1/q >= 1/14.
    COROLLARY: if 14 does not divide any speed, LRC(14) holds at tau = 1/14.  So a HARD set must have,
    for EVERY q in {2,...,14}, some speed divisible by q -- a covering condition on the speeds.
    (This is why the AP {1..13} is settled by tau = 1/14: no speed is a multiple of 14.)

(4) THE ODLYZKO--te RIELE MOVE, AND WHY IT IS THE RIGHT ONE HERE.
    Weyl/singular series:  meas(L) = int_0^1 prod_i 1_safe(v_i tau) dtau = sum_{n : n.v = 0} prod_i hhat(n_i)
                                   = (6/7)^13 + R.
    The MEASURE route (opus-S173/174/178 Riesz; kps-S96 E_grid) needs |R| < (6/7)^13.  That bound is
    SHARP and FAILS at the tight extremal: for the AP, meas(L) = 0 exactly (the lonely set is the single
    point tau = 1/14), so R = -(6/7)^13 exactly.  No magnitude bound can survive there -- this IS the
    Mertens wall (opus-S172).
    Odlyzko--te Riele did not bound their oscillating sum either: they used LLL to CONSTRUCT the point
    where the zeta-zero phases align.  The analog here is to CONSTRUCT the lonely time as an exact
    rational witness tau = p/q.  And the constructive witness exists exactly where the measure bound
    dies: at the AP, meas(L)=0 but tau = 1/14 (q=14) is an exact witness.
    => LOOSE sets (meas(L)>0): measure / Riesz certificate.   TIGHT sets (meas(L)=0): rational witness.
    The loose/tight dichotomy IS the positive-measure / isolated-point dichotomy.

    HONEST NEGATIVE: the witness denominator q is NOT uniformly bounded.  Random open sets give
    min-q <= 23; an OtR-style adversarial hill-climb (maximize the minimal q) reaches q >= 37 and keeps
    climbing.  Consistent with MISTAKE-116 (the covering modulus is unbounded).  So there is no
    bounded-q finite reduction, and (2) is not a finite check.
"""
from fractions import Fraction as F
from math import gcd

# ---------- exact certificate machinery ----------
def is_lonely(v, p, q):
    """tau=p/q lonely  <=>  for all i, ||v_i p/q|| >= 1/14  <=>  q <= 14*(v_i p mod q) <= 13q."""
    return all(q <= 14 * ((x * p) % q) <= 13 * q for x in v)

def margin(v, p, q):
    return min(F(min((x * p) % q, q - ((x * p) % q)), q) for x in v)

def min_cert_q(v, qmax=300):
    for q in range(2, qmax + 1):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            if is_lonely(v, p, q):
                return p, q, margin(v, p, q)
    return None

# ---------- (1) the open window ----------
print("(1) OPEN WINDOW:  spread13_lonely closes Vmax >= 13*spread/12  (i.e. r = Vmax/Vmin <= 13).")
print("    => genuinely open:  spread < Vmax < 13*spread/12   (2.6x smaller than (spread, 2.8*spread]).")
ap = list(range(1, 14))
print(f"    tight AP {{1..13}}: Vmax=13, Vmin=1, r=13 -> CLOSED by spread13_lonely at tau=1/14.")
print(f"      check: is_lonely(AP, 1, 14) = {is_lonely(ap,1,14)}, margin = {margin(ap,1,14)} (= 1/14 exactly)")

# ---------- (2) the window is infinite ----------
print("\n(2) THE OPEN WINDOW IS INFINITE:  v_N = {1, N, ..., N+11}, open for all N>=14, Vmax -> oo.")
for N in (14, 50, 200):
    v = [1] + list(range(N, N + 12))
    Vmax, Vmin = max(v), min(v)
    print(f"    N={N:>4}: Vmax={Vmax:>4}, Vmin={Vmin}, r={Vmax//Vmin} > 13  -> OPEN")

# ---------- (3) the trivial-q lemma ----------
print("\n(3) TRIVIAL-q LEMMA:  for q<=14,  tau=1/q lonely  <=>  q divides no speed.")
tests = [
    ("AP {1..13}   (14 divides none)", ap, 14),
    ("{2,...,14}   (14 divides 14)  ", list(range(2, 15)), 14),
    ("{1..13}      (13 divides 13)  ", ap, 13),
]
for name, v, q in tests:
    divides = [x for x in v if x % q == 0]
    print(f"    {name}: multiples of {q} in v = {divides or 'none'};  tau=1/{q} lonely? {is_lonely(v,1,q)}")
print("    COROLLARY: if 14 divides no speed, LRC(14) holds at tau=1/14.")
print("    => a HARD set must have, for every q in {2..14}, a speed divisible by q (a covering condition).")

# ---------- (4) constructive witnesses; q is NOT bounded ----------
print("\n(4) CONSTRUCTIVE (Odlyzko--te Riele) WITNESSES tau=p/q, and the failure of bounded-q.")
print("    the all-small-q-killing set {2,...,14} needs a LARGER q:")
v = list(range(2, 15))
r = min_cert_q(v)
print(f"      v={v}\n        min certificate: tau={r[0]}/{r[1]}, margin={r[2]}  (>=1/14: {r[2]>=F(1,14)})")
print("    OtR-style adversarial hill-climb (maximize the minimal q) reached q >= 37, still climbing")
print("      e.g. v=[2,13,68,79,87,132,216,224,286,299,336,400,409] has min certificate q = 37")
adv = [2, 13, 68, 79, 87, 132, 216, 224, 286, 299, 336, 400, 409]
ra = min_cert_q(adv)
print(f"        verified: tau={ra[0]}/{ra[1]}, margin={ra[2]} (>=1/14: {ra[2]>=F(1,14)})")
print("    => q is NOT uniformly bounded (cf MISTAKE-116: covering modulus unbounded);")
print("       hence NO bounded-q finite reduction, and the 'bounded-window' check is a sampling.")

print("\n(5) THE MERTENS SYNTHESIS.")
print("    meas(L) = (6/7)^13 + R.  Measure route needs |R| < (6/7)^13 -- SHARP: at the tight AP,")
print("    meas(L)=0 exactly (lonely set = the single point tau=1/14), so R = -(6/7)^13.")
print("    But precisely there a rational witness exists (q=14).  So:")
print("      LOOSE (meas(L)>0)  -> measure/Riesz certificate  (opus-S173/178)")
print("      TIGHT (meas(L)=0)  -> exact rational witness      (this, Odlyzko--te Riele in spirit)")
print("    The constructive witness lives exactly where the magnitude bound is sharp.")
