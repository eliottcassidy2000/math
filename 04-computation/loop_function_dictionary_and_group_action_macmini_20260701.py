#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S79 -- THE LOOP-FUNCTION DICTIONARY (expanded) + THE GROUP-ACTION REFORMULATION.

Builds on opus-S13's circle-map dictionary (AGL(1)+PSL2Z+Szego) and kps-S7's Verblunsky/Ramanujan atoms.
Two deliverables:

  (I)  DICTIONARY of maps between points on the loop R/Z, organized by the GROUP each generates, with the
       group relations verified. ~24 entries across 6 families (isometry, endomorphism, Mobius/Verblunsky,
       flow/dynamics, configuration maps, cocycles/measures).

  (II) THE CLEAN NEXT STEP (synthesis): the LRC WITNESS SEARCH IS A GROUP ACTION. At modulus q the dilation
       group G = (Z/q)* acts diagonally on the residue configuration c=(v_1,...,v_13) mod q. A rational time
       a/q is lonely (>= r) iff a.c lands in the SAFE BOX B_r = [rq, q-rq]^13. So:
         M(S) >= r  <=>  the orbit G.c meets the safe box B_r   (for SOME modulus q).
         covering-min lower bound  <=>  every covering config's dilation-orbit meets B_r.
       This unifies: S77 safe-band frame (the box) + kps units/Ramanujan (G's regular orbit = the atoms;
       Ramanujan sums = G's Fourier) + S78 MSS finiteness (q, |G| bounded) + opus AGL/PSL2 (G in the
       dictionary) + the Dedekind sum (G's cocycle) + Verblunsky/Blaschke (the continuous disk completion).
       The character-sum count |G.c cap B_r| via Ramanujan sums is the analytic engine (kps-S7).
"""
import numpy as np
from math import gcd, comb
from fractions import Fraction as F

def frac(x): return x - np.floor(x)

print("="*94)
print("(I) THE LOOP-FUNCTION DICTIONARY -- maps R/Z -> R/Z (and config maps), by group")
print("="*94)
entries = [
 ("ISOMETRY O(2)",[
   ("rot_c        x -> x+c",                 "rotation = observer/origin shift (inhomogeneous LRC)"),
   ("refl (iota)  x -> -x",                  "antipode/complement; the killer -1 symmetry (t<->1-t)"),
 ]),
 ("ENDOMORPHISM / AFFINE AGL(1)",[
   ("dil_v        x -> v x   (v in Z)",      "a RUNNER of speed v; the x |-> vx map"),
   ("unit-dil_a   x -> a x   (a in (Z/q)*)", "INVERTIBLE dilation = the WITNESS SEARCH generator"),
   ("aff_{v,a}    x -> v x + a",             "the ax+b group; runner v seen from moving observer a"),
   ("times-n      x -> n x   (endo)",        "Bernoulli/doubling endomorphism (expanding)"),
 ]),
 ("MOBIUS / VERBLUNSKY (disk aut = PSL2R boundary)",[
   ("blaschke_A   z -> (z-A)/(1-Abar z)",    "Verblunsky/OPUC atom-pusher (|A|->1 = atomic = lonely)"),
   ("cayley       disk <-> half-plane",      "the standard disk/H transfer"),
   ("mobius       x -> (ax+b)/(cx+d)",       "fractional-linear; PSL2 on the boundary circle"),
   ("gauss        x -> {1/x}",               "continued-fraction shift = PSL2Z generator (S<->T)"),
   ("mediant      (a/b,c/d)->(a+c)/(b+d)",   "Farey/Stern-Brocot; builds t*=[0;n-1,n]"),
 ]),
 ("FLOW / DYNAMICS",[
   ("flow_v,t     x -> x + v t",             "the runner FLOW; LRC = avoid arc(-r,r) for all v at one t"),
   ("3-gap return interval-exchange",        "three-distance first-return (Steinhaus); the S73 gaps"),
   ("winding      rotation number",          "the a in a/q = winding of the witness"),
 ]),
 ("CONFIGURATION MAPS (multiset -> multiset)",[
   ("scale_c      {v_i} -> {c v_i}",         "dilate all runners (scale-invariance THM-522)"),
   ("insert/del   {v_i} -> {v_i} +/- {w}",   "vertex insertion/removal (Mode A; remove-large)"),
   ("conv/Mink    mu*nu / A+B",              "convolution of point measures / Minkowski sum on loop"),
   ("W_N          Atkin-Lehner involution",  "modular involution on X_0(N) (apex X_0(14))"),
   ("ramanujan_hat c -> (c_q(k))_k",         "the config's Fourier = Ramanujan sums (kps-S7)"),
 ]),
 ("COCYCLES / MEASURES (functions ON the group)",[
   ("dedekind s(a,q)=sum((ka/q))((k/q))",    "the PSL2Z cocycle = the LRC MARGIN (HYP-3768)"),
   ("verblunsky  alpha_j (Szego recursion)", "reflection tower; |alpha_j|=1/(n-1-j) for AP (opus-S13)"),
   ("saw ((x))=x-floor-1/2",                 "iota-ODD sawtooth (Dedekind summand)"),
   ("dist ||x||",                            "iota-EVEN distance-to-integer (the loneliness)"),
 ]),
]
k = 0
for fam, items in entries:
    print(f"\n  [{fam}]")
    for name, note in items:
        k += 1
        print(f"    {k:2d}. {name:34s} {note}")
print(f"\n  => {k} loop-functions. GROUP STRUCTURE: isometry O(2) + dilation = AFFINE AGL(1,R/Z);")
print(f"     Mobius/Gauss/mediant generate PSL2(Z) (CF descent); Blaschke = disk-aut PSL2(R) (Verblunsky);")
print(f"     the arithmetic apex = Gamma_0(14) / PSL2(7) (Klein quartic). Runners = the DILATION subgroup.")

# ---- verify the group relations (they operate group-like) ----
xs = np.linspace(0, 1, 13, endpoint=False)
rot = lambda a: (lambda x: frac(x+a)); dil = lambda v: (lambda x: frac(v*x)); refl = lambda x: frac(-x)
aff = lambda v,a: (lambda x: frac(v*x+a))
R = {
 "dil compose  dil_a dil_b = dil_ab": np.allclose(dil(3)(dil(5)(xs)), dil(15)(xs)),
 "affine law   a_{v,a}a_{w,b}=a_{vw,vb+a}": np.allclose(aff(3,0.1)(aff(5,0.2)(xs)), aff(15, frac(3*0.2+0.1))(xs)),
 "refl-conj    refl dil refl = dil": np.allclose(refl(dil(5)(refl(xs))), dil(5)(xs)),
 "semidirect   dil rot = rot' dil": np.allclose(dil(5)(rot(0.1)(xs)), rot(frac(5*0.1))(dil(5)(xs))),
}
print("\n  GROUP RELATIONS verified:", "; ".join(f"{k2}={v}" for k2,v in R.items()))

print()
print("="*94)
print("(II) THE CLEAN NEXT STEP -- THE LRC WITNESS SEARCH IS A GROUP ACTION (dilation orbit hits safe box)")
print("="*94)
n = 14; r = F(1, n)              # LRC threshold; covering-min uses r=n/Phi6
def normq(x, q): return min(x % q, q - (x % q))
def orbit_hits_box(S, q, r):
    """does the (Z/q)* dilation orbit of S mod q meet the safe box (all coords >= r q)?"""
    band = r * q                 # half-width; safe iff normq >= band
    units = [a for a in range(1, q) if gcd(a, q) == 1]
    best = 0; best_a = None
    for a in units:
        m = min(normq(a * (v % q), q) for v in S)
        if m > best: best = m; best_a = a
    return best, best_a, F(int(best), q)

print("Safe box B_r = {c : ||c_i||_q >= rq for all i}. Relative volume = (1-2r)^13:")
for rr, lbl in [(F(1,14), "LRC 1/14"), (F(14,183), "covering-min 14/183")]:
    vol = (1 - 2*float(rr))**13
    print(f"  r={lbl}: (1-2r)^13 = {vol:.4f}  (the orbit must hit a box of this relative volume)")

print("\nConstruction {1..12,182} -- the (Z/q)* orbit reaches the box at the RIGHT modulus q:")
S = list(range(1, 13)) + [182]
for q in [14, 15, 91, 183]:
    best, a, Mq = orbit_hits_box(S, q, F(1,14))
    hit = "HITS box(1/14)" if Mq >= F(1,14) else "misses"
    print(f"  q={q:4d}: best dilation a={str(a):>4}, min-dist={best}, M-at-q={Mq}={float(Mq):.4f}  {hit}")

print("\nThe orbit view of a few covering 13-sets at q=183 (does the dilation orbit clear 1/14?):")
def covers(S,n=14): return all(any(v%qq==0 for v in S) for qq in range(2,n+1))
tests = {
  "construction": list(range(1,13))+[182],
  "AP {1..13} (non-cov)": list(range(1,14)),
  "GW {1..11,13,24}": [1,2,3,4,5,6,7,8,9,10,11,13,24],
}
for name, S in tests.items():
    # best over a RANGE of moduli (the witness modulus is adaptive)
    bestM = F(0); bq=None
    for q in range(n, 400):
        _,_,Mq = orbit_hits_box(S, q, F(1,14))
        if Mq > bestM: bestM = Mq; bq = q
    print(f"  {name:24s} covering={covers(S)}  best over q<400: M={bestM}={float(bestM):.4f} at q={bq}  {'>=1/14' if bestM>=F(1,14) else '<1/14'}")

print()
print("SYNTHESIS: M(S) = max_q [ (1/q) max_{a in (Z/q)*} min_i ||a v_i||_q ] = max_q (1/q)(orbit-best at q).")
print("  => covering-min lower bound = EVERY covering config's (Z/q)* dilation-orbit meets the safe box B_{1/14}")
print("     at some modulus q<=91^12 (MSS, S78). The atoms (kps) = the regular orbit = units; the Ramanujan")
print("     sums c_q(k) = the group Fourier => |orbit cap box| is a CHARACTER SUM (the analytic engine).")
print("DONE.")
