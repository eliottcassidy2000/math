#!/usr/bin/env python3
"""
lrc14_far_element_J_recursion_klein_S271.py
===========================================
klein-2026-07-12-S271  (owner: push consec-maximizes / three-gap toward proof)

THE CRUX (Route A [A], THM-711/716): consec_k minimizes J = E[N(7-N)] over k-sets, where
N = N_empty = #{empty sectors among the 7 arcs [s/7,(s+1)/7)} of the phases {frac(e_i x)}.
Operative requirement: J >= 432/91 = 4.7473 (consec_9 = 4465/882 = 5.0624, margin +0.315).
Compact (diam<=18) exhaustively verified (kps cont.32). OPEN: the TAIL (wide sets) rigorous.

THE PUSH: a FAR-ELEMENT RECURSION. Add a decorrelated element w to a k-set E: w lands in a
uniform-random sector independent of E, emptying one previously-empty sector w.p. N/7. The
algebra (total-variance) gives, in the far/decorrelated limit,

      J_{k+1}(E u {far w}) = (5/7) J_k(E) + (6/7) mu_k(E)          [CLAIM]

with mu' = (6/7) mu. If exact, this (a) quantifies "far elements raise J", (b) reduces the
k=9 tail to the k=8 base + a mu bound -- a rigorous step (the finite-far error is THM-687/688).

This script:
 (1) computes J(E), mu(E) on a fine grid (approx, to check the 5/7, 6/7 coefficients);
 (2) verifies the recursion J(E u {W}) ~ (5/7)J(E)+(6/7)mu(E) as W -> large (decorrelated);
 (3) checks the tail consequence: does the k=8 base + mu give J_9 >= 432/91 for wide 9-sets?
"""
import math
from fractions import Fraction

def phases_empty_count(E, x):
    """N_empty = # of the 7 sectors [s/7,(s+1)/7) containing NO phase frac(e*x)."""
    occ=[False]*7
    for e in E:
        f=(e*x)%1.0
        occ[int(f*7)%7]=True
    return 7-sum(occ)

def J_mu(E, Ngrid=200000):
    """E[N(7-N)] and E[N] over x on a fine grid (approx)."""
    sN=0.0; sJ=0.0
    for k in range(1,Ngrid):
        x=k/Ngrid
        N=phases_empty_count(E,x)
        sN+=N; sJ+=N*(7-N)
    return sJ/(Ngrid-1), sN/(Ngrid-1)

thr=432/91; consec9=4465/882; consec8=291/49
print("="*72)
print("(1) baseline J, mu (fine grid)")
print("="*72)
for nm,E in [("consec8 {1..8}",list(range(1,9))),("consec9 {1..9}",list(range(1,10)))]:
    J,mu=J_mu(E)
    print(f"  {nm:16s} J={J:.4f}  mu={mu:.4f}   (canon: consec9=4465/882={consec9:.4f}, consec8=291/49={consec8:.4f})")
print(f"  threshold 432/91 = {thr:.4f}")

print()
print("="*72)
print("(2) FAR-ELEMENT RECURSION: J(E u {W}) vs (5/7)J(E)+(6/7)mu(E) as W grows")
print("="*72)
for nm,E in [("consec8 {1..8}",list(range(1,9))),("2AP {2,4,..,16}",list(range(2,17,2)))]:
    J8,mu8=J_mu(E)
    pred=(5/7)*J8+(6/7)*mu8
    print(f"  E={nm}: J(E)={J8:.4f} mu(E)={mu8:.4f} => predicted J(E u far)=(5/7)J+(6/7)mu={pred:.4f}")
    for W in [97,997,9973,99991]:
        J9,mu9=J_mu(E+[W])
        print(f"     W={W:6d}: J(EuW)={J9:.4f} (pred {pred:.4f}, err {J9-pred:+.4f})  mu(EuW)={mu9:.4f} (pred {6/7*mu8:.4f})")

print()
print("="*72)
print("(3) TAIL CONSEQUENCE: does k=8 base + mu give J_9 >= 432/91 for wide 9-sets?")
print("="*72)
print("  A wide 9-set = compact k'-cluster C + far element(s). In the far limit, adding 1 far w to")
print("  a compact 8-set C gives J_9 -> (5/7)J_8(C) + (6/7)mu_8(C).")
print("  If the k=8 base holds (J_8(C) >= 291/49 = 5.9388) and mu_8(C) >= mu_min, then:")
# what mu_min is needed? J_9 >= thr <=> (5/7)J8 + (6/7)mu8 >= 432/91
# with J8 >= 291/49: (6/7)mu8 >= 432/91 - (5/7)(291/49) = 432/91 - 1455/343
need = thr - (5/7)*consec8
mu_need = need*7/6
print(f"     need (6/7)mu8 >= 432/91 - (5/7)(291/49) = {need:.4f}  => mu8 >= {mu_need:.4f}")
print(f"     (if J_8(C) >= 291/49; and mu8 is typically ~1.4-1.6, so mu8 >= {mu_need:.4f} is EASILY met)")
# verify on the actual wide family consec8 + far
J8,mu8=J_mu(list(range(1,9)))
Jwide,_=J_mu(list(range(1,9))+[99991])
print(f"     check: consec8 J_8={J8:.4f} mu_8={mu8:.4f}; wide {{1..8,99991}} J_9={Jwide:.4f} >= 432/91={thr:.4f}? {Jwide>=thr}")
print(f"     recursion floor: (5/7)(291/49)+(6/7)(1.4) = {(5/7)*consec8+(6/7)*1.4:.4f} vs threshold {thr:.4f}")
print()
print("If the recursion is exact (2), the k=9 TAIL reduces to [k=8 base J_8>=291/49] + [mu_8 bounded")
print("below, easy] => J_9(wide) >= 432/91 with margin. The far-element decorrelation RAISES J via the")
print("affine map J -> (5/7)J + (6/7)mu (a contraction toward the decorrelated fixed point).")
print("\ndone.")
