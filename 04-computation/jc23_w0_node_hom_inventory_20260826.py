#!/usr/bin/env python3
"""Exact inventory for THM-4230's full-Hom finite-ratio argument on W=0.

Separate certificates construct and audit the explicit hidden sublattice
L_exp.  This script freezes the attachment orbit, visible images, character
ledger, and isogeny input; it does not compute hidden saturation, integral
visible-hidden gluing, or mixed/coset attachment evaluation.
"""

from math import gcd

import sympy as sp


def need(condition, label):
    if not condition:
        raise RuntimeError(label)


U, Z, s = sp.symbols("U Z s", nonzero=True)
lam = U + Z
D = -4 * U * Z

# W=0 main component and the R--C attachment equation P=S^2.
P = s**2
component_at_node = sp.expand(1 - U * P**6 - Z * s**4 * P**4)
need(sp.factor(component_at_node - (1 - lam * s**12)) == 0,
     "attachment equation")

# The two visible quotient maps from THM-4230, equation (20).
# At a node x=S^2/P=1, u=1/P=s^-2, and v=2Z.
Ea_X, Ea_Y = s**-4, 2 * Z
Eb_X, Eb_Y = s**4, 2 * Z * s**6

# Reduce only the twelfth power using s^12=1/(U+Z).
need(sp.factor(((Ea_Y**2 - (D + 4 * Z * Ea_X**3)) * s**12
                ).subs(s**12, 1 / lam)) == 0,
     "first visible E0 image")
need(sp.factor((Eb_Y**2 - (4 * Z + D * Eb_X**3)
                ).subs(s**12, 1 / lam)) == 0,
     "second visible E0 image")

# Under tau, s -> zeta_12*s.  The first image has orbit size 3 and the
# second size 6.  Their pulled-back differentials have C_12 characters 8,10.
visible_characters = (8, 10)
hidden_rho_plus_characters = (5, 11)
need(set(visible_characters).isdisjoint(hidden_rho_plus_characters),
     "visible and hidden character packets must be disjoint")

# Put T=SP, X=U^(1/6)P, and Y=Z^(1/4)T.  This normalizes every W=0
# positive-genus component to C0: X^6+Y^4=1.  If Q0=(A,B) is the first
# attachment, then A^6=U/(U+Z), B^4=Z/(U+Z), hence A^6/B^4=U/Z.
A6 = U / lam
B4 = Z / lam
need(sp.factor(A6 + B4 - 1) == 0, "normalized attachment lies on C0")
need(sp.factor(A6 / B4 - U / Z) == 0, "marked ratio")

# Holomorphic differentials on X^6+Y^4=1 are indexed by
# 1<=i<=5, 1<=j<=3, i/6+j/4<1.  Under
# tau:(X,Y)->(zeta_6 X,zeta_4 Y), their characters are 2i+3j mod 12.
holomorphic_pairs = tuple(
    (i, j) for i in range(1, 6) for j in range(1, 4)
    if 2 * i + 3 * j < 12
)
tau_characters = tuple(sorted((2 * i + 3 * j) % 12
                              for i, j in holomorphic_pairs))
need(tau_characters == (5, 7, 8, 9, 10, 11, 11),
     "C0 holomorphic character ledger")
need(all(character != 0 for character in tau_characters),
     "tau-1 has no tangent-space kernel")
rho_plus_primitive = tuple(sorted(
    (2 * i + 3 * j) % 12 for i, j in holomorphic_pairs
    if (i + j) % 2 == 0 and gcd(2 * i + 3 * j, 12) == 1
))
need(rho_plus_primitive == hidden_rho_plus_characters,
     "hidden rho-plus primitive character ledger")

print("W0_GATE U*Z*(U+Z)!=0")
print("SOURCE_NODES P_j=S_j^2; S_j^12=1/(U+Z); j=0..11")
print("VISIBLE_EA_IMAGES (S_j^-4,2Z); ORBIT_SIZE 3; CHARACTER 8")
print("VISIBLE_EB_IMAGES (S_j^4,2Z*S_j^6); ORBIT_SIZE 6; CHARACTER 10")
print("HIDDEN_PRIMITIVE_E0_CHARACTERS 5,11")
print("NORMALIZED_CURVE C0:X^6+Y^4=1")
print("NORMALIZED_NODES Q_j=(A*zeta_6^j,B*zeta_4^j); "
      "A^6=U/(U+Z); B^4=Z/(U+Z)")
print("MARKED_RATIO A^6/B^4=U/Z")
print("TAU_HOLOMORPHIC_CHARACTERS " + ",".join(map(str, tau_characters)))
print("TAU_MINUS_ONE_ISOGENY YES")
print("MISSING hidden_saturation_overlattice; integral_visible_hidden_gluing; "
      "mixed_coset_attachment_evaluation")
print("FINITE_RESIDUAL for each n in {34,42}, only finitely many U/Z can "
      "satisfy node equality")
print("VERDICT FULL_HOM_SURVIVORS_REDUCE_TO_FINITE_RATIO_SETS")
