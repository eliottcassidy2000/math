#!/usr/bin/env python3
"""Hostile minimal checks for THM-4341.

The positive arithmetic identities are checked together with counterexamples
to stronger geometric, projective-closure, order, parity, and indexing claims.
This imports neither proposed certificate.
"""

from fractions import Fraction as F
from math import floor
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, label):
    if not bool(condition):
        raise AssertionError(label)


def odd_tail_data(m, r):
    eps = r & 1
    degree = m - r + eps
    genus = (degree - 1) // 2
    delta = r // 2
    return eps, degree, genus, delta


# Core odd-m identities.
for m in range(3, 42, 2):
    g = (m - 1) // 2
    for r in range(1, m):
        eps, degree, genus, delta = odd_tail_data(m, r)
        need(degree & 1, "odd branch degree")
        need(genus == (m - r) // 2, "tail genus")
        need(delta == r // 2 and genus + delta == g, "delta partition")
        need(genus == (m - r) // 2, "complement delta")
        need(F(r, m-r) * F(m-r, r) == 1, "reciprocal slope")


# There is no geometric isomorphism/birationality behind the numerical swap.
# For m=5, reciprocal orders 1 and 4 give genera 2 and 0.
data_51 = odd_tail_data(5, 1)
data_54 = odd_tail_data(5, 4)
need(data_51[2] == 2 and data_54[2] == 0, "different-genus reciprocal witness")


# "Projective closure is smooth" is false in ordinary P^2.  The quintic tail
# y^2=z(z^4-c) homogenizes to Y^2 W^3-Z^5+cZW^4; [0:1:0] is singular.
Z, Y, W, c = sp.symbols("Z Y W c")
plane_quintic = Y**2 * W**3 - Z**5 + c * Z * W**4
infinity = {Z: 0, Y: 1, W: 0}
need(plane_quintic.subs(infinity) == 0, "plane infinity point")
need(all(sp.diff(plane_quintic, v).subs(infinity) == 0 for v in (Z, Y, W)),
     "ordinary plane closure is singular")


# Even m is structurally different: the normalized branch degree is even and
# genus+horizontal-delta is one short.  Its two infinity places make attachment
# incidence a necessary extra sidecar; the polynomial alone does not determine
# whether the missing unit is realized as a graph cycle.  At r=m/2 reciprocity
# has a fixed point and the proposed odd centered orientation is zero.
for m in range(4, 24, 2):
    g_even = m // 2
    for r in range(1, m):
        eps = r & 1
        degree = m-r+eps
        need(degree % 2 == 0, "even-m branch degree")
        genus = (degree-2)//2
        delta = r//2
        need(genus + delta == g_even-1, "even-m missing graph unit")
    rfix = m//2
    need(2*rfix-m == 0 and F(rfix,m-rfix) == 1, "even reciprocal fixed point")


# Only the excess beyond k is reciprocal.  Full orders do not have the stated
# constant product when k is nonzero.
chi, d, lam, k = F(1,2), 7, F(1,2), 5
ex1, ex2 = chi*d*lam, chi*d/lam
ord1, ord2 = k+ex1, k+ex2
need(ex1*ex2 == (chi*d)**2, "excess product")
need(ord1*ord2 != (chi*d)**2, "full-order product counterexample")

# Comparing separately minimal clearing bases also changes units.  A product
# identity needs one common pre-clearing valuation normalization.
m, r = 5, 1
n = m-r
d_r, d_dual = 2*n, 2*r
minimal_ex1 = chi*d_r*F(r,n)
minimal_ex2 = chi*d_dual*F(n,r)
need((minimal_ex1, minimal_ex2) == (F(1), F(4)), "minimal-base excesses")
need(minimal_ex1*minimal_ex2 != (chi*d_r)**2,
     "separate minimal bases have no common-d product")


# A bare odd square is only an in-block rank.  Square 1 recurs for every g at
# r=g or g+1; the triangular block offset is indispensable globally.
square_one_rows = []
for g in range(1, 8):
    m = 2*g+1
    for r in (g, g+1):
        h = (abs(2*r-m)+1)//2
        need((2*h-1)**2 == 1, "repeated odd-square label")
        N = g*(g-1)//2+h
        square_one_rows.append((g,r,N))
need(len({N for _,_,N in square_one_rows}) == 7,
     "triangular block offset separates repeated square labels")
need(F(1,2) == F(3,6), "reduced slope loses total scale")


# Scope assumptions are load-bearing.
z = sp.symbols("z")
char3 = sp.Poly(z**3-1,z,modulus=3)
need(sp.gcd(char3,char3.diff()).degree() > 0, "inseparable characteristic hostile")
zero_lead = sp.Poly(z*(z**4),z,domain=sp.QQ)
need(sp.gcd(zero_lead,zero_lead.diff()).degree() > 0, "c=0 hostile")


print("PASS THM4341 HOSTILE MINIMAL CHECKS")
print("NUMERICAL_CORE=odd_delta_genus_swap_and_slope_reciprocity_PASS")
print("NO_GEOMETRIC_ISOMORPHISM=m5_r1_genus2_vs_r4_genus0")
print("PLANE_CLOSURE_HOSTILE=degree5_infinity_singular;use_smooth_projective_model")
print("EVEN_BOUNDARY=genus_plus_delta=g-1;two_infinity_places_need_attachment_sidecar;fixed_r=m/2")
print("ORDER_FIREWALL=only_excess_in_common_valuation_units_has_constant_product")
print("INDEX_FIREWALL=odd_square_is_in_block_only;triangular_g_offset_and_orientation_required")
print("SCOPE=char0_or_char_not_2_and_not_dividing_m-r;c_nonzero")
