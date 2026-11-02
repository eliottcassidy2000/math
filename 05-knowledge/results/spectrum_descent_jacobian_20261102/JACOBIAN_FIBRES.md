# Algebraic inverse, integer image, actual source

**PROVED linear-algebra comparison; CITED centralizer theorem; no planar
Jacobian-conjecture proof or equivalence with tournament completeness.**

## 1. The exact rank-one loss in a fixed block response

For a nonempty actual block B with profile(H,E,p,Q,s,t,F_B), put a=H-s and
b=F_B-H. Its inherited map is

`M_B(F,q,g,h)=(HF+pq+ag+bh, Eq+Qh, (H-t)h, Hh).`              (J1)

H and E are positive. Over Q the matrix has rank3, with one-dimensional
kernel. If d=gcd(H,a), its integer kernel is precisely

`Z * (-a/d,0,H/d,0).`                                        (J2)

**Proof.** A kernel vector has h=0 from the fourth row, q=0 from the second,
and HF+ag=0 from the first. Dividing by the gcd gives the primitive generator.
The third row is then automatic. Every nonempty integer response fibre is
one affine translate of (J2). QED. This describes ambient integer fibres;
intersection with the actual source is a separate step.

Instead of overwriting g, retain it:

`L_B(F,q,g,h)=(HF+pq+ag+bh, Eq+Qh, g, Hh).`                  (J3)

This triangular matrix has determinant H^2 E. It is an explicit linear
polynomial automorphism over Q. Given (u,v,w,z), its inverse is

`h=z/H, q=(v-Qh)/E, g=w, F=(u-pq-ag-bh)/H.`                   (J4)

The overwritten output (H-t)h is recovered using the known block and h.
Equivalently, adding old g as a fifth observable repairs (J1)'s rank loss.
One independent scalar is necessary over Q; this is not a one-bit bound or
a minimal carrier theorem for literal metagraph operations.

The integer image has index H^2 E, despite the rational inverse. Formula
(J4) supplies three successive divisibility guards. Actual marked states
must additionally obey positivity, incidence, and recursive source membership.
These are not consequences of the determinant.

## 2. Exact integral guards for D

For D, H=5,E=32. Remove the unit g pivot in (J3). The remaining matrix is

`[5 p b; 0 32 Q; 0 0 5].`

Its first determinantal divisor is1. Its second is

`d2=gcd(5,pQ-32b)`,

because the other two-minors include25 and160 and the remaining relevant
ones are multiples of5. Its determinant is800. Thus the full four-coordinate
Smith factors are

`(1,1,d2,800/d2),  d2 in {1,5}.`                             (J5)

All-minor enumeration for all16 cuts gives12 profiles of type(1,1,1,800)
and4 of type(1,1,5,160). These are exact finite certificates of the general
formula, not fits to modular samples. The corresponding lifted singleton
map has determinant2 and Smith factors(1,1,1,2).

## 3. Actual sources leave a hole in an integer fibre

For S0 at height5, the target (15,10,0,5) has these three predecessors
among the necessary phase and cut/interior guards:

`(5,5,5,5),  (7,5,3,5),  (9,5,1,5).`                        (J6)

They lie consecutively on the primitive response-kernel line(-1,0,1,0),
with step twice that primitive vector because the odd-F guard removes the
intervening parity class. The endpoints are actual, constructed respectively
by D_0 S1 and D_1. The middle is impossible because H=7 never occurs.
All three map to the same target under the integer formula.

This proves an exact nonconvexity of actual source membership on a response
fibre. Even its parity sublattice, positivity and phase inequalities do not
fill the hole. The companion (9,7,3,5) obstruction in [DESCENT.md](DESCENT.md)
shows that excluding the two forbidden scalar counts is still insufficient.

## 4. What transfers from the older Jacobian work

[THM-2230 / planar-jacobian-response-fiber-and-exact-target-shear-quotient](../../../01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md)
uses the Cheng–McKay–Wang theorem, extended to characteristic-zero fields by
[Moskowicz, Theorem2.3](https://arxiv.org/abs/1802.04685v2). **Assuming** P
already has a constant nonzero-Jacobian mate, the operator

`D_P(Q)=Jac(P,Q)`

has kernel k[P]. Consequently a nonempty response fibre is Q0+k[P]. This
classifies fibres over the operator's image; it does not construct the
assumed mate or turn every requested response into an attained one.

The source-to-target comparison is literal at the level of response maps:

| Object | Map | Invisible directions | Additional gate |
|---|---|---|---|
| Polynomial complement Q | Q -> Jac(P,Q) | k[P], under the mate hypothesis | requested response belongs to the image; polynomial/effectivity conditions |
| Tournament block state | (J1) | primitive line (J2) | integer congruences and actual-source realization |
| Lifted tournament state | (J3) | none over Q | lattice index H^2E, then source realization |

The comparison preserves the distinction between a fibre and its image. It
does not map polynomial pairs to tournaments or transfer a conjecture between
them. Incoming [Session142's bridge](../session142_tournament_odd_spectrum_20261102/jacobian_bridge/PROOFS.md)
already shows that even unimodular normalized one-arc maps can invert a
forbidden count outside the Boolean source, and that positive grading needs
a separate legality condition. Here (J6) strengthens the operational lesson:
the relevant discrete source can be punctured inside a feasible integer fibre.

## 5. Why the amplifier coordinates help, but do not close the problem

For the empty cut, a=H, b=0, p=E-H, Q=0. On the invariant g=0 slice, using
c=q-F turns (J3) into the diagonal action

`(c,q,h) -> (Hc,Eq,Hh).`

The change from (F,q,h) to (c,q,h) is integral unimodular and the determinant
remains H^2E. At the initial seed, use c=q-F-g until the first block resets g.
This explains the uniform family algebraically while preserving its paid
endpoint contribution. It does not remove lattice or source-entry debt.

A promising next target is an operation atlas whose inverse branches have
**proved coverage of actual source sets**, with a decreasing well-founded
clock. The required coverage is weaker than unique inversion and stronger
than an ambient Jacobian criterion. The compressed grammar inverse makes
that distinction executable, but uniform coverage of the two prime lanes
is still unproved.
