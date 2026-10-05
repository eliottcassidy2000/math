# Level11: a modular dessin, a marked tournament, and the golden torsion clock

2026-10-05. **PROVED:** the elementary permutation, stabilizer and matrix
statements below. **CITED:** the modular-curve/newform identification and
the Frobenius characteristic equation. **FINITE-EXACT:** an independent
literal reader of the nine3-torsion points, and the declared finite groups.
No global Collatz coverage, modularity from numerical coincidence, or
identification of a field with a residue ring is asserted.

Artifacts: [script](../../04-computation/experiments/level11_dessin_golden_torsion_20261005.py)
and [output](level11_dessin_golden_torsion_20261005.out).

## 1. Inheritance and the useful distinction

The earlier [level11 short note](level11_short_20260922.md) already imports
the normalized eigenform

    f(tau)=eta(tau)^2 eta(11tau)^2
          =q product_(n>=1)(1-q^n)^2(1-q^(11n))^2,

its elliptic curve E:y²+y=x³−x²−10x−20, and the Hecke recurrence at2.
The current [five/eight/nine transfer](five_eight_nine_transfer_20261005.md)
proves the golden F9 clock and explains why it does not transfer a Collatz
guard. [THM-2768, modular C2–C3 quotients](../../01-canon/theorems/THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks.md)
supplies the closest group mechanism: keep the modular generators and
specify the extra quotient relation. The old exceptional degree11 action
appears in [the source interplay, M3](procgen_sources_20260923_interplay.md);
it must not be substituted for the degree12 action below.

The separate **PROVED** prime-absorption classification in
[crossings, section1.2](collatz_crossings_20260926_potential_and_seeds.md)
says that2,3,11 are exactly the primes p for which2p lies below the next
odd square above p. That assertion is about an inequality. No map from it
to modular level or ramification has been established.

Our board is **modular quotient / marked affine chart / eta differential /
local Frobenius / torsion coefficient field / coordinate field / source guard**.
The canonical hostile is equal clock order without conjugacy; the useful
sidecars are a chosen projective point, a basis, and the coefficient modulus.

## 2. The2–3–11 dessin is on the eta curve itself

**CITED.** Elkies identifies E with X0(11), associates the eta product, and
records its degree12 Belyi map with passport(2^6,3^4,11·1).
[Primary source](https://people.math.harvard.edu/~elkies/nature.html),
the conductor11 entry. Jones–Zvonkin independently describe the modular
dessin, its curve and its regular cover in sections11.3–11.4 of
[their paper](https://www.labri.fr/perso/zvonkin/Research/Klein-11.pdf).

Here is a direct finite model. On Omega=P1(F11)=F11 union{infinity}, put

    S(z)=−1/z,   T(z)=z+1,   R(z)=−1/(z+1)=S T(z).

Composition acts right to left. S²=R³=1 and S R=T. Their cycle types are

    S: 2^6,    R: 3^4,    T: 11·1.                  (1)

For instance S pairs0 with infinity; R cycles(0,10,infinity); T fixes
infinity and cycles the eleven finite elements. The fixed-point tests for
S and R are z²=−1 and z²+z+1=0. Their discriminants−1 and−3 are nonsquares
modulo11, so neither has a fixed point. This proves (1), without inferring
the cycle types from group orders.

The group generated is PSL2(F11). Translations T^b and their S-conjugates
are the two elementary unipotent families, which generate SL2(F11) by
row reduction; quotienting ±I gives the projective group. Its order is
11(11²−1)/2=660. The stabilizer of infinity has order55, so its coset
action has degree12. Passing from the modular group to this action is
exactly reduction modulo11 followed by forgetting the full level basis
and retaining its line. Equivalently the corresponding subgroup is
Gamma0(11), the inverse image of that stabilizer.

Treating the twelve symbols as dessin edges, the two vertex permutations
have6 and4 cycles and the face permutation has2 cycles. Hence

    chi=6+4−12+2=0,   genus=1.                       (2)

The face permutation is(SR)^−1=T^−1; its lengths equal those of T.
Our order of the two vertex colors exchanges the conventional0/1 branch
labels if needed. The standard j-map has order3 branching over0 and order2
over1728; normalizing or exchanging those two branch values changes no
surface or ramification multiplicity.

This gives an actual chain of maps:

    modular generators -> P1(F11) monodromy
        -> the dessin on X0(11)
        -> the one-dimensional space of holomorphic differentials,

in which the eta eigenform supplies the normalized differential
2*pi*i*f(tau)d tau=f(q)dq/q. The last identification uses the cited
modular-form input; it is not a deduction from a shared integer11.
The unmarked surface forgets the chosen map to the sphere, and the finite
permutation action forgets an analytic coordinate and q-normalization.
Retain the Belyi map and normalized cusp parameter to recover those data.

There are two nearby objects that must remain distinct. The regular action
on660 symbols has330,220,60 cycles, giving genus26; it is X(11), covering
X0(11) with degree55. The exceptional eleven-point action has an A5
stabilizer and produces Klein's planar degree11 dessins. Neither is this
twelve-point genus1 quotient. Also Delta(2,3,11) is infinite: its reciprocal
sum1/2+1/3+1/11=61/66 is below1. The finite group PSL2(F11) is a quotient,
not that triangle group itself.

## 3. Marking one point exposes an intrinsic Paley tournament

Retain infinity and restrict to the other eleven symbols. Its stabilizer
B consists of affine maps

    z -> a² z+b,  a nonzero, b arbitrary in F11.     (3)

It has55 elements. Define x->y exactly when y−x is a nonzero square.
Because−1 is nonsquare, every distinct pair has exactly one direction:
this is a tournament. Every map (3) preserves it, since it multiplies all
differences by a square. Thus the same Borel subgroup used by the modular
quotient acts on an actual eleven-vertex tournament.

The sidecar is essential. The full twelve-point projective action is
2-transitive and cannot preserve a tournament: it includes a permutation
interchanging any specified ordered pair. More concretely S takes1->2 to
10->5, whose difference6 is nonsquare. Keeping a projective point supplies
an affine chart; choosing squares as positive fixes the orientation gauge.
The two opposite orientations are both B-invariant. This assertion gives
an explicit automorphism subgroup; it does not use or claim a separate
classification of the full automorphism group.

The dessin is therefore not itself a tournament. Its twelve edge labels,
two vertex colors, rotations and two faces differ from the eleven affine
vertices and quadratic pair predicate. The map is to a marked chart of
the common projective action, with the marked point and gauge retained.

## 4. At prime2, the eta curve's3-torsion is the golden F9 module

Reduce E modulo2:

    E2: y²+y=x³+x².

There are four affine F2 points and the point at infinity, so a2=−2.
**CITED:** the Frobenius pi:(x,y)->(x²,y²) satisfies

    pi²+2pi+2=0.                                    (4)

For a prime ell different from the characteristic, its action on E[ell]
has characteristic polynomial T²−a2*T+2 reduced modulo ell.
See Sutherland's [endomorphism lecture, Theorems7.1,7.18,7.19](https://math.mit.edu/classes/18.783/2015/LectureNotes7.pdf)
for the arbitrary-characteristic torsion and trace/determinant statements,
and his [Schoof lecture, section9.2](https://math.mit.edu/classes/18.783/2015/LectureNotes9.pdf)
for the restricted characteristic equation. We use no odd-characteristic
coordinate formula from the latter lecture.

At ell=3 this becomes T²−T−1, irreducible over F3. Consequently

    F3[pi] = F3[phi]/(phi²−phi−1) = F9,
    phi maps to pi on E[3].                         (5)

For any nonzero P in E[3], P and pi(P) are a basis: dependence would
give an eigenvalue in F3, contradicting irreducibility. The explicit map

    a+b*phi -> aP+b*pi(P)                           (6)

is an additive bijection F9->E[3], intertwining multiplication by phi with
Frobenius. It depends on the choice of P. In that basis the matrix is
M=[[0,1],[1,1]]. Thus pi^4=−I and pi has order8 on E[3]; every nonzero
point has orbit8. Quotienting points by sign gives four cyclic subgroups,
and the projective period is4. This is an exact connection between the
prime2, the3-torsion coefficient field, the level11 eta curve and the
inherited golden clock.

The inherited Hecke companion A=[[-2,-2],[1,0]] has A^4=−4I. Reducing it
modulo3 gives A=J M J with J=[[0,1],[1,0]]. This is a matrix-level second
reader of (5). The Frobenius trace sequence and the Hecke prime-power
coefficient sequence have different initial conditions; we do not identify
their scalar outputs merely because they share a characteristic equation.

The field types matter. F9 in (5) is an algebra of operators on a torsion
group, not the coordinate field of those points. A direct tangent calculation
on E2 gives x(2P)=x(P)^4+1. Hence the nonzero3-torsion x-coordinates satisfy

    x^4+x+1=0.                                     (7)

This polynomial has no linear factor over F2 and is not divisible by the
only irreducible quadratic x²+x+1, so it is irreducible. Its four roots lie
in F16. The full eight points require F256: pi has order8 on each, while
pi^4 is their nontrivial negation (x,y)->(x,y+1). The literal reader below
constructs all these objects in one finite coordinate field.

**Hostile to naive lifting:** modulo9, A and M both have order24, but their
traces are7 and1 and their determinants are2 and8. Therefore they are not
conjugate over Z/9Z. An equal clock period does not retain the operator,
and the exact identification (5) does not automatically lift to9-torsion.
Likewise Z/9Z is not F9. The modulo9 inverse-ray translation decoder has
its own input and guard; (5) supplies no arithmetic Collatz transition.

## 5. Exact controls and scope

Run from the repository root:

    python 04-computation/experiments/level11_dessin_golden_torsion_20261005.py
    python -O 04-computation/experiments/level11_dessin_golden_torsion_20261005.py

The script writes no files. Its declared universe is all11^4 two-by-two
matrices filtered by determinant1; the complete660-element projective
group and all its twelve-point actions; all55 marked stabilizer maps on
55 arcs; all256² affine coordinate pairs over F256; all nine3-torsion
points and81 additions; and matrix orders bounded by8/24 at moduli3/9.
The irreducibility and point-division tests use a path independent of the
matrix clock. The finite-field polynomial X^8+X^4+X^3+X+1 is validated
by checking every nonzero element's255th power, not assumed from a label.

The controls recover225 points over F256, nine3-torsion points including
the identity, an eight-cycle on the nonzero points, the four x-coordinates,
the passport, and both genus calculations. The positive map is (6); the
negative controls are the unmarked Paley inversion and the mod9 trace/
determinant mismatch. None changes the original-source, binary valuation
or endpoint-payment obligations of a Collatz certificate.
