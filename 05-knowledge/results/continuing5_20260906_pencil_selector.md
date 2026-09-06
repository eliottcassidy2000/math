# Middle resultant roots select the exact anchored coefficient fibre

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The new result gives an ordered endpoint selector for every
fixed coefficient pair having a strict reference point. It also provides
exact controls refuting zero-based fibres and nesting of the two channels.
It does not claim that every coefficient pair has a reference point or prove
the remaining original-phase response sign.

## 1. Inheritance and the missing order

The closest proved mechanism is the degree-eight quotient-algebra decoder
and compact-fibre corollary in
`05-knowledge/results/continuing4_20260906_moments_packet.md` with its audit.
The incoming `long_frontier_sep06_fibre_domain.md` identifies the exact
intervals at(x,y)=(84,35) using ordered interlacer roots and their resultants.
The present theorem extracts the general selection operation from that
single fibre. Its canonical hostile is a positive determinant with the wrong
matrix inertia. The corrected near miss is treating the D-boundary that
dominates at(84,35) as a universal channel nesting. The underused sidecar is
the inertia of the derivative of the Hankel pencil in its fifth coefficient.

The live objects are: the two moment functionals; their common denominator;
the slope matrix; the ordered resultant zeros; the weak model boundary; and
the original carried response. The chosen method cards are retaining the
product/shared labels and typing a residual before calling it one item.

Source: both native moment matrices. Target: two ordered quartic root lists.
Map: take determinant after identifying the slope inertia at a strict
reference point. Preserved: exactly the PSD interval, including endpoints.
Lost if one retains only the determinant sign: which positive component
contains the admissible matrices. Needed sidecars: a strict reference,
the two-positive/two-negative slope inertia, and nonnegative z. Cheapest
tests: actual channel reversal and a positive lower endpoint, supplied below.

## 2. The exact selector

Fix x,y>=0 and define

    F(v)=v^5-13v^4+55v^3-xv^2+yv, B_z(v)=F(v)-z,
    C(v)=v^4-12v^3+45v^2-(2/3)xv+(3/7)y,
    D(v)=v^4-11v^3+36v^2-(5/12)xv+(1/7)y.

Let H_A(z) be the ordinary5-by-5 Hankel of the formal moments of A/B_z,
through moment8, for A=C,D. Suppose a real z_*>0 has both H_C(z_*),H_D(z_*)
positive definite. The same statements hold for a reference z_*=0 when its
two Hankels are positive definite. This is an explicit hypothesis, not an
assumed supplier for every x,y.

For A=C,D let

    R_A(z)=Res_v(B_z,A)=det H_A(z).

These are monic quartics with four real roots, counted with multiplicity,

    a_(A,1)<=a_(A,2)<z_*<a_(A,3)<=a_(A,4).

Then, for arbitrary real z,

    H_A(z)>=0  iff  a_(A,2)<=z<=a_(A,3).              (1)

Consequently the complete weak two-interlacer model fibre is exactly

    [ max(0,a_(C,2),a_(D,2)), min(a_(C,3),a_(D,3)) ]. (2)

The reference guarantees this interval is nonempty. If z_*>0 it has
positive length. The endpoints include weak/canceled/multiple-root cases
through the inherited exact decoder. No assertion that the lower endpoint
is0, or that the same channel supplies both endpoints everywhere, is made.

## 3. Proof: a fixed inertia identifies the middle positive component

Write m_j for either formal moment sequence. Moments0 through4 are
independent of z. The exact z coefficients of m5,...,m8 are

    C: (1,14,130,904+4x/3),
    D: (1,15,147,1067+19x/12).

These are inherited recurrence identities. Thus H_A(z) is affine in z.
Its derivative has a zero first row and column and the remaining block

    S(a,b,c) = [ 0  0  0  1 ]
               [ 0  0  1  a ]
               [ 0  1  a  b ]
               [ 1  a  b  c ].

Its determinant is1 for every a,b,c. Along the path replacing(a,b,c)
by(t a,t b,t c), the block never becomes singular. At t=0 it is the
reversal matrix, with two positive and two negative eigenvalues. Hence
the full slope has inertia(2 positive,2 negative,1 zero) in both channels.
This is an analytic all-parameter argument, not a numerical eigenvalue test.

Congruence by H_A(z_*)^(-1/2) gives

    H_A(z) ~ I+(z-z_*) T_A,
    T_A=H_A(z_*)^(-1/2) H'_A H_A(z_*)^(-1/2).

The real symmetric T_A has the same inertia. Each of its two positive
eigenvalues contributes one zero of the determinant below z_*; each of
its two negative eigenvalues contributes one above. Its zero eigenvalue
contributes the constant factor1. Thus there are precisely four finite
real determinant zeros with the stated ordering. Every eigenvalue of the
affine normalized matrix is nonnegative exactly between the largest lower
zero and the smallest upper zero. This proves(1), including equality and
coincident endpoint roots. The determinant is monic: its leading z^4
coefficient is m0 det S=1.

To identify the determinant, first suppose B_z has five distinct roots b_i.
Partial fractions give w_i=A(b_i)/B'_z(b_i). If V is the Vandermonde matrix
with entries b_j^i, then H_A=V diag(w_i) V^T. Since
product B'_z(b_i)=(-1)^10 det(V)^2=det(V)^2,

    det H_A=product_i A(b_i)=Res_v(B_z,A).

The entries and resultant are polynomial in x,y,z over Q, so equality on
the dense simple-root set proves the identity everywhere. The sign is
positive for degree5; no absolute value or square root is substituted.
Equivalently, since A is monic quartic, R_A=product_(A(a)=0)(F(a)-z),
with multiplicities. The strict reference implies A has distinct positive
roots by the positive residue representation, but this alternate product
is not needed to choose the middle component.

Finally intersect(1) for the two channels with z>=0 and apply the proved
degree-eight iff decoder. This proves(2). At a selected boundary the Hankel
rank loss equals the multiplicity of that root of its determinant: in the
normalized affine matrix it is precisely the number of eigenvalue factors
that vanish. The corresponding cancellations and weak root contacts remain
part of the exact model, rather than being thrown away as singular cases.

## 4. Exact endpoint switches and a fibre that does not reach zero

The companion certificate stores all six monic resultant quartics, their
four isolated real roots, and the five positive reference minors in each
channel for these three pairs. The middle-root intervals have rational
endpoints; all displayed ordering comparisons are exact.

| Fixed(x,y) | Strict reference z | Upper C endpoint | Upper D endpoint | Active upper channel |
|---|---:|---|---|---|
| (155/2,9) | 1/200 | (9/869,2/193) | (124/507,79/323) | C |
| (155/2,37/4) | 1/5 | (130/489,109/410) | (232/897,217/839) | D |
| (311/4,21/2) | 1/5 | (298/867,309/899) | (177/530,176/527) | D |

For the last pair, C's lower endpoint is in(-4695/1259,-4516/1211),
whereas D's lower endpoint is in(163/2101,9/116). The complete model fibre
therefore has a strictly positive lower endpoint between those last two
rationals and an upper endpoint in(177/530,176/527). In particular it
contains z=1/5 and does not contain z=0. Propagating the model domain from
z=0 loses this entire fibre; a sign theorem on a containing coefficient
box may still legitimately use z=0 without asserting it is admissible.

The channel nonnesting has especially simple actual controls:

    (x,y,z)=(155/2,9,1/10): H_D>0,
        det H_C=-21973980341/10210252500<0;
    (x,y,z)=(155/2,37/4,13/50): H_C>0,
        det H_D=-11640479453184647/52276492800000000<0.

The surviving positive-definite channel already ensures five simple real
B roots; nonnegative coefficients put all of them in[0,infinity), and z>0
makes them positive. Thus these are genuine one-interlacer root models,
not nonreal-root surrogates. Neither individual channel implies the other.

The determinant alone also cannot supply(1). To the right of the largest
quartic root it is positive again, but the normalized slope calculation
shows two negative eigenvalues. The exact source probes both exterior
positive determinant components and rejects their Hankels. The middle-root
selection, not just resultant positivity, restores the predicate.

## 5. Reproduction and next question

The source independently builds the native recurrence, slope matrices,
Hankel determinants and literal resultants. Its finite universe is three
strict reference fibres, all six complete quartics, both exterior components,
and the two nonnested-channel controls. It uses exact rational Sylvester
minors and Sturm isolation, with no floating signs. All119 always-active
gates pass. The companion JSON retains every rational coefficient and interval.

    python -B 04-computation/continuing5_20260906_pencil_selector.py
    python -B -O 04-computation/continuing5_20260906_pencil_selector.py

An optional positional path writes the JSON elsewhere. The independent
referee must reconstruct the polynomials and matrix signs by a different
route; the universal inertia proof does not follow from this finite bank.

The remaining useful operation is to follow the selected middle roots as
(x,y) moves, retaining which channel is active and any degenerating strict
reference. Interlacer nesting and zero-based anchoring are unavailable
shortcuts. The separately proved rectangle sign argument may apply on a
larger containing box without needing either shortcut. No full phase-sign
theorem, global domain nonemptiness, or unproved determinant sufficiency is
claimed here.

Independent [proof and exact referee](continuing5_20260906_pencil_selector_audit.md) passes.

Filed checkpoint provenance: the [raw-byte manifest](continuing5_20260906_manifest.json)
pins the final report, source and output. Reviewed candidate-report hashes
above refer to the pre-promotion bytes. Source-location defaults and local
links were made portable where necessary; all emitted outputs were replayed
as raw bytes. The independent audit supplies the stated promotion basis.
