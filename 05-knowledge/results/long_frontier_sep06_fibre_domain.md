# Exact admissible intervals on the four-anchor Laurent fibre

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This identifies the entire weak-interlacing domain of the coefficient
fibre used by [the all-root continuation](long_frontier_sep06_residue_tail.md).
It is a domain theorem, not itself a sign theorem for the doubled response.

## 1. Exact statement

For real f define

    F(v)=v^5-13v^4+55v^3-84v^2+35v,
    B_f(v)=F(v)-f,
    C(v)=v^4-12v^3+45v^2-56v+15,
    D(v)=v^4-11v^3+36v^2-35v+5.

Let rho be the unique root in (4,41/10) of

    R_C(f)=f^4+6f^3-285f^2-778f+7125,

Let kappa be the unique root in (43/10,22/5) of

    R_B(f)=3125f^4+125758f^3-4825473f^2-30354226f+211485225,

and let L=14sqrt(2)-16. The exact hierarchy is

    0<L<19/5<4<rho<41/10<43/10<kappa<22/5<5.

Then:

1. B_f has five nonnegative real roots **if and only if** 0<=f<=kappa.
2. B_f has five nonnegative real roots and C weakly interlaces B_f
   **if and only if** 0<=f<=rho.
3. B_f has five nonnegative real roots and D weakly interlaces B_f
   **if and only if** 0<=f<=L.
4. The simultaneous C,D-interlacing domain is exactly [0,L].

All five B roots are simple throughout the two interlacing domains. The B roots
are strictly positive except for one simple zero at f=0. C strictly
interlaces B on [0,rho), including the whole simultaneous domain;
at f=rho it shares only its smallest root with B. D strictly interlaces
B on [0,L), and at f=L it shares only its smallest root with B.

On the larger B-only domain, the first two B roots coalesce into one
double positive root at f=kappa; the other three remain simple.

Thus the D determinant boundary found in the residue continuation is an
attained exact boundary of actual admissible shapes. It is not merely
a scalar necessary bound on a possibly empty class. The C-only fibre
extends to rho>4, strictly beyond the D boundary.
The whole B-only domain lies in [0,5]. Therefore the independently audited
[all-root sign theorem](long_frontier_sep06_residue_tail.md) applies to every
nonnegative-root B with these four anchors, without an interlacer assumption.
This implication requires the separate sign theorem; it is not proved
from domain geometry alone.

## 2. Inheritance and connection

The inherited fixed actual point is f=1, with e1,e2,e3,e4=(13,55,84,35).
The new [residue continuation](long_frontier_sep06_residue_tail.md)
obtains the D bound through its fifth moment determinant. Root independently
recognized that determinant as a resultant and factored D into two
quadratics. The resulting exact root values turn the determinant bound
into a complete domain classification.

The live concepts are coefficient fibres, rational residue measures,
resultants, ordered root intervals, weak boundary collisions, and the
original response phase. The closest proved mechanism is the positive
residue envelope of [the C-only first-phase theorem](long_frontier_sep06_anchor.md).
Its repeated (0,0,3,5,5) shape is the canonical hostile to discarding D.
The corrected near miss is treating one positive moment determinant as
sufficient for interlacing. Here sufficiency comes from explicit root
ordering and sign changes, with all five B roots constructed.

The source-to-target map takes the original coefficient f to the four
values B_f at the fixed interlacer roots. It preserves the full monic
coefficient law and root order. The resultant alone forgets which root
collides first; that ordered-root sidecar is retained below. Targeted
searches at HEADa4aebc127e found no prior exact domain statement or the
same radical boundary in the relevant canon and Laurent/anchor reports.
No external priority claim is made.

## 3. C gives a complete root construction

Let c1<c2<c3<c4 be the roots of C. The following four disjoint rational
intervals contain one root each, as follows from alternating endpoint
signs of the quartic. Exact rational interval Horner evaluation gives
the displayed F bounds:

| Root | Isolating interval | Enclosure of F at that root |
|---|---|---|
| c1 | (136/373,101/277) | (4,41/10) |
| c2 | (461/262,783/445) | (-7,-6) |
| c3 | (2114/537,2425/616) | (14,15) |
| c4 | (2245/378,3807/641) | (-19,-18) |

Four distinct sign changes exhaust the degree, proving that these are
all roots and they are simple and positive. The exact resultant is

    Res_v(B_f,C)=product_(j=1)^4 (F(c_j)-f)=R_C(f).     (1)

The four F values lie in disjoint intervals in the table. Therefore
rho=F(c1) is indeed the unique root of R_C in (4,41/10), with no
unselected conjugate or branch ambiguity.

For 0<f<rho, the five signs at 0,c1,c2,c3,c4 are

    B_f(0)<0, B_f(c1)>0, B_f(c2)<0,
    B_f(c3)>0, B_f(c4)<0.

Together with B_f(v) tending to positive infinity as v tends to
positive infinity, these give one B root in each of the five intervals
(0,c1),(c1,c2),(c2,c3),(c3,c4),(c4,infinity). Degree exhaustion
gives five simple positive roots and strict C interlacing.
At f=0, replace the first open interval by the simple zero at0;
B_0'(0)=35, and the other four intervals remain strict.

At f=rho, take the continuous limits of these real roots. No boundary
except c1 is a root of B_rho. Exact interval evaluation on the first
box gives F'(c1)<-6. Thus the two adjacent roots cannot both collapse
to c1: that would be a multiple B root. All five roots remain simple
and nonnegative, and the second B root is c1, because its derivative
is negative. This proves weak C interlacing at the endpoint.

Conversely, if B_f has nonnegative roots then f is their product and
f>=0. If C weakly interlaces B_f, its smallest root lies between the
first two B roots, so B_f(c1)>=0. Hence f<=F(c1)=rho. This proves
both directions and the stated C boundary behavior.

## 4. D selects an exact subinterval

The fixed quartic factors as

    D(v)=(v^2-6v+1)(v^2-5v+5).

Its four ordered positive roots and corresponding F values are

| Root d_j | F(d_j) |
|---|---|
| 3-2sqrt(2) | -16+14sqrt(2)=L |
| (5-sqrt(5))/2 | (15-15sqrt(5))/2<0 |
| (5+sqrt(5))/2 | (15+15sqrt(5))/2>4 |
| 3+2sqrt(2) | -16-14sqrt(2)<0 |

These are exact radical identities. For 0<=f<L, B_f at the four
d_j has signs +,-,+,-. As before the signs construct strict
interlacing with five positive B roots (with the zero exception at f=0).
Alternatively the C construction already supplies the five simple roots,
since L<rho. Continuity at f=L gives weak D interlacing. Only d1
can be shared, because the other three F values differ from L; all B
roots stay simple by the C construction.

If D weakly interlaces any nonnegative-root B_f, the same necessary
smallest-root sign gives B_f(d1)>=0, hence f<=L. This proves the
exact D domain and its simultaneous intersection with the C domain.

For the strict comparison L<19/5, square the equivalent positive
inequality 70sqrt(2)<99: 9800<9801. The other displayed order bounds
follow from the elementary radical intervals checked in the source.

## 5. Why the fifth moment determinant found this boundary

For distinct roots a_i of a monic degree-five B, write the partial
fractions D/B=sum_i w_i/(v-a_i), where w_i=D(a_i)/B'(a_i).
Let nu_j=sum_i w_i a_i^j. The five-by-five moment matrix factors as

    (nu_(i+j))_(i,j=0)^4 = Vandermonde(a)*diag(w)*Vandermonde(a)^T.

Since product B'(a_i)=Vandermonde(a)^2 in degree five, taking the
determinant gives the exact polynomial identity

    det(nu_(i+j)) = product D(a_i)=Res(B,D).

Both sides extend across repeated roots as polynomial identities.
For this fibre they equal

    (f^2-15f-225)(f^2+32f-136).                     (2)

Positivity of this single determinant does not generally characterize
interlacing. The ordered-root table identifies the relevant first
collision and proves that its allowed side really is admissible.
This is the additional information needed to turn (2) into an iff.

## 6. The full nonnegative-root domain, without either interlacer

The derivative

    F'(v)=5v^4-52v^3+165v^2-168v+35

has four simple positive roots eta1<eta2<eta3<eta4. Four disjoint
sign-changing boxes and exact interval evaluation give:

| Critical point | Isolating interval | Enclosure of F at that point |
|---|---|---|
| eta1 | (131/472,68/245) | (43/10,22/5) |
| eta2 | (564/401,737/524) | (-10,-9) |
| eta3 | (1521/457,1957/588) | (26,28) |
| eta4 | (4391/815,4655/864) | (-63,-62) |

The derivative has signs +,-,+,-,+ on its five monotonicity intervals.
Also F(0)=0, and F tends to positive infinity at positive infinity.
For 0<f<F(eta1), its five monotone branches each cross the level f
once, giving five simple positive B_f roots. At f=0 the first root
is the simple zero at0. At f=F(eta1), the first two crossings merge
at eta1; since F' has a simple root there, this is exactly a double
B root. The other three crossings remain strict.

Conversely, if B_f has five nonnegative roots then f>=0. By Rolle's
theorem with multiplicities, its derivative roots interlace them.
In particular eta1 lies between the first two B roots, so
B_f(eta1)>=0 and f<=F(eta1). This proves the full iff.

The exact resultant of B_f and F' is R_B(f). It is 5^5 times the
product of F(eta_j)-f, so the four enclosures above identify
kappa=F(eta1) as its unique root in (43/10,22/5). The endpoint
argument retains repeated B roots at the discriminant boundary.

## 7. Verification and scope

The exact source checks the four C root boxes and signs, their interval
F enclosures and derivative boundary, both resultant identities, the
D factorization and radical values, the full B critical-point table and
discriminant, and the moment determinant through independent formal
series division. The full interval theorem is proved
by sign changes, degree exhaustion and the displayed endpoint arguments;
no shape scan is promoted to a universal statement.

The [independent analytic/source audit](long_frontier_sep06_fibre_domain_audit.md)
passes without correction. The [source](../../04-computation/long_frontier_sep06_fibre_domain.py)
has42 always-active exact gates; normal and optimized output agree with
the [frozen output](long_frontier_sep06_fibre_domain.out). The auditor
separately reconstructed the critical-point boxes, image intervals and
discriminant, and audited all iff and endpoint arguments.

```bash
python3 -B 04-computation/long_frontier_sep06_fibre_domain.py
python3 -B -O 04-computation/long_frontier_sep06_fibre_domain.py
```

Raw LF source SHA256:
dacdf6a3323cce90e14002a03bcbd87b105790717fdf243db754c06ef1b93f29.
Raw LF output SHA256:
41d5f6db13170eac6926051272ad429326e0b3e0ce8fb9c70ac5ec62d828e6bf.

All admissible-domain claims retain the
four original coefficient anchors and the prescribed contiguous rows.
They do not assert that a free coefficient shape is an actual binomial
row, nor infer any doubled-response sign without the separate consumer.
