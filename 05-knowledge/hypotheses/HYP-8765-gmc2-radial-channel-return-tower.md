# HYP-8765 — the radial-channel return tower

**Status:** OPEN, sharply formulated; exact support evidence and two cutoff
falsifications.

**Owner:** codex-2026-07-21-S83.

**Target:** the strong two-dimensional nullcone statement NC2, hence GMC(2).

## Exact moment model

For

\[
P=\sum_{i=1}^k c_iZ^{a_i}W^{b_i},\qquad q_i=a_i-b_i,
\]

the (m)-th moment is

\[
\mathbb E[P^m]=
\sum_{\substack{r\in\mathbb N^k\\|r|=m,\ q\cdot r=0}}
\frac{m!}{\prod_i r_i!}
\left(\sum_i a_ir_i\right)!
\prod_i c_i^{r_i}.                                  \tag{1}
\]

Thus a balanced relation has two essential coordinates:

- charge (q_i=a_i-b_i), which decides whether it survives angular averaging;
- radial height (h_i=a_i+b_i), which decides its factorial weight after balance.

The old first-return quotient retained only the first coordinate and therefore
could not see cross-channel cancellation (MISTAKE-211).

## Corrected conjectural bound

For a positive-negative pair ((i,j)), define its primitive return length

\[
r_{ij}=\frac{q_i+|q_j|}{\gcd(q_i,|q_j|)},
\qquad
R(S)=\max_{q_i>0>q_j}r_{ij}.
\]

> **Radial-channel tower conjecture.** For a (k)-monomial support (S),
> every positive-negative pair satisfies
> \[
> (c_ic_j)^N\in
> \sqrt{\langle M_1,\ldots,M_{(k-1)r_{ij}}\rangle}
> \]
> for some (N\), where (M_m=\mathbb E[P^m]). Consequently
> (M^*(S)\le(k-1)R(S)).

If true, these pairwise radical inclusions cut the moment-null variety down to
the union of the two charge-one-sided coordinate subspaces and prove NC2.
The factor (k-1) is the rank of the integer relation lattice
\(\ker(q:\mathbb Z^k\to\mathbb Z)\); it replaces the false universal factor two.

## Exact evidence and falsification pressure

Exact rational/symbolic saturation sweeps give:

- every two-sided trinomial support of total degree at most eight closes at
  (2R): (8579/8579) supports;
- among four-monomial supports of degree at most four, (434/441) close at
  (2R), all seven failures close at (3R), and (441/441) close at (3R);
- among five-monomial supports of degree at most three, (64/65) close at
  (2R), and (65/65) close at (4R).

Two exact fakers show why the larger channel count is necessary.

### Return two, death at six

\[
P=W(a+cU)+Z(b+dU),\qquad U=ZW,qquad R=2.
\]

With (c=d=1), choose (a+b=A), (ab=B), where

\[
A=-\frac{15}{2}+\frac{3\sqrt3}{2}i,qquad
B=9-3\sqrt3 i.
\]

For (H=U(U^2+AU+B)),

\[
L(H)=L(H^2)=0,qquad
L(H^3)=3888i(3\sqrt3-2i)\ne0.
\]

Since \(\mathbb E[P^{2j}]=\binom{2j}{j}L(H^j)\), moments through
(2R=4) vanish and the obstruction first appears at (3R=6).

### Return three, death at nine

\[
P=aW+bW^4+ZW^2+Z^2,
\]

with charges ((-1,-4,-1,+2)) and (R=3). Let (a) be any root of

\[
f(a)=a^4-15a^3-240a^2-945a-1260
\]

and put (b=-a^2/12-a/2-1). Then (M_3=M_6=0), while the ninth-moment
factor, reduced modulo (f), is

\[
360(3560a^3+38530a^2+137613a+174552).
\]

Its resultant with (f) equals

\[
8335804094167484712960000\ne0,
\]

so every root dies at (M_9\). This refutes (2R) but matches (3R).

## Proposed proof mechanism

Localize at a pair product (c_ic_j\). At the levels

\[
r_{ij},2r_{ij},\ldots,(k-1)r_{ij},
\]

replace raw moments by successive resultants or logarithmic cumulants so that
decomposable lower relations are removed. Rank-one balanced flows decompose
into positive-negative transports. The desired remaining coefficient matrix
is a confluent factorial-Hankel/Vandermonde matrix indexed by radial channels.
Strict total positivity of the factorial Hankel kernel (((r+s)!)) is the
candidate nonzero-determinant input.

An equivalent Newton-face formulation places every monomial at
((q_i,h_i)). Maximize average height subject to average charge zero. Words
on the exposed balanced face share one factorial level; lower words are ordered
by face defect. The missing theorem is that successive defects give a triangular
channel filtration after earlier circuit ideals are removed. Pure face
dominance without this elimination is false because multinomial entropy can
populate an (O(\sqrt m)) boundary layer.

## Assumption challenge and Tournament Analysis

Candidate vertices considered were monomials, charges, primitive circuits,
radial channels, Newton-face defects, and proof obligations. Charges alone
preserve the angular predicate but destroy height and phase; circuits alone
still collide inside scalar equations. The least lossy working vertices are
radial channels/circuit obligations.

A tournament on runners or charges is **not** meaningful here: there is no
canonical antisymmetric binary observable whose orientation controls (1).
Forcing one would discard the many-to-one cancellation that is the problem.
If Tournament Analysis is used, it should be only as a proof-obligation
priority tournament (edge (A\to B) when discharging (A) is logically prior
to (B)); it is navigation, not mathematical evidence. The tie path is
`exact moment identity → localized circuit → channel determinant → pair radical
→ one-sided nullcone → GMC(2)`.

## Relation to THM-2014

THM-2014 proves a full infinite-degree slice by exactly the separation missing
from naive face dominance. For (P=aZ+b(U)+cW) with (d=\deg b\ge2), the
charge-zero endpoint beats the (k)-th charged channel by
(m^{-(2d-3)k}/(k!)^2), uniformly even for (k\asymp m). This establishes a
real entropy-versus-factorial theorem and identifies the general boundary:
channels within two degree units of the exposed face require resummation rather
than termwise dominance.
