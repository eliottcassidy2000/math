# Exact binomial radical differentials in every exponent

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The theorem classifies algebraic exactness for the stated binomial family in
all exponents. Its proof is elementary and retains the actual differential,
normalization and coefficient field. This is a recovered classical-integration
mechanism in the workspace's exactness language, not an external-priority
claim. No general Jacobian-mate or planar Jacobian conclusion is asserted.

## 1. The complete criterion

Fix integers A>=0 and r>=1 and complex a,b with ab nonzero. Put u=x-p for a
fixed p, and define

    N=u^A(a*u^r+b),       y^2=N,       K=C(u)(y),
    eta=du/y.

The residual binomial has r distinct nonzero roots. In particular N is never
a square under these hypotheses, even when u^A has a large square factor.
The actual smooth projective curve of K is connected.

**Theorem.** The following are equivalent:

1. eta has a primitive in K.
2. eta has a primitive in a finite algebraic extension of K.
3. Either `(A,r)=(0,1)`, or

       A=r+2+2r*ell,       ell an integer >=0.           (1)

The primitive in the positive cases has a formula rational in a,b,u,y; no
additional radical of a or b is necessary. Consequently every nonzero member
of each two-dimensional space

    span{1,u},
    u^(r+2+2r*ell) span{1,u^r},  r>=1,ell>=0,           (2)

is exact, with the coefficient-zero endpoints treated separately as pure
powers. No claim that (2) exhausts all exact pencils of unbounded degree is
made. In degrees at most eight, it recovers all six spaces of the proved
[exact_pencils theorem](planar_jc48_sep08_exact_pencils.md).

On the staircase (1), the genus is `floor((r-1)/2)`, independent of ell.
Increasing ell raises the total primitive pole degree to `r(2ell+1)` while
preserving the radical curve. The general rejected inputs do not all have
that genus: their branch count is

    r + (A mod2) + ((A+r) mod2).                        (3)

For example A=1,r=2 has genus one. Formula (3) prevents incorrectly assigning
the staircase genus to every binomial with the same r.

The closest proved mechanisms are the all-degree odd-purity and pole-space
arguments in [genus_capacity](planar_jc48_sep08_genus_capacity.md), together
with the coefficient-rational construction in exact_pencils. The canonical
hostile is the logarithmic quadratic `u^2+1`. The corrected near miss is
preservation of genus or even the whole radical field while changing the
weight of the differential: `u^5(u^3+1)` is exact, `u^7(u^3+1)` is not.
The least-used sidecar is the decomposition of the derivative operator by
exponents modulo r.

The live concepts are the actual weighted radical field; regular functions
on its complete affine model; exponent residue classes; lowest/highest
coefficient compatibility; classical logarithmic terms; and parameter descent.
The source-to-target map is inversion with its differential factor retained,
followed by the odd part under the radical involution. It loses additive
constants only. Retaining just the curve would lose the numerator exponent;
the same-field hostile is the cheapest test of that mistake.

## 2. Classical context and a necessary distinction

The primary antecedent is P. Tchebichef,
[*Sur l'intégration des différentielles irrationnelles*](https://www.numdam.org/item/JMPA_1853_1_18__87_0.pdf),
Journal de mathématiques pures et appliquées, first series18 (1853), pp87–111.
I inspected §I, pp87–88, where algebraic and logarithmic terms are separated,
and §VIII, pp106–108, which addresses rational-exponent binomial differentials.
That classical integration problem includes logarithms. The criterion here
requires an algebraic primitive, so elementary integrability alone is
insufficient. The operator proof below is self-contained; no unread version
of a classical theorem is imported as an iff for the stronger predicate.
The map to the classical notation is
`u^(-A/2)(b+a*u^r)^(-1/2)du`.

The inherited [leading_exactness](planar_jc48_sep08_leading_exactness.md)
and [boundary_exactness](planar_jc48_sep08_boundary_exactness.md) distinguish
this same algebraic-exactness gate from elementary expressions or formal
antiderivatives. Targeted exact and synonym searches of current canon and
session results recovered those mechanisms and the low-degree sparse cases,
not a separate current all-exponent binomial criterion. This is not a claim
of external novelty; the classical source is explicitly credited.

## 3. Low degree and odd total degree

The total degree is n=A+r, with no leading cancellation because a is nonzero.
If n=1 the unique pair is A=0,r=1 and

    B=2y/a,        dB=du/y.                            (4)

If n=2, the only pairs are (0,2) and (1,1). Each polynomial has two distinct
roots. At its two infinity points, with v=1/u and y=v^-1*w, the function w
has nonzero limits of opposite sign. Hence

    eta=-dv/(v*w)

has nonzero simple-pole residues and cannot be exact.

For odd n>=3, the proved all-degree purity theorem says an exact polynomial
must have only one distinct root. If A=0, there are r=n>=3 distinct roots.
If A>0 there is the root u=0 and at least one distinct nonzero residual root.
Thus no such binomial is exact. Multiplicity-two logarithms are already part
of that purity theorem's hypotheses and proof; they are not silently allowed.

It remains to classify even n=2k>=4.

## 4. A complete affine primitive space

Use the exact field change

    X=1/u,       Y=y/u^k,
    Y^2=P(X):=a+b*X^r,
    eta=-X^j dX/Y,       j=k-2>=0.                     (5)

This changes neither the field nor exactness. In particular the square factor
in u^A is absorbed by an actual rational field coordinate, not discarded from
the differential. The curve `Y^2=P(X)` is nonsingular affine: P and its
derivative are coprime since a,b are nonzero and characteristic is zero.

The differential in (5) is regular everywhere on this affine curve. At a
simple zero of P, Y is a local parameter and dX/Y is a unit; at X=0 the
nonzero a prevents a pole, and j is nonnegative. Therefore a rational
primitive has no affine poles. Its membership in the full normal affine
coordinate ring gives the unique form

    A0(X)+Y*B(X),       A0,B in C[X].

The differential is odd under `Y -> -Y`. Taking the odd part of a primitive
preserves its derivative, so an exact differential has a primitive Y*B(X)
with polynomial B. There is no a priori degree truncation or residue-only
claim in this reduction. Conversely a polynomial B with the stated derivative
is already a primitive in the actual field.

The equation for B is

    L(B):=(a+b*X^r)B' + (r*b/2)X^(r-1)B = -X^j.        (6)

For each q>=0 the complete monomial image is

    L(X^q)=a*q*X^(q-1)+b*(q+r/2)*X^(q+r-1),             (7)

where the first term is zero if q=0. In particular L has no nonzero
polynomial kernel: a polynomial of degree M has image of degree M+r-1,
with nonzero leading coefficient b(M+r/2) times its own leading coefficient.

## 5. Highest and lowest coefficients give the iff

Suppose a solution B exists. Its forced degree is

    M=j-r+1 >=0.                                      (8)

The operator sends an input exponent class h modulo r into output class h-1
modulo r. Decomposing B into these classes, every class except the one
congruent to j+1 would have zero image; injectivity from (7) makes those
classes zero. This applies also to r=1, when there is just one class.

Let q0 be B's smallest nonzero monomial exponent. If q0>0, its first term
in (7) gives a nonzero coefficient `a*q0*[X^q0]B` at degree q0-1. No smaller
input monomial exists to cancel it with a high term. Also
`q0-1<=M-1=j-r<j`, so it cannot be the target coefficient. This is a
contradiction. Therefore q0=0. The surviving residue class is consequently
zero, and

    j+1=r(ell+1),       M=r*ell,       ell>=0.

Using j=k-2 and 2k=A+r, this is exactly (1). Both positivity and the modulus
are necessary: merely satisfying a congruence with a negative ell would not
provide a polynomial primitive.

Conversely assume (1). Set

    B(X)=sum_(h=0)^ell c_h X^(r*h),
    c_ell=-2/[r*b*(2ell+1)],
    c_(h-1)=-2*a*h*c_h/[b*(2h-1)],   h=ell,...,1.        (9)

The highest coefficient in (6) is -1. At every lower potential degree
rh-1, the two contributions are

    r*[a*h*c_h+b*(h-1/2)*c_(h-1)],

which cancel by (9). There is no bottom derivative term because h=0.
Thus L(B)=-X^j exactly. In the original coordinates a primitive is

    y*R(u),
    R(u)=u^-k sum_(h=0)^ell c_h u^(-r*h),
    k=r*(ell+1)+1.                                    (10)

One may check (10) directly by `N R'+N'R/2=1`; this independent coordinate
check is included in the source. At ell=0 it is precisely the sharp family
primitive `-2y/(r*b*u^(r+1))` from genus_capacity.

All coefficients in (9) are rational in a,b, with only nonzero integer and
b-power denominators. Formulas (4) and (10) are therefore rational in the
original coefficient field and y. If a,b are nonzero rational functions of
parameters over any characteristic-zero field, the same formulas provide
rational parameter descent; no square root of a parameter is introduced.
The necessary condition can also be checked after extending the constant
field to an algebraic closure, so the exact formula does not hide a
constant-field exception.

For completeness, allowing a further finite algebraic extension does not
weaken this exactness criterion. In characteristic zero the u-derivation
extends uniquely to such an extension L/K and commutes with trace. If
`dC=eta` in L, then `Tr_(L/K)(C)/[L:K]` has derivative eta and lies in K.
This proves the equivalence of the first two statements without imposing an
extra field restriction on the meaning of an algebraic primitive.

## 6. Genus, whole pencils and the mate boundary

For the exact staircase, A and r have the same parity and n is even. There
are r simple residual branch points and an additional branch at zero just
when r is odd. The genus is consequently `floor((r-1)/2)`. Equivalently it
is the genus of the fixed model `Y^2=a+bX^r` in (5). At fixed r,a,b, changing
ell multiplies N by a square power of u. The curve remains the same but the
exact differential changes its weight. Formula (9) increases the odd
primitive space in precisely the compatible steps.

The primitive in (10) has poles only over the high root u=0. If r is odd
there is one such point with pole order A-2; if r is even there are two,
each of order A/2-1. In either case the total pole degree is

    A-2=r(2ell+1).

For arbitrary inputs the finite branch count is r+(A mod2), and infinity
adds a branch exactly when A+r is odd, proving (3). In particular, for even
r and odd A the genus is one larger than the staircase value. The explicit
A=1,r=2 control catches a missing infinity branch.

The nonzero-coefficient generic members of (2) are now exact. At either
coefficient-zero endpoint they are pure powers `N=c*u^q` with q nonnegative
and q!=2. Their separate primitive is

    sqrt(N)*[2*u^(1-q)/(c*(2-q))].                     (11)

For the linear pencil the constant endpoint is included by q=0. This proves
the all-member statement, without claiming that a formula with 1/b remains
regular when b=0. Through degree eight the exact generic pairs are precisely

    (A,r)=(0,1),(3,1),(5,1),(7,1),(4,2),(5,3),

matching the previously proved six spaces. Large-degree sparse exact pencils
are supplied by (2); a complete large-degree pencil classification is not
claimed.

There is also an exact, limited Jacobian consumer. For the pure quadratic
row `F=N(u)t^2`, a rational mate exists exactly for the binomials in the
theorem: necessity is the proved leading_exactness gate, while the primitive
`yR` gives the rational mate

    G=-R(u)/(2t),       J(F,G)=N R'+N'R/2=1.             (12)

The linear exception uses R=2/a. For a polynomial with this same leading
coefficient and additional lower t-terms, only the necessary exactness gate
is paid. No claim of sufficiency for those lower terms, global regularity
of G, or a Keller source is made.

## 7. Hostiles and precise failure boundaries

* **Algebraic versus logarithmic.** For N=u^2+1, the elementary expression
  `log(u+y)` differentiates to du/y. Nevertheless the two infinity residues
  are nonzero, so there is no primitive in any finite algebraic extension.
  The first failed implication is replacing algebraic exactness by elementary
  integration; the classical logarithmic term is the missing datum.
* **Same field, different weight.** The three polynomials
  `u^5(u^3+1)`, `u^7(u^3+1)`, and `u^11(u^3+1)` define the same elliptic
  function field. The first and third are exact, the middle one is not.
  At A=7,r=3, equation (6) forces B to have degree one; its nonzero linear
  coefficient creates an uncancellable constant term. Genus and field
  isomorphism alone do not preserve the weighted differential.
* **Low exponents.** A=0 is exact only at r=1. At r=1 the other exact values
  are odd A>=3. At r=2 the exact values are positive multiples of four.
  These checks retain the degree-one exception, infinity logarithms, and
  the first possible pole-bearing primitives.
* **Nonzero coefficients.** Dropping b!=0 changes the root inventory: for
  example A=0,r=3,b=0 gives the exact pure cubic although that exponent pair
  is rejected in the theorem. This is a degeneration outside its domain,
  not a counterexample. Square factors in u^A within the stated domain are
  retained by (5) and cannot be dropped from eta.

These controls identify the strongest survivor: a complete sparse operator
criterion for algebraic primitives, with explicit coefficient descent. They
do not turn a necessary leading-field gate into a sufficiency theorem for an
arbitrary two-variable polynomial.

## 8. Reproduction and finite scope

Run from the worktree root:

    python3 04-computation/planar_jc48_sep08_binomial_exactness.py
    python3 -O 04-computation/planar_jc48_sep08_binomial_exactness.py

The declared rectangle is every A=0..48 and r=1..12, 588 rows. Each even-degree
row at least four is checked against the full unrestricted polynomial
coefficient system (6), with two nonzero rational choices of a,b: 584 systems
in all. Their ranks are computed by exact rational Gaussian elimination, with
no restriction to the predicted exponent residue class. Low and odd degrees
retain the independent residue/pole-capacity tests. The rectangle has 74
accepted rows; no inference beyond the analytic theorem is made from that
finite range.

The source separately differentiates all 48 symbolic primitives r=1..8,
ell=0..5 in both the inverted and original coordinates. It checks coefficient
rationality, branch count and primitive pole degree, all monomial operator
images q=0..12 for r=1..8, the six low-degree pencil rows, the coefficient
endpoints, and the named logarithmic and same-field hostiles. Always-active
exact checks remain enabled under Python optimization. The source imports no
other producer.

Normal and optimized runs both pass **1,654 gates**, with byte-identical
628-byte outputs. Source and output are frozen for independent audit; the
primary is now promoted after the full independent audit.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `04-computation/planar_jc48_sep08_binomial_exactness.py` | 6562 | `ebc2df4e64b7abf29f9a1caf939cbe29329c94dfc7361590d0e61c50873d15bb` |
| `05-knowledge/results/planar_jc48_sep08_binomial_exactness.out` | 628 | `e70e6daab742bfa00fc7cc26daaacc9fb62f20df90395c01a99fc30550618683` |

The semantic rectangle hash is
`1f94751f29386df2e09b9058d23d90adff53ee3efad97e700d4aa0080af615c7`;
the gate-label hash is
`31c4523a4ff4b3ea29844a22a7c1ed19fe5b99fa54b64593fdcab3da7fe4927a`.
No prior frozen artifact has been edited.

The [independent complete audit](planar_jc48_sep08_binomial_exactness_audit.md)
accepts the full iff, parameter descent, boundaries and both1654-gate
replays. A separate original-coordinate Laurent calculation recovers
all1170 exact systems, including odd-degree failures. Classical
integration is credited; no external-priority claim is made.
