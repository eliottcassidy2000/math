# Same-source modular observations: exact precision and sharp paid forcing

Session73; source commit `92e7f12f3b95840cff92c5e6cf40b2ae8952cae4`.
**PROVED scoped:** P1–P3 and W1–W3 below. **FINITE-EXACT:** `check.py` controls.
No native ROOT positivity, Collatz convergence, or all-word linear exception
budget is claimed. These are elementary arithmetic/orthogonal-projection
specializations, not new general CRT, least-squares, or electrical theorems.

## Result first

Two complementary repairs are necessary when moduli change:

1. **Keep denominator and intermediate guards.** For a fixed odd Collatz
   word of length h and total valuation S, endpoint odd-index precision q
   pulls back to the *minimal* source odd-index period
   `2^S q/gcd(3^h,q)`, with an image-divisibility guard. The physical clock
   remains h+S. This extends the inherited dyadic precision law to every q.
2. **Keep observations on the same source.** All separate q- and m-residue
   weighted current sums can lose exactly `lcm(q,m)-q-m+gcd(q,m)` joint
   directions, even when CRT makes the paired residue an exact address.
   With uniform mass on one complete residue frame this is also the exact
   unresolved squared evaluation norm. Thus positive modular marginals plus
   nearly constant energy need not force one source positive.

A supplied weighted Gram matrix gives the sharp interval for the source
value, and hence the strongest possible positive certificate from those
observations and that energy budget. Replacing the marginals by genuinely
same-source joint cells removes the cycle ambiguity; aliases beyond the
joint modulus still require the inherited weighted-renewal tail payment.

## Inheritance and scope audit

- `weighted_renewal_moving_modulus_20261026/README.md`, R1–R4: supplied weights,
  the Hardy cut budget, same-clock refinement, paid passive identity, and
  source-alias extraction. Their AUDIT says internal, not independent, review.
  R4 supplies the one-cell sufficient test; W1 identifies its exact optimal
  boundary and combines several nonorthogonal observations.
- `clock_even_guarded_composition_20261026/CLOCKS-MACHINES.md`, C2/C3: division
  consumes precision; branches and integer clocks must remain. P1 identifies
  the exact least period for a fully specified affine guarded word.
- `native_return_necklaces_20261026/README.md`, N6: the affine Collatz cylinder
  and clock h+S are inherited. P2 derives its arbitrary-modulus image/fibre law.
- `session72_even_graph_unifying_carriers_20261026/session72-response/REPORT.md`,
  Theorem B: pair-marginal kernels are graph cycle spaces. W2 applies that
  mechanism to CRT cells; the rank argument is inherited, not rediscovered.
- `collatz_source_refinement_floor_20261005.md`, SF4: a fixed evaluation norm
  is an inverse-Gram quadratic form, not an averaged spectral statistic.
  It already relates this to THM-4015 / character-sensitive Foster transference.
  W1 adds the sharp *data-conditioned feasible interval* and W3 computes the
  residue-frame threshold; no novelty for inverse-Gram evaluation is claimed.

Canonical hostiles: cancelled division with a surviving branch guard; two
incompatible CRT residues; the alternating rectangle on residues modulo 2
and 3; a positive supplied source mass mistaken for a proof of unknown ROOT
mass. Closest corrected near miss: R4's energy/alias test controls a signed
observable, but changing precision cannot supply its directional signal.
Least-used sidecar: the weighted same-source incidence matrix rather than
separate modular distributions.

Concept board: **source mass / guarded denominator / integer clock / joint
incidence / evaluation energy / signed forcing**. Anchor is paid pointwise
forcing; niche is exact affine precision; wildcard is electrical cycle debt.
This uses the repository meta-pattern “retain incidence when margins lose
targets,” not a new pattern promotion.

## P1. Exact least precision for affine words with congruence guards

Let an integer word have stages

    x_j=(a_j x_(j-1)+b_j)/d_j,    d_j>0,
    x_(j-1) ≡ r_j (mod h_j),     h_j>0,

and require every intermediate x_j to be integral. Retain output observations
`x_j mod m_j`, where m_j>=1; m_j=1 means integrality only. Fix one legal
integer source x*. Guards here are the stated *residue classes*, not arbitrary
Boolean unions of residues. Define without cancelling fractions

    A_0=1, B_0=0, D_0=1,
    A_j=a_j A_(j-1),
    B_j=a_j B_(j-1)+b_j D_(j-1),
    D_j=d_j D_(j-1).

Then precisely the sources `x*+Q Z` have all the same specified guard residues
and output residues, where

    Q = lcm_j {
       D_(j-1) h_j / gcd(A_(j-1),D_(j-1) h_j),
       D_j m_j / gcd(A_j,D_j m_j) } .                 (P1)

Use gcd(0,N)=N. This includes actual reset slopes a_j=0. It is a least
positive period on the integer domain; any separate positivity or height
guard must be retained and is not made periodic by this theorem.

**Proof.** Substitution gives x_j=(A_j x+B_j)/D_j. Relative to x*, the
jth output difference is A_j delta/D_j. It is integral and zero mod m_j
exactly when D_j m_j divides A_j delta. The incoming guard similarly requires
D_(j-1) h_j divides A_(j-1) delta. For integers A,N>0 in absolute divisor
convention, N|A delta iff N/gcd(A,N)|delta. Intersecting these subgroups of Z
is the lcm in (P1). These conditions also prove all intermediate integrality,
so necessity and sufficiency hold simultaneously. QED.

**Minimal hostile.** The legal word x -> x/2 -> 2(x/2), with even-input guard
and final residue modulo 3, has Q=6. Its cancelled endpoint is x, whose
unguarded period is 3. Sources 2 and 5 agree at the cancelled endpoint mod3,
but 5 fails the first guard. Final rational simplification may not erase an
intermediate denominator/guard. Zero-slope reset likewise erases no earlier
guard observation.

## P2. Arbitrary-modulus Collatz cylinders and exact composition

Write U(n)=(3n+1)/2^v2(3n+1) for the **odd-accelerated** map, not the
shortcut parity-step map called T in CLOCKS-MACHINES. Fix k_1,...,k_h>=1 and
S=sum k_j. The empty word is also allowed with
h=S=0, c=d=1. In inherited N6 notation,

    n=c+2^(S+1)t,   U^h(n)=d+2*3^h t,   t>=0,

where c in [0,2^(S+1)) is the unique odd cylinder representative and d is
odd positive. Set z=(n-1)/2, w=(U^h(n)-1)/2, c0=(c-1)/2, d0=(d-1)/2.
Then

    z=c0+2^S t,       w=d0+3^h t.                    (P2)

For q>=1 let g=gcd(3^h,q). The endpoint requirement w≡r modq is nonempty
**iff** g divides r-d0. In that case t occupies one residue modulo q/g,
obtained by inverting 3^h/g modulo q/g, and the legal source index occupies
one residue modulo

    Q_word(q)=2^S q/gcd(3^h,q).                       (P3)

For q/g=1 the inverse calculation is vacuous and every t qualifies. Thus
(P3) is the exact least source congruence period preserving the fixed word
and endpoint residue, not merely a sufficient overestimate. In the odd-source
n coordinate the corresponding period is twice (P3). Intersect with
t>=0 for actual positive sources. The physical clock is still h+S.

**Proof.** Equation (P2) is N6's exact cylinder, now in both odd-index
coordinates. Solve the linear congruence 3^h t≡r-d0 modq by dividing its gcd.
The remaining multiplier is a unit, giving one class modq/g. Distinct branch
sources differ by 2^S times an integer, so the displayed period is minimal.
Every residue class of t has a nonnegative representative. QED.

**Composition.** Split a word into u,v, with (h_u,S_u),(h_v,S_v). First pull
q back through v, then pull that precision through u:

    Q_u(Q_v(q)) = Q_(uv)(q).                         (P4)

Indeed on the prime 3 the two gcd cancellations subtract h_v and then h_u
from v_3(q), truncated at zero; on prime2 the S contributions add; every
other prime is unchanged. The actual affine residue and the image guard
must also be propagated. Equality of modulus sizes alone is not a valid
composition certificate. Native clocks add exactly:
(h_u+S_u)+(h_v+S_v)=h_(uv)+S_(uv).

**Original source weights.** If source z has weight μ_z, a target cell gets
`sum μ_z` over precisely its legal preimage. The g-to-one image law on a
complete t modq frame is a cardinality fact, not a uniform posterior law.
No factor g may replace an actual nonuniform weighted sum.

## P3. Joint cells use CRT compatibility, not an artificial product clock

For x≡a modq and x≡b modm, set g=gcd(q,m), L=lcm(q,m). Their intersection
is empty iff a and b differ modg. Otherwise it is one residue modL. This
follows by writing x=a+qt and solving qt≡b-a modm. For a list of such
requirements the pairwise gcd compatibility conditions are necessary and
sufficient, by repeated CRT (or prime-power nesting). Guards with several
allowed cells must retain their actual unions and intersections.

Under a source weight μ, overlap of two test indicators is the *actual*
μ-mass of this intersection. It is not the product of their masses except
where independence has separately been established. These precision changes
are observation refinement at the same physical timestamp. Combining f_t
and f_u as if they were one f requires a separate transport identity and
clock bookkeeping; this packet never does so.

## W1. Sharp pointwise forcing from a weighted observation Gram matrix

Let X be finite or countable, with supplied masses μ_x>0, and fix a named
source m. Let g_1,...,g_k be real L²(μ) tests on the *same* signed observable
f at one specified physical time. Indicators of guarded source cells are
allowed. Write

    y_i=<g_i,f>_μ,   C_ij=<g_i,g_j>_μ,   v_i=g_i(m),
    ||f||_μ² <= G.

For redundant tests use the Moore–Penrose inverse C^+; equivalently delete
dependencies and invert the resulting positive Gram matrix. Compatibility
requires y∈range C. Put

    E0=y^T C^+ y,
    a=v^T C^+ y,
    κ=1/μ_m-v^T C^+ v >=0.                            (W1)

Then the data are feasible iff y∈range C and E0<=G. For feasible data the
**entire possible interval** for the source value is

    f(m) ∈ [ a-sqrt(κ(G-E0)), a+sqrt(κ(G-E0)) ].       (W2)

If κ=0 this is the singleton {a}. If κ>0, both endpoints and every point
between are attained by actual real L²(μ) functions having exactly y and
energy at most G. Therefore these data force f(m)>0 iff

    a>0 and a²>κ(G-E0),                              (W3)

with the obvious κ=0 specialization a>0. Equality allows f(m)=0; it is not
strict positivity. No positivity/nonnegativity restriction on the unknown
f is assumed; adding such information can shrink the interval.

**Proof.** Evaluation is represented by h=1_{m}/μ_m, with norm²1/μ_m. Let V
be the span of the g_i and P its orthogonal projection. The least-norm
solution of the observation equations is f0=Σ_i(C^+y)_i g_i; every solution
is f0+u with u perpendicular to V. Thus ||f||²=E0+||u||². Also f0(m)=a and
||(I-P)h||²=κ. Cauchy–Schwarz gives (W2). If κ>0 choose
u=t (I-P)h, where |t|<=sqrt((G-E0)/κ); this realizes the entire interval.
If κ=0 all residuals have zero evaluation. Since μ_m>0, every null relation
among the g_i vanishes at m, so v is in range C and the redundant formula is
well defined. QED.

**One-cell reduction.** If a cell contains m with mass M=μ_m+E and observed
sum S, then a=S/M, E0=S²/M, κ=E/(μ_m M). For feasible data with E>0,
(W3) simplifies exactly to

    S>0 and S²>E G.                                  (W4)

This recovers weighted-renewal R4's test with the *actual* alias mass E,
and proves its minimax necessity from this information alone. Any failure
of (W4) permits an f(m)<=0 with the same S and budget G. Bounds E<=Ebar or
G<=Gbar give safe but potentially nonsharp sufficient tests; replacing E by
a convenient prior's alias tail without preserving μ is invalid.

**Important noncircularity boundary.** μ_m>0 is a hypothesis, not a conclusion.
This theorem can force a signed ROOT observable against a separately supplied
positive source measure. It cannot prove positivity of an unknown ROOT mass
by taking that unknown mass itself as μ_m and assuming it positive. SF4's
atom-positivity problem is logically different.

## W2. Separate modular current marginals have an exact CRT cycle kernel

For this section take the finite source model X=Z/LZ, q,m>=1,
g=gcd(q,m), L=lcm(q,m), with any positive supplied
weights μ_x. Observe *all* weighted current sums in residues modq and modm.
Represent x by an edge between (x modq) on the left and (x modm) on the right.
CRT identifies the resulting graph with g disjoint copies of
`K_(q/g,m/g)`, containing L edges and q+m vertices.

Writing z_x=μ_x f(x), the observations are the vertex-edge incidence map
Dz. Its exact rank and kernel dimension are

    rank D=q+m-g,
    dim ker D=L-q-m+g=g(q/g-1)(m/g-1).                (W5)

Call the lost-law dimension β=L-q-m+g; it is not generally the evaluation
residual κ in(W1). The weights do not alter these ranks because diag(μ) is
invertible. If any weights vanish, restrict to the actual supported-edge graph
and recompute its components; the complete-bipartite rank need not remain.
If the named source has zero mass, W1 is unavailable.  By the
inherited Session72 incidence argument, alternating cycle vectors form the
kernel; a spanning forest plus one exact current coordinate per chord repairs
it optimally. The full collection of *same-source joint* cell sums instead
gives each edge z_x exactly. Separate marginals cannot be multiplied to obtain
those joint sums.

**Proof.** CRT gives an edge exactly when the two residues agree modg, hence
the complete bipartite components stated. Multiply right rows of D by -1.
A row dependence of this oriented incidence matrix must be constant on each
component, so its row-nullity is g. Subtract from q+m and then from L to get
(W5). An alternating rectangle provides a kernel vector whenever both shore
sizes exceed1; spanning trees supply an independent fundamental-cycle basis.
This is a specialization of the inherited theorem, not a new rank mechanism.

**Minimal nonnested hostile.** For q=2,m=3 and uniform original masses1/6,
the signed f=(1,1,0,-1,-1,0) has every mod2 and mod3 weighted marginal zero,
although f(0)=1. Its paired address modulo6 is perfectly exact. The failed
implication is from exact *state* addressing to exact *linear-law* recovery.
Nested moduli have κ=0 in this complete finite frame because the finer marginal
already specifies each source. On an infinite source set, aliases beyond L
remain even after joint recovery and must still be controlled.

### Weighted electrical form, with the original masses retained

Let B be oriented incidence and W=diag(μ). The observation Gram matrix after
right-row sign change is the graph Laplacian `B W B^T`. For source edge e,
let b_e be its incidence column and define

    R_eff(e)=b_e^T (B W B^T)^+ b_e.

This is the standard effective resistance with conductance μ_e. Formula(W1)
becomes

    κ_e=1/μ_e-R_eff(e).                              (W6)

Equivalently its dimensionless unresolved fraction is
`1-μ_e R_eff(e)`. A bridge has κ_e=0; an edge on a cycle has κ_e>0.
For the latter statement, a nonzero alternating cycle current through e
has zero marginals and nonzero evaluation, hence (I-P)h is nonzero; for a
bridge, summing marginals on one side determines its current. This is a
direction-specific payment, not a trace or average-resistance guarantee.

## W3. Exact uniform-frame price of concealing one source

In W2 with the original uniform masses μ_x=1/L, fix source0. The orthogonal
projection of h=L1_{0} onto the marginal observation space is

    Ph = q 1_{x≡0 modq} + m 1_{x≡0 modm}
                        - g 1_{x≡0 modg}.           (W7)

Indeed subtracting this from h has zero sum in every q- and m-cell; within
its component each q-cell has L/q points and each m-cell L/m points.
Consequently

    κ=L-q-m+g.                                      (W8)

This equals the total missing-law rank for this symmetric frame, although
that equality is not a general weighted identity.

Each choice of L here is a separate finite illustrative model with masses1/L.
The family as L varies is **not** a refinement of one fixed original source
measure. At fixed measure, adding tests can only decrease the unresolved
evaluation norm; changing a convenient prior is not a lawful refinement.

Suppose every marginal is that of the constant observable c>0, meaning the
raw q-cell sums are c/q and the raw m-cell sums are c/m. Its
least-energy compatible solution is f0=c with E0=c². If κ>0, the minimum
energy among all compatible observables satisfying f(0)=0 is exactly

    G_zero=c²(1+1/κ),                               (W9)

attained by `f=c-(c/κ)(h-Ph)`. Above this threshold a compatible negative
source value is possible. Thus strict positivity is forced exactly when
G<G_zero (with feasibility G>=c²). When κ=0 the marginals already give
f(0)=c regardless of the remaining energy allowance.

For q=2,m=3, κ=2. The unit-root, zero-marginal minimum-energy residual is

    u=(1,1/2,-1/2,-1,-1/2,1/2),  ||u||²=1/2.

Hence f=1-u=(0,1/2,3/2,2,3/2,1/2) has every marginal of the constant1
but source value0 and energy3/2. For growing coprime q,m, κ=(q-1)(m-1)
and G_zero tends to1: even energy arbitrarily close to the constant minimum
cannot by itself certify that particular source. The explicit zero witness
is nonnegative, so this equality-boundary hostile also survives a mere
nonnegativity restriction. A negative witness necessarily leaves that cone.


## W4. A modulus bank has separate alias and joint-incidence payments

**Fixed original measure, PROVED.** Take any supplied finite positive atomic
measure μ on a finite/countable integer source set, with μ_m>0, and a finite nonempty
bank of moduli q_i. Let V be the span of the separate residue-cell indicators
and U the space of functions of their *joint* residue address. Then V⊆U.
Write C for the joint cell containing m, M=μ(C), and E=M-μ_m. Let P_V,P_U be
the weighted orthogonal projections. For h=1_m/μ_m,

    κ = ||h-P_V h||²
      = ||h-P_U h||² + ||P_U h-P_V h||²
      = (1/μ_m-1/M) + (1/M-v^T C_Gram^+ v).           (W10)

Both terms are nonnegative. The first is the source-alias debt
`E/(μ_m M)`; the second is the separate-marginal/joint-incidence debt.
The identity follows from the nested projections and `P_U h=1_C/M`.
It uses the exact support and original masses, not a replacement uniform law.
For guards, use their actual joint source partition; cells of zero mass are
ignored, but the named source must retain positive mass.

Adding the *single named joint-cell test* 1_C to V removes the second debt
for this source. Indeed h-1_C/M is perpendicular to every function of the
joint address, while 1_C/M belongs to the enlarged observation span. This
is a local pointwise repair; recovering the entire source law can require
many more joint coordinates. Actual data and residual energy still enter
(W2); no sign appears from the debt split alone.

**Finite fixed-N uniform formula, PROVED.** Keep X=Z/NZ and the same original
masses1/N throughout, with every q_i dividing N. Put L=lcm_i q_i. The
separate-marginal rank is

    R = Σ_(nonempty I) (-1)^(|I|+1) gcd(q_i : i∈I).   (W11)

Here repetitions in the bank are allowed; they cancel correctly. Its
unresolved source energy is

    κ=N-R=(N-L)+(L-R),                               (W12)

where N-L is alias debt and L-R is joint-incidence debt. Full joint residue
address has L classes; it determines an individual source exactly iff L=N.
If N=L, this extends W2/W3 to an arbitrary finite modulus bank. If the q_i
are pairwise coprime, R=1+Σ_i(q_i-1), although L=Π_i q_i.

**Proof.** Let P_q be uniform conditional expectation onto q-residue cells.
Direct CRT counting gives `P_q P_r=P_gcd(q,r)`; equivalently their finite
Fourier supports are the characters whose periods divide q and r. These
orthogonal projections commute. The projection onto the sum of their ranges
is `P=I-Π_i(I-P_qi)`. Expand this product. Each nonempty term is the projection
P_gcd, of trace/rank gcd; trace gives(W11). Also
`P_q h=q 1_{x≡m modq}`, so evaluation at m, or translation invariance of the
projection diagonal, gives ||Ph||²=R. Formula(W12) follows from(W10), since
the joint cell has mass1/L. QED.

The sharp constant-marginal zero-root threshold remains
`G_zero=c²(1+1/(N-R))` whenever N>R. An exact extremizer is
`f=c-c(h-Ph)/(N-R)`. It is nonnegative: translation invariance makes every
diagonal of the residual projection equal, and Cauchy–Schwarz bounds
`|(h-Ph)(x)|<=N-R`. Thus this hostile does not rely on negative f. It is a
real-valued example, not automatically a Boolean ROOT indicator.

For the separate illustrative model N=30 and bank{2,3,5}, R=8 and κ=22;
constant1 marginals hide a zero source at energy23/22. On that *fixed* N,
adding the genuine joint modulus30 test gives R=30 and κ=0. One must not
compare these with masses1/6 at N=6 as though the source measure were unchanged.
For the fixed-N statement, adding any observation space grows P and therefore
can only decrease κ; it never worsens the available certificate.

## Composition contract and remaining obligations

The safe pipeline is:

1. Keep the original named source, weight, full native word and integer clock.
2. Pull each endpoint/guard query back with P1/P2, including image restrictions.
3. Intersect same-source queries with CRT compatibility. Sum original weights;
   do not pretend that separate marginal observations are intersection data.
4. Form the exact weighted Gram matrix for the observations actually available.
5. Supply an independent energy bound G for the actual signed observable and
   apply W2/W3's sharp feasible interval. Precision does not create G or sign.

Changing observational precision is not adding physical transitions. P4 pays
division and composes clocks exactly; W1 handles the actual overlapping tests.
Neither equation licenses dropping named shared ports in the tournament
carrier: its joint-incidence loss is analogous, not an identified native
operator. No estimate here bounds the tournament exception-graph size.

**New useful boundary:** distributing a precision budget among several
separate modular response tables can leave a large cycle kernel despite
exact CRT addressing. The repair is a joint test/sidecar, or a quantitative
κ-dependent energy budget. A claimed uniform forcing proof across these finite models must therefore
pay the joint-law loss as well as the changing source mass. For a genuinely
fixed original measure, the general Gram theorem and the two-debt identity
below apply; no loss is caused merely by adding observations.

Reproduction and exact finite universes are in `check.py` and its receipts.
`SOURCE_PINS.json` records commit/blob/SHA256 for recovered sources;
`MANIFEST.json`/`SHA256SUMS` seal this packet. No source status was promoted.


## Incoming correction at `64a7985e705`: polynomial degree aliases modulo p

This bounded correction concerns **HYP-9419(4)**,
`05-knowledge/hypotheses/HYP-9419-unifying-candidates-three-layer-composition-rank-laws-dissection-obstructions-pattern-clocks.md`,
and §9(4) of its source
`05-knowledge/results/three_quantities_fortyfour_five_dudeney_fermat_clock_arrows_20261026.md`.
Their proposed exclusion of every polynomial sequence of degree>=3 from every
tournament-pattern clock is **REFUTED**. No inference is made about the
unrelated unproved rank/depth/dissection clauses. The incoming status remains
an audit-pending source, not an independently certified package.

The explicit polynomial `x^4` has nonzero residues `{1,2,4}` modulo7, so its
circulant is exactly Paley7. Modulo3 the same polynomial has nonzero pattern
`{1}`, giving the directed triangle. Both lie below120, so the broad
“FINITE-EXACT to120” headline is not true for all degree>=3 polynomials.
The source note's §7 census names **cubes**; those tests cannot license the
quantifier over all polynomials.

**PROVED corrected monomial criterion.** For an odd prime p and integer d>=1,
the nonzero values of `x^d modp` give a tournament connection set iff

    gcd(d,p-1)=2 and p≡3 (mod4).                       (I1)

Indeed the multiplicative group F_p^* is cyclic. Its dth-power image has
index gcd(d,p-1), so half-size requires that gcd to be2. Then the image is
the quadratic-residue subgroup; it is disjoint from its negative iff -1 is
not a square, equivalently p≡3 mod4. These are necessary and sufficient for
the inherited skew-half-set criterion. This uses the standard finite-field
cyclic-group fact, not a new classification of general polynomial images.

Odd d has the stronger immediate hostile to orientation: `(-x)^d=-x^d`, so
the nonzero image contains1 and−1 for every modulus>=3. The cube census is
consistent with that elementary obstruction; even degrees behave differently.
At a fixed prime, positive exponents with the same residue mod(p−1) give
the same function on all residues, so high polynomial degree alone cannot
be a modular obstruction. Even more generally, for any fixed modulus m and
D>=3, `x²+m x^D` has degree D but exactly the square pattern modulo m.
These are observation aliases, not an identification of their integer dynamics.

**Additional composite-unit boundary.** The incoming source §7 describes
square-oriented composite moduli as those for which −1 is a nonsquare in the
unit group. That property is not sufficient for the full square image:
modulo15, −1=14 is not a unit square, but 6 and9 are both squares and are
opposites. Nonunit residues must remain in the connection set. This does not
contradict THM-4809 C's explicitly odd-prime square criterion.

Exact controls enumerate every prime p<=199 and d=1,...,60, comparing actual
power-value sets with (I1); they also check all odd degrees1,...,59 at every
modulus3,...,120 for the ±1 hostile and the degree-alias construction at the
Paley primes<=120. No Fibonacci/primes pattern claim is decided here.
