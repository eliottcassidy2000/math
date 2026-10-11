# Renewal laws, source floors, and odd-zeta refinement tails

2026-10-11. **PROVED:** fixed-depth natural-density renewal, density zero for
fixed-cap finite-bank families, rational atomic source floors with odd-zeta
tail constants, and the elementary integer-separation tests
below. **PROVED with CITED Matveev input:** every positive expense moment of
the Haar block law is finite. **FINITE-EXACT:** the companion controls.
**OPEN:** a source-specific renewal budget proving universal Collatz.
The renewal and density-zero conclusions repair and independently derive
claims of the [October5 minimum-map note](collatz_minimum_map_renewal_20261005.md);
they are not claimed as new discoveries or literature-priority results.

## Inheritance and what changed after the pull

Anchor: hard excursions at running minima. Niche: the integer/nonvanishing
gate in the [odd-zeta audit](odd_zeta_parallels_petersen_20261001.md).
Wildcard: shape and branch data in the
[planar triangle-carrier theorem](planar_triangle_carrier_20261011.md).
Board: **native clock / ordered carry / renewal mass / arithmetic separation /
source atom / marked terminal**. Every connection below is tested against
these retained coordinates.

The closest proved mechanism is
[THM-4512, coefficient-descent cylinders](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
The hostile is the integer1 in the contracting word(2); contraction alone
does not strictly descend. The corrected near miss is replacing an actual
odd excursion by a shortcut-map coefficient exit. The underused sidecars
are a cofinal ordinary-height threshold and uniform tail control before
passing from finitely many cylinders to a countable block law.

The October8–9 integration adds a stronger warning:
[source incidence and named terminals](collatz_fundamental_integration_20261009.md)
survive neither covariance quotients nor ambient Haar coalescence. Its
state-selected block technique motivates adaptive budgets here, but its
noncoalescence theorems are not imported as Collatz descent statements.
META-PATTERNS used: **Search the statement before the method**, **Separate
local support from bounded-height coverage**, and **Type every analogy**.

## R1. Exact renewal in the right clock

Use the accelerated odd map

\[
 U(n)=(3n+1)/2^{a(n)},\qquad a(n)=v_2(3n+1),\qquad \alpha=\log_2 3.
\]

For a positive word w=(a1,...,a_l), put A_j=sum_{i<=j}a_i, A=A_l,
P=3^l, Q=2^A, and B=sum_{j=0}^{l-1}3^{l-1-j}2^{A_j}.
The exact native cylinder is

\[
 n\equiv (Q-B)P^{-1}\pmod {2Q},\qquad U^l(n)=(Pn+B)/Q.       \tag{1}
\]

Its density **among odd integers** is1/Q. These finite cylinder identities
prove that the valuation digits under odd2-adic Haar measure are independent
with probability2^{-a}; no assumption about random integer orbits is needed.

Let C be the first coefficient-descent code:

\[
 2^{A_j}<3^j\quad(1\le j<l),\qquad 2^A>3^l.               \tag{2}
\]

It is prefix-free in valuation digits. If sigma is its stopping length,
survival to time t implies A_t<t alpha. Markov's inequality gives

\[
 \Pr(\sigma>t)\le \kappa^t,\qquad
 \kappa=\frac{\sqrt3}{2\sqrt2-1}<1,                       \tag{3}
\]

because E[2^{-a/2}]=1/(2sqrt2-1). Thus C has total mass1 and its successive
blocks under Haar measure are iid, with word mass2^{-A(w)}. The code-length
tail also has all polynomial moments. Endpoint3-adic information does not
alter this dyadic statement: multiplication by3 is a2-adic unit.
For a fully rational implementation, kappa<19/20 follows by squaring the
positive inequality1444sqrt2<2049; use(19/20)^t as a looser bound.

**Fixed-depth integer renewal theorem.** For every fixed k and every subset
E of C^k, the natural density among positive odds whose first k *actual*
running-minimum excursions belong to E exists and equals

\[
                   \sum_{(w_1,...,w_k)\in E}2^{-\sum_i A(w_i)}.       \tag{4}
\]

Sources with fewer than k such excursions are left outside E. Their density
is zero; this says nothing about any particular source.

**Proof.** Fix a tuple. Concatenation occupies one cylinder (1). All proper
prefixes of each block increase its own starting value, since their slopes
exceed1 and their carries are positive. The endpoint decreases precisely
when that block starts above B/(Q-P). Every block start is an increasing
affine function of the original source with positive slope. Hence, beyond
a finite ordinary-height threshold, this tuple gives exactly the first k
actual excursions, and all its states exceed any prescribed finite bank.
Conversely an actual occurrence has that native cylinder. This proves (4)
for a singleton and for finite unions.

For the countable union, choose a finite dictionary D subset C of mass at
least1-epsilon. The finitely many tuples D^k have total mass at least
(1-epsilon)^k and satisfy the claim cofinally. The complement has upper
density at most1-(1-epsilon)^k. Let epsilon tend to0. This supplies the
tightness missing from a bare assertion that each word has finitely many
exceptions. It also covers words that are actual first descents but not
first coefficient descents: their aggregate density at any fixed block
index is zero. No universal equality of the two stopping times is used. QED.

**Clock and conditioning controls.** Source13 has U-word(3), endpoint5.
The shortcut map instead goes13->20->10 and has coefficient exit after two
steps, at an even endpoint. These are different stopping objects. Moreover,
conditioning a block on fitting a bit budget changes its distribution:
among one-letter code words with a<=3, a=2 has conditional probability2/3,
not its unconditional probability1/4. On a finite residue space, each
specified tuple whose total A plus the final odd bit fits has its exact
cylinder frequency; this is not unconditional independence after selection.

## R2. Fixed caps and finite banks really have density zero

For a contracting word define

\[
 q(w)=\left\lceil\frac{l}{A-l\alpha}\right\rceil.
\]

Fix Q0 and a finite bank B0. Consider sources whose actual running-minimum
chain reaches B0 and whose preceding complete excursion blocks all satisfy
q<=Q0. This includes F_q and the usual fixed two-tier families: put all
small allowed starting states into the finite bank.

There is a code word with q>Q0 and positive Haar mass. Choose l>Q0, set
A_j=floor(j alpha) for j<l and A_l=ceil(l alpha). All successive differences
are positive; the terminal excess is in(0,1), so q>l>Q0. Let p>0 be the
mass of all such hard words. By R1 the density of sources whose first k
actual blocks avoid them is(1-p)^k. Reaching a fixed bank before those k
blocks has density zero, by the same finite-dictionary proof. Therefore
the bank family's upper density is at most(1-p)^k for every k, and is0.

This proves the qualitative conclusion without a moving-k error estimate.
The earlier explicit integer power-law exponent is **not established** by
this argument. A weak law of large numbers supplies no exponential error
rate, and an exceptional later minimum can have many source preimages.
Neither a one-exception-per-word bound nor a fitted exponent pays those
two bills. A killed iid renewal is a useful separate model; transferring
its exponent to integer starting intervals remains an additional problem.

The minimum map is partial until existence of first descent is proved.
On the ROOT basin its graph is a tree. Globally it is an acyclic partial
functional graph; saying that every chain ends at1 would already assert
Collatz. Every unit does have infinitely many incoming minimum-map edges:
take its arbitrarily large one-step inverse predecessors. This local
branching does not make all vertices belong to the ROOT tree.

## R3. Unbounded hard expenses nevertheless have every Haar moment

There is a useful positive alternative to a fixed cap. Put C0=10^10.
The explicit Matveev estimate already used in
[BC1](collatz_boundary_compiler_20261004.md#2-the-sharp-all-length-carry-bound)
gives, for A*=ceil(l alpha),

\[
 2^{A_*}/3^l-1>\exp[-C_0(1+\log(2l))].                 \tag{5}
\]

The imported statement is the positive-real algebraic-number estimate in
[Jeong–Park, Theorem2, p.172](https://www.mathos.unios.hr/mc/index.php/mc/article/view/5358/1593),
which attributes it to Matveev, Corollary2.3. Here the degree is1, there
are two logarithms, heights may be1 and2, and B=2l. The resulting constant
is below C0. Multiplicative independence makes the form nonzero. This is
a quantitative theorem about2 and3, not a consequence of odd-zeta
irrationality.

Writing delta*=A*-l alpha in(0,1), the inequality
2^{delta*}-1<=2 delta* implies

\[
 q(w)\le1+K_0 l^{C_0+1},\qquad K_0=2\,6^{C_0},             \tag{6}
\]

for every contracting word of length l, since a larger A decreases q and
e<3. The integer constant is explicit but far too large for a practical
controller; the immediate gain is the qualitative all-moment theorem.
Together (3) and (6) show

\[
                    \mathbb E_{\rm Haar\ block}[q^s]<\infty
                    \quad\text{for every real }s>0.       \tag{7}
\]

They also give a stretched-exponential upper tail with the explicit, very
weak exponent1/(C0+1). No square-root exponent or sharp constant is claimed.
This supplies a proved rate-versus-cost mechanism where qualitative
irrationality alone supplies only nonvanishing.

For example choose T_i=ceil(log((i+2)^2/epsilon)/(-log kappa)), i>=1.
The union bound shows that, with Haar probability at least1-epsilon, every
block satisfies l_i<=T_i and hence q_i<=1+K0 T_i^{C0+1}. A budget growing
as a power of log(i+2) can therefore cover arbitrarily high Haar mass,
whereas every fixed budget eventually fails almost surely. This is a
statement about the complete iid code, **not** a ROOT certificate or an
integer coverage theorem. Fixed-depth density convergence also does not
imply convergence of unbounded expense moments over finite integer ranges;
that would require uniform integrability under those integer laws.
Alternatively, choose the least integer T_i with
(19/20)^T_i<=epsilon/(i+2)^2. This conservative version needs only rational
arithmetic and gives a reproducible adaptive schedule.

## R4. An arithmetic escape certificate and the odd-zeta comparison

The hard block(1,1,3) acts by G(n)=(27n+19)/32, with fixed point19/5.
For k actual repetitions,

\[
                 5G^k(n)-19=(27/32)^k(5n-19).             \tag{8}
\]

If all these are positive odd integers above the fixed point, their left
side is an integer at least6. Thus

\[
                 (27/32)^k(5n-19)\ge6.                  \tag{9}
\]

This is a source-specific finite limit on consecutive repetitions. Every
finite repeated word has positive integer realizations, but its nested
2-adic cylinders have unique limit19/5, not a positive integer. In fact
their residues satisfy5r_k=19 modulo2^(5k+1). Compactness preserves the
2-adic address while losing the ordinary-integer domain.

More generally, for a finite collection of actual affine maps sharing a
rational fixed point a/b and slopes in(0,r] with r<1, a nonfixed integer
orbit staying on one side has |b n_j-a|>=1 and
|b n_j-a|<=r^j|b n_0-a|. This gives a finite exit bound. Different anchors
require their translation cocycle; replacing them by a common scalar
contraction is not valid.
An integer fixed point is excluded from this nonvanishing argument: a
cycle word can have zero residual and must retain its named cycle target.

The source of the zeta analogy is the integer-separation gate, not a
numerical identity between zeta and Collatz. Apéry-style arguments clear
denominators and force a nonzero integer to have absolute value below1.
The [repo's odd-zeta audit](odd_zeta_parallels_petersen_20261001.md) found a
proposed gain insufficient after denominator costs. Here (8) identifies
the exact integer first, and (9) pays the entire contraction. It proves
escape from this repeat grammar, not termination of the subsequent orbit.
The missing next step is a compatible separation quantity across switches.

## R5. What the shape analogy and the five papers contribute

The papers are sources of constructions, not assumed Collatz theorems.
These contracts supplement the October5 probes rather than rerun them.

| Source and mechanism | Proposed target and map | Preserved / lost / decisive test |
|---|---|---|
| [Hoelscher](https://arxiv.org/html/2102.07083v1): subset-sum completeness and product generating functions | Compress accumulated excursion costs into sums | Preserves cost; loses word order and native addresses. Words1132 and2113 have equal l,A but distinct carries and residues. A semicomplete cost set alone is not a guard cover. |
| [Liu–Wang–Zhang](https://arxiv.org/html/2608.08703v1): quadratic transforms and Hankel/Somos identities | Look for a smaller recurrence for the block law | Preserve the exact generating function and specialization conditions. The inherited finite Hankel probe rejects its fitted constant Somos-4 parameters; it does not prove absence of every algebraic or finite-order recurrence. Keep zero denominators and finite-head scope. |
| [Fox–Schildkraut](https://arxiv.org/html/2608.20533v1): a dependency graph permits a local-lemma avoidance construction despite bounded gaps | Design adaptive bad-block avoidance events | A positive probability of some avoiding path is preserved; a specified integer is lost. R3 provides an elementary growing-budget version, with the domain explicitly Haar. No AP-free-to-Collatz map is supplied. |
| [Riordan–Scott](https://arxiv.org/html/2609.15893v1): root swaps and an order bijection convert extendable frontier edges into copies | Count extendable certificate ports, with their source | Preserve incidence, not just degree. Their expectation averages the root; a dense component plus an isolated designated root defeats a source-specific inference. No artificial tournament is needed. |
| [Toth](https://arxiv.org/html/1801.08014v1): nested prime intervals with independently guaranteed children | Compare nested legal Collatz cylinders with their limits | Preserve compatible prefixes; lose the carrier of the limit. Repeated113 gives the explicit19/5 hostile. The [earlier Mills packet](seam_mills_20260925.md) retains its fractional carry for exactly this reason. |

The simplex/spheroid/cuboid picture is useful as three certificate domains:
a simplex records mixtures and extremal choices, an ellipsoid records a
quadratic budget, and a box records simultaneous coordinate bounds. They
are not interchangeable without an inclusion and its constants. For
example ||x||_infinity<=||x||_2<=sqrt(d)||x||_infinity; the upper loss is
attained at the all-ones vector. A Euclidean ball is polar self-dual; a
simplex is combinatorially self-dual; in dimension3 a cube is dual to an
octahedron. A circle and triangle can be homeomorphic as boundaries while
having different face and curvature data. There is no intrinsic universal
ordering making a circle the most complex shape.
The field matters too: the complex vector(1,i) has quadratic sum0 but
Hermitian norm squared2. Real ellipsoid positivity cannot be silently
inserted into the complex polynomial Jacobian equations.

The planar companion makes one triangle connection exact: its three marked
branch points on a sphere determine genus and a rational-flow obstruction.
This is different from rounding a polygon. In both lanes the correct
small object is a **marked carrier**, with its branch/guard and target.

## R6. Odd zeta values occur as exact source-refinement tail constants

There is a concrete connection stronger than the earlier irrationality
analogy. For an integer r>=1 put a probability on the positive odd source
2m+1, m>=0, by

\[
 \nu_r(m)=\frac1{(m+1)^r}-\frac1{(m+2)^r}.                \tag{10}
\]

Every atom is positive and rational, and the total is1 by telescoping.
This gives **universal coverage of source addresses**, without assigning
ROOT status to any unresolved address.

Fix m and a dyadic modulus q>m in the index coordinate. The cylinder
C_q(m) consists of indices m+kq, k>=0. Write its off-source mass as
E_q(m)=sum_{k>=1}nu_r(m+kq). Then

\[
 \boxed{\quad E_q(m)\le\frac{r+1}{q^{r+1}},\qquad
          \lim_{q\to\infty}q^{r+1}E_q(m)=r\zeta(r+1).\quad}       \tag{11}
\]

**Proof.** The fundamental theorem of calculus gives

\[
 \nu_r(m+kq)=r\int_{m+1}^{m+2}(kq+t)^{-r-1}\,dt.
\]

Since k>=1, summing upper and lower bounds on this integral yields

\[
 \frac{r\zeta(r+1)}{(q+m+2)^{r+1}}
 \le E_q(m)\le\frac{r\zeta(r+1)}{q^{r+1}}.
\]

The integral test gives zeta(r+1)<=1+1/r, proving the rational bound;
squeezing proves the limit. No irrationality claim is used. The condition
q>m removes smaller congruent indices and is essential. QED.

For r=2 the rational prior is

\[
 \nu_2(m)=\frac{2m+3}{(m+1)^2(m+2)^2},\qquad
 E_q(m)\le3/q^3,\qquad q^3E_q(m)\longrightarrow2\zeta(3).
\]

More generally even r exposes the odd zeta value zeta(r+1). Thus odd-zeta
constants arise here from summing the remote arithmetic-progression aliases
of a rational, telescoping source measure. They do not require irrational
weights at individual sources. One additional binary address bit doubles
q and divides the **proved tail bound** by8 when r=2. This is a precise
three-bit scaling statement, rather than an identification of unrelated
three-part patterns.

The designated source survives every refinement with conditional floor

\[
 \nu_r(\{m\}\mid C_q(m))
 \ge\frac{\nu_r(m)}{\nu_r(m)+(r+1)q^{-r-1}}
 \longrightarrow1.                                    \tag{12}
\]

This is a source-specific, refinement-stable lower bound in an explicit
atomic measure. It is **not** the ROOT-supported weight W or lambda of
the earlier papers. Those weights can vanish at an unresolved source;
nu_r deliberately remembers the source before convergence is settled.

Let R_q(m) be the nu_r mass in C_q(m) of sources that actually reach1.
An independent proof of

\[
                      R_q(m)>(r+1)/q^{r+1}              \tag{13}
\]

would force the designated source to reach1, because all the other atoms
together have at most that mass. The remaining obligation is the lower
bound (13), not (12). In particular a finite bank omitting m cannot pass
this test merely by summing its certified neighbours. A global identity
or paid source-aware argument is still needed. This is a weaker universally
covered carrier on which to strengthen the proof mechanism, with the
missing predicate explicit.

There is also a quantitative version for any score0<=V_j<=B. Let
M_q=sum_{j in C_q(m)}nu_r(j)V_j. A truthful lower measurement L_q<=M_q
gives the refinement-stable source bound

\[
 V_m\ge\max\left(0,
       \frac{L_q-B(r+1)q^{-r-1}}{\nu_r(m)}\right).         \tag{14}
\]

With exact M_q, the unclipped expression differs from V_m by at most
B(r+1)/(nu_r(m)q^{r+1}), uniformly in the entire score field. For the
[inherited nonroot score W<=1/3](collatz_source_refinement_floor_20261005.md),
take r=2 and q>m>=1: the subtracted bill is simply1/q^3. This is a
source-specific lower-bound compiler with an explicit error at every
refinement. It does not manufacture L_q; supplying a sufficiently large
independently justified measurement remains the substantive Collatz target.

## Next productive targets and reproduction

1. Combine (7) with a proved source-weighted, refinement-stable estimate.
   The r=2 atomic law in R6 gives a rational source floor and an8-fold
   improvement of its alias bound per extra bit. Seek an independent
   ROOT-mass lower bound (13); Haar expense moments do not supply it.
2. Search for a switch-compatible integer separation functional extending
   (8). Retain anchor changes as signed affine data before scalarizing.
3. For growing block depth, pay both consumed dyadic precision and the
   preimage multiplicity of low exceptional minima. Only then compare an
   integer coverage exponent with the killed-renewal exponent.

The [checker](../../04-computation/experiments/collatz_renewal_source_20261011.py)
and [exact output](collatz_renewal_source_20261011.json) use a complete
bounded word universe and independent literal replay, including clock,
ROOT, commutation, conditioning and compactness hostiles. They test the
formulas; the quantified results rest on the proofs above.
The normal and optimized runs agree byte for byte, with2,485 always-active
exact gates, including source-score reconstruction and the missing-target
hostile. No independent peer audit is claimed.

```text
python3 -B 04-computation/experiments/collatz_renewal_source_20261011.py
python3 -B -O 04-computation/experiments/collatz_renewal_source_20261011.py
```
