# Exact ROOT-current folding and a joint contact obstruction

2026-10-17, clock_duality_oct15. **PROVED, scoped and self-audited:** the
odd-source kernel, algebraic degree obstruction, four-direction contact
classification, specialization/retraction implications, and the joint
linear-retraction obstruction over the full base field. **FINITE-EXACT:**
21 native ROOT-current snapshots, two complete induced contact graphs, all
570 boundary extensions, and the complete 888-parameter specialization atlas.
**OPEN:** uniform ROOT-current bounds, uniform small-uncertainty entry,
six-colour forcing, and planar JC. No independent referee or novelty claim.

## 1. Inheritance and live concepts

Anchor: an independent uniform bound on the already certified ROOT current.
Niche: mixed-coordinate forcing beyond the Cartesian spindle cube.
Wildcard: the prime463 as a marked conductor clock.
The board is **literal odd source / integer clock / signed covariance /
contact coefficient / fixed ports / compatible specialization**.

Closest mechanisms: [D1, refinement-preserving ROOT brackets](escape_renewal_idoneal_20261017d.md)
and [K1, the marked cubic Fourier formula](keller_marked_cubic_fourier_20261017d.md).
Canonical hostiles: the old spindle cube has no mixed-coordinate edges; a
nonempty relative six-core still allows every boundary pattern; norm and
Lipschitz contraction alone do not pay a uniform tail. Corrected near miss:
converse-merged nodes and free presentation orbits are different objects.
Least-used sidecars: the coefficient of the new quadratic generator and
the SAME specialization of a shared cubic rotation in several attachments.

The startup inherited aca9619bba. The later fetch through8817698205 brought
the S39 audit and [THM-4665's EI17 finite-field tests](../../01-canon/theorems/THM-4665-ei17-has-no-faithful-realization-over-the-split-and-conjugation-fixed-primes-of-q-zeta55-tested-so-every-q-zeta55-realization-has-bad-reduction-there.md).
Its circle-intersection elimination prompted the complete local test in
section8. Its faithful-embedding condition is different from our retractions:
nonadjacent points may collapse here, while ports remain fixed. The incoming
audit's corrected signs and quantifiers are retained; its restricted Collatz
mod-p equidistribution is not promoted to a source-uniform current estimate.
THM-4663/4664 enter only through the marked-coordinate mechanism, not as
colouring or Collatz theorems. Recovered stronger results include THM-781's all-order
Hamiltonian-path fibre inverse, HYP-6825's three-sorted incidence description,
and [G1's complete six-leaf boundary relation](golden_taxicab_boundary_20261011.md).
Method cards used: controlled forgetting with a retained sidecar; compute the
repair quotient before the residual defect; re-evaluate after specialization.
No new meta-pattern is promoted.

The coefficient method has earlier stronger precedents than a finite spindle
test: [B8, cubic adelic work](cubic_adelic_20261016e.md) already states the
two-coordinate degree gate, and [G1--G3, the cubic hub](anchor_deficits_cubic_hub_20261017a.md)
already retain the middle coefficient of a cubic-square norm and its
denominator19. Section4 supplies the general degree statement for reference;
the new geometric work is the quadratic contact, its complete retraction
atlas, and the incompatibility of three shared contact conditions.

## 2. Fold doubling fibres without losing their signs

Use the shortcut map T(n)=n/2 for even n and (3n+1)/2 for odd n. Let f be
any real profile with f(2n)=-f(n), and let w(n)=3[n=2 mod3]-1.
Write chi(n)=0,1,-1 for residues0,1,2 modulo3, and

    R_f(N)=sum_(N<n<=4N) w(n) f(n).

For the native current, f(n)=(-1)^clock(n) when the authenticated terminal
is stopped ROOT1, and f(n)=0 otherwise. The native graph starts with all
halving edges and closes retained paths, so this profile obeys the required
doubling identity on the whole finite source universe. Its clock is the
actual number of shortcut steps; it is not replaced by a word length or
an endpoint residue.

**R1 (PROVED).** For every N>=1,

\[
 R_f(N)=-3\sum_{\substack{m\le2N\\m\ \mathrm{odd}}}\chi(m)f(m)
       +\sum_{\substack{2N<m\le4N\\m\ \mathrm{odd}}}w(m)f(m). \tag{1}
\]

For each odd m<=2N, exactly two consecutive members of its doubling chain
lie in (N,4N]. Their combined coefficient is
(-1)^a(w(2^a m)-w(2^(a+1)m))=-3chi(m), independently of a.
Odd m in (2N,4N] contribute just their first member. This proves (1).
It is valid even when f vanishes at unresolved sources.

Thus each odd source has only two threshold events:

\[
 K_N(m)=w(m)[N\ge\lceil m/4\rceil]
       +(-3\chi(m)-w(m))[N\ge\lceil m/2\rceil].        \tag{2}
\]

For m=1,3,5 mod6 the two coefficients are respectively (-1,-2),(-1,1),(2,1).
This is an exact compression of all doubling addresses; the marked odd
source and its signed clock survive. In particular, putting R_f(0)=0,

\[
 \begin{split}
 R_f(N)-R_f(N-1)={}&w(4N-3)f(4N-3)+w(4N-1)f(4N-1)\\
 &+(-3\chi(2N-1)-w(2N-1))f(2N-1).                  \tag{3}
 \end{split}
\]

The three absolute coefficient values sum to4. When |f|<=1 this recovers
the inherited four-unit increment bound, rather than improving that bound.

**Exact energy kernel.** On X<=N<2X, set t_0(m)=ceil(m/4), t_1(m)=ceil(m/2)
and c_0(m)=w(m), c_1(m)=-3chi(m)-w(m). Then

\[
 G_X(m,n)=\sum_{i,j=0}^1c_i(m)c_j(n)
       (2X-\max(X,t_i(m),t_j(n)))_+,
 \quad \sum_{N=X}^{2X-1}R_f(N)^2=\sum_{m,n\ \mathrm{odd}}f(m)G_X(m,n)f(n). \tag{4}
\]

Only m,n<=8X are needed. The kernel counts overlaps of marked threshold
relations; it is PSD because it is K^T K. Removing signs, clocks or source
positions would change the contracted quadratic form. Formula (3) computes
all rows in linear time; it does not bound (4) for Collatz signs.

**Hostile.** f(n)=chi(n) has f(1)=1 and f(2n)=-f(n), but gives R_f(N)=-3N.
It violates the odd Collatz relation. Hence normalization, doubling
self-similarity and four-Lipschitz control do not imply the desired
X^(2+o(1)) energy bound. The odd transition is an essential remaining
arithmetic constraint, exactly where a free orthogonality argument would
otherwise have erased it.

## 3. What the native current actually does

The fixed first-reset policy is tested at X=2^10,2^14,2^18 and stages
1,2,4,8,12,20,28. Every snapshot has a complete observed-edge audit,
including exterior excursions, and retains source endpoint/clock hashes.
Write (1) as R=I+B. Selected normalized energies are:

| k | stage | ||I||²/X² | ||B||²/X² | <I,B>/X² | ||R||²/X² |
|---:|---:|---:|---:|---:|---:|
|10|4|2.475022|3.585965|-2.644469|0.772049|
|14|4|5.562330|4.045576|-4.301327|1.005252|
|18|4|0.378666|0.913585|0.113718|1.519686|
|18|28|0.516032|0.771572|0.173301|1.634206|

The sign change is an exact finite rejection of a universal nonpositive
I/B covariance rule for this policy. At k18 stage28, R ranges from162 to1070;
at k10 stage28 it ranges from-76 to54. The already certified current need
not itself decrease monotonically when ROOT sources are added.

The open sufficient target remains an independently proved uniform bound on
(4), for a declared controller policy, together with uniform entry into
||U||_2<=eta||R||_2 for some eta<1. The preceding refinement brackets then
bound the actual current energy. A finite matrix identity, a random-sign
model, and a census of eventual ROOT labels do not discharge this target.
The component correlations must be retained jointly; neither a fixed
negative sign nor independence survived these exact tests.

## 4. Why high degree can prevent contacts

Let K be a subfield of C closed under conjugation. Suppose |u|=1 and
[K(u):K]>2. For a,b in K,

\[
 |a+ub|^2=1\ \Longrightarrow\
       a=0,\ |b|=1\quad\hbox{or}\quad b=0,\ |a|=1.   \tag{5}
\]

Multiply the norm equation by u. It is a quadratic polynomial over K;
its leading and constant coefficients are b*conj(a) and a*conj(b).
The degree hypothesis forces both to vanish. Consequently a two-coordinate
sum of a finite seed over K is exactly its Cartesian product graph, with
injective addresses. This explains a whole class of failures of generic
high-degree rotations, beyond the previous finite spindle probe.

More generally, if [K(u):K]>2r, a polynomial sum sum_(j=0)^r a_j u^j
has norm1 only if exactly one coefficient is nonzero and that coefficient
has norm1. Clearing the Laurent norm polynomial gives degree<=2r; its
outermost coefficient is nonzero whenever two nonzero endpoints remain.
This statement needs the degree bound; the cubic-square experiment at r=2
is outside it and remains an independent finite check.

## 5. Build a mixed edge and classify every other one

Take d,e in K nonzero with |d|²!=|e|², retain u as above, and put

\[
 Z=d+ue,\quad N=Z\bar Z,\quad s^2=1-4/N,\quad
 v=\tfrac12 Z(-1+s).                                \tag{6}
\]

Assume 0<N<4 in the chosen complex embedding and s not in K(u), with
conjugation s->-s. Then |v|=|Z+v|=1. The intended mixed edge is built in.

**G1 (PROVED).** The map (a,b,c)->a+ub+vc from K^3 is injective. Its
unit vectors are exactly

\[
 (q,0,0),\ (0,q,0),\ (0,0,q),\ (dq,eq,q),\qquad q\bar q=1. \tag{7}
\]

**Proof.** Injectivity follows first from s not in K(u), then from u not
in K. If c=0 use (5). For c!=0 set A=a/c, B=b/c. The s-coefficient of
the norm equation vanishes only if (A+uB)/Z is real. Clearing u gives
coefficients of a polynomial of degree at most2. Its outer coefficients
imply A=d*conj(B)/conj(e). The middle coefficient then gives

    (|d|²-|e|²)(conj(B)/conj(e)-B/e)=0.

Thus B=te and A=td with real t in K. The remaining norm condition is
|c|²(1+N t(t-1))=1. Also N is not in K, since otherwise its expression
|d|²+|e|²+d*conj(e)/u+conj(d)*e*u would give another quadratic equation
for u. Hence t=0 or1 and |c|=1. This proves (7).

The retained coefficient is doing the same logical job as at+2b in the
marked cubic Fourier proof: a preliminary elimination leaves an exact
vanishing test, whose omission would change the target predicate.
The calculations are different; no Fourier bound is imported into geometry.

## 6. Two exact contact graphs, and why they still extend

Let K=Q(zeta15,sqrt(-11)), degree16; eta=(1+sqrt(-3))/2 and
rho=(5+sqrt(-11))/6. The spindle is

    M={0,1,eta,1+eta,rho,rho*eta,rho*(1+eta)}.

Use the old cubic w^3+2w²-4w-4=0, 1<w<2, and
u=(w+3-sqrt(-11))/(w+3+sqrt(-11)). It has degree3 over K:
w is recovered rationally from u, and a cubic cannot acquire a root over
the degree16 field K. In (6) take d=1+eta, with either e=1 (straight)
or e=-eta (tilted). Both choices satisfy 0<N<4. For h=w+3 in (4,5),

    N_straight=(7h²+11-2h sqrt33)/(h²+11),
    N_tilted=(h²+77-2h sqrt33)/(h²+11).

Positivity and N<4 follow directly from 4<h<5<sqrt33.
Both extensions have degree96. Nonsquareness of s² is certified at a
split31-place of K(w): zeta15=9,w=25, and j=19 or12. The corresponding
(N,s²) values are (23,17) and (16,24);17 and24 are nonsquares modulo31.

Let A be the54 distinct points in ports+M. Each finite graph is A+uM+vM,
with all eight original ports kept first. Formula (7) enumerates every edge.
A separate all-pair finite-field sieve followed by exact characteristic-zero
norm checks verifies the complete induced graph independently of (7).

| construction | vertices | axis edges | mixed edges | origin mixed edges |
|---|---:|---:|---:|---:|
| straight |2646|15470|65|8|
| tilted |2646|15470|65|8|

Each origin block has343 vertices and1625 edges, so the new eight contacts
are real. Nevertheless all570 original six-colour boundary patterns extend,
and both full graphs have independently replayed four-colourings.
The entire2638-vertex exterior is still a relative six-core.

There is now an exact explanation stronger than a SAT result. On the whole
K-linear space V=K+uK+vK, define

\[
 P_{straight}(a+ub+vc)=a-\eta b-\eta c,\qquad
 P_{tilted}(a+ub+vc)=a+\bar\eta b-c.                 \tag{8}
\]

Both maps fix K pointwise and send all four families (7) to unit vectors:
the fourth maps to 1-eta or eta-1 times q. Hence they are **unit-graph
retractions onto K**, preserving every boundary relation whose ports lie
in K. They are not field homomorphisms. Their finite images have825 vertices/
4328 edges and977 vertices/5474 edges respectively. Every one of the570
boundary patterns also extends on these images, giving a second set of
source colourings by pullback.

Changing e from1 to-eta also only rescales a coordinate in the full K^3
four-direction graph. It defeats one chosen additive-colour formula, but
does not change the full graph's isomorphism type. The port-fixing
retraction, rather than the failed formula, is the decisive obstruction.

## 7. Enumerate the repair possibilities before enlarging the graph

For any d,e in G1, a pair of base-field units r,t with

\[
                   |d+re+t|=1                              \tag{9}
\]

gives the retraction a+ub+vc -> a+rb+tc, fixing the original first axis.
Thus the full-space chromatic number and all boundary extension relations
at ports in K equal those of the base field whenever (9) is witnessed.
This is a relative equivalence, not a claim that K is six-colourable.

The exact catalog uses all34 nonzero spindle differences for d,e, unequal
squared norms, and the geometric feasibility condition
(|d|-|e|)²<=4. Its specialization bank is all14 oriented spindle unit
directions for r,t; the output d+re+t may be any exact base-field unit.
**All888 eligible ordered parameter pairs have a witnessed retraction.**
This closes the single-coupler strategy throughout this declared seed bank,
subject to G1's field and degree hypotheses, without building888 large graphs.

The atlas also retains *every allowable r* in the bank as a14-bit mask.
There are236 distinct masks and2222 incompatible pairs of mask
representatives. Of these representatives,248 pairs share the same d.
For example catalog rows369 and378 have masks1 and8192. They share

    d=rho*(1+eta)-(1+eta),
    e= rho*eta-1   and   -e.

Each single attachment specializes; no common r from the declared bank
specializes both. This is **FINITE-EXACT bank incompatibility**, not proof
that a common r in K does not exist. Section8 upgrades this test using the
whole norm-one residue group and three shared couplers. Separate positive
witnesses cannot be combined without a shared r.

The natural metagraph here has contact parameters as vertices and their
allowed marked rotations as incidence data. Disjointness is symmetric, so
it is a compatibility graph/hypergraph, not a tournament. Keeping the mask
makes the finite compatibility test exact; keeping one witness per vertex
would lose the obstruction.

## 8. A genuinely ternary obstruction over the full base field

Keep K,u,eta,rho from section6 and write D=1+eta. Choose the following
three displacement pairs from the atlas:

| atlas row | d | e | certified interval for N=|d+ue|² |
|---:|---|---|---|
|1| -D | -1 | (3.3798,3.3799) |
|162| rho-D | rho*eta-D | (2.5865,2.5867) |
|315| eta-rho | (rho-1)D | (1.5174,1.5175) |

The intervals are rational, from a rational isolating interval for
w^3+2w²-4w-4 on(1,2) and rational square-root bounds for3,11,33.
They certify0<N_i<4, not merely approximate numerical placement. Define
v_i by(6) for each pair and let

    V = K + uK + v_1 K + v_2 K + v_3 K.

This is an actual complex vector space inside an extension of degree at
most384 overQ. No independence of the v_i, exact compositum degree, or
classification of all joint unit vectors is assumed.

**G2 (PROVED).** There is no K-linear unit-graph homomorphism V->K fixing1
(and therefore fixingK and all original ports).

**Proof.** Such a map would give r=P(u) and t_i=P(v_i), with

\[
 r\bar r=t_i\bar t_i=1,\qquad
 (d_i+e_i r+t_i)\overline{(d_i+e_i r+t_i)}=1
 \quad(i=1,2,3).                                      \tag{10}
\]

These conditions follow from actual unit vectors u,v_i,d_i+e_i u+v_i
inV; they do not assume that these are its only unit vectors.
At the prime29, use the residue field F29[t]/(t²-2), with

    sqrt(-3) -> 10t,   sqrt(-11) -> 3t,   zeta15 -> 2+25t.

These satisfy the defining relations, with zeta15^5=(-1+10t)/2.
Frobenius sends every generator to its complex conjugate. Thus the chosen
place ofK is conjugation-fixed. Every norm-one element a ofK is a unit
there: v(a)+v(bar(a))=2v(a)=0. The coefficients d_i,e_i have denominator
dividing12 and are also integral there. Consequently every solution of(10)
must reduce to a solution over the complete norm-one group of size30.

Enumerating every residue r and every residue t_i gives the following
allowed-r masks (lexicographic order on(a,b), a²-2b²=1):

    625600743, 962028976, 503717383.

Their intersection is empty. The independent checker uses a different
integer encoding of all841 field elements and checks all2700 contact
equations, including every failed entry. This disproves(10) overK, hence
disproves the proposed linear graph homomorphism. Notice the integrality
argument: this is not a bounded-height search or a failure within the
fourteen-direction bank. There is no bad-denominator escape for r,t_i at
this conjugation-fixed place. QED.

**Pairwise positive controls overK.** Each pair of the three displayed
contact systems is simultaneously solvable overK, with the same r for that
pair. Put e_*=(rho-1)D, a unit. Explicit assignments are:

| pair | r | t for its first row | t for its second row |
|---|---|---|---|
|1,162| -1 | 1 | rho*eta |
|1,315| e_* | rho | 1 |
|162,315| -e_* | -e_* | rho*(1-eta) |

All nine assigned numbers have norm1; all six output expressions in(10)
also have norm1. They are independently checked in Q(sqrt(-3),sqrt(-11)).
The positive claim concerns the displayed equations, not a classification
of full two-coupler graphs or a proof that their every edge retracts.

Thus pairwise feasibility is genuinely insufficient for this marked
equation system. Keeping all allowable shared rotations detects the
ternary obstruction; keeping a single witness, or only binary adjacency,
does not. The no-retraction conclusion removes a proved obstruction to
the former strategy. It **does not** imply a missing six-colour boundary
pattern, absence of nonlinear graph maps, or a chromatic lower bound.
The next consumer is still a full joint graph and its port relation.

## 9. Binary information and eliminated witnesses

The tournament metagraph's [THM-781 path-fibre inverse](../../01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md)
recovers all presentations from Hamiltonian paths modulo automorphisms.
A merged node alone is a quotient; the observer fibre is the reconstruction
coordinate. The same discipline appears here in three different forms:

* the source Gram contraction retains clocks and literal odd-source positions;
* the contact graph retains which coordinates changed, enabling (7)--(9);
* the specialization atlas retains the shared-rotation mask, detecting the
  pairwise-solvable but jointly impossible system(10).

Binary constraints are complete when their hidden vertices are retained.
Eliminating them can create a higher-arity relation. The already proved
six-leaf star is the clean hostile: its boundary accepts every colour
partition except the all-six-distinct one. Every binary projection is full;
even every five-port projection is full. Replacing the six-ary relation by
pairwise compatibility loses the only obstruction. This is why the570-row
boundary relation remains the consumer, rather than a degree sequence or
an arbitrary tournament orientation.

The [463 conductor companion](prime463_conductor_clock_20261017e.md) adds an
all-order positive example of this research move: retain one ideal and its
coefficient valuation, then derive the entire class-group tower instead of
enumerating it. It supplies no current estimate or geometric obstruction.

## Receipts and remaining targets

[Manifest](../../04-computation/experiments/marked_contact_current_20261017e/manifest.json)
records commands and hashes. The ROOT controls exhaust729 small ternary-valued
odd-source profiles. Native runs authenticate all21 endpoint/clock snapshots;
normal and optimized outputs match. The two contact graphs have independent
all-pair exact norm checks, complete positive colour witnesses, and separate
combinatorial replays. The specialization atlas retains every source pair,
all14-bit rotation masks, and witnesses for each set bit.
The independent local replay checks the entire selected residue obstruction
and every characteristic-zero pair witness. Rational interval certificates
authenticate the actual cubic embedding; the residue place is used to
obstruct specialization intoK, not to specialize the cubic u itself.

Still OPEN: a uniform bound on the signed ROOT kernel (4), uniform entry into
the uncertainty regime, a joint coupler with an actual missing six-colour
boundary pattern, and planar JC. The next geometric comparison must use a
common specialization and actual cross-edges; the next current comparison
must use signed source correlations and the odd transition guard.
