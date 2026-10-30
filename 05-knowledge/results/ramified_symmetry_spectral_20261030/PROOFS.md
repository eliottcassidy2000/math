# Ramified symmetry, spectral budgets and actual tournament operations

2026-10-30. **PROVED scoped deductions from the named inherited/CITED
inputs; FINITE-EXACT controls separate. Same-agent proof review.**
No literature-priority claim. Throughout K is the REAL skew adjacency
matrix: K_xy=1 for x->y, -1 for y->x, and zero on the diagonal. Norm
means its Euclidean operator norm, not the norm of the 0/1 adjacency
matrix or a metagraph transition operator.

## S1. The inherited ramified source and the operation boundary

Fix an odd prime p and an integer a with v_p(a^2+4)=1. Necessarily
p=1 mod4 and p>=5. Use the literal ring and clock

    O=Z_p[mu], mu^2=a mu+1, delta=2mu-a,
    R_k=O/(p^k), N=p^(2k), g=mu^4.

The inherited [LOGFREE_PRIME_TOURNAMENTS, Sections1--7](../session110_overnight_20261030/LOGFREE_PRIME_TOURNAMENTS.md)
proves that some recursively balanced g-invariant Cayley tournaments
T_k have

    sqrt(pN) <= ||K_k|| <= 200 sqrt(pN).                (1)

The lower bound holds for EVERY bank state; the selected balanced
states are prime for k>=2. Its proofs retain simple ramification,
the full child-sum constraints, new-conductor orbit periods, and the
cited Lovett--Meka discrepancy input. This packet does not re-prove
that construction or present an implementation of its optimizer.

Every source is regular and contains the actual affine group

    A_k={x -> g^t x+b : t mod p^k, b in R_k},
    |A_k|=p^(3k).                                      (2)

This is a CERTIFIED SUBGROUP of the full graph group. Prime tournament
structure does not make (2) the full group: S3b below proves a universal
larger group, and Session110 retains additional scalar hostiles in its
other depth conventions. The first-level product also has extra symmetries.
The new attachment choices below are made AFTER selecting a fixed
spectrally controlled source; the bank signs themselves are not replaced
by independent random bits.

## S2. A complete affine cycle law at every depth

For 1<=m<=k there are exactly two cycle types among the elements of
A_k of order p^m. Write their multiplicities L_m and U_m:

    L_m=(p-1)p^(3m-2),
    U_m=(p-1)(p^2+1)p^(3m-3).                          (3)

The U_m elements have no fixed points and every cycle has length p^m,
so they have N/p^m cycles. The L_m elements have fixed points. Put
r=2k-2m+1. Their cycle inventory is

    length1: p^r;
    length p^j, 1<=j<m: (p^(r+2j)-p^(r+2j-2))/p^j;
    length p^m: (N-p^(2k-1))/p^m.                     (4)

In particular their TOTAL number of cycles is

    2N/p^m-N/p^(2m).                                  (5)

Together with the identity these types exhaust A_k. Both the group
element order and its fixed-point branch are needed; order alone
does not specify its action on cuts or joint named ports.

**Proof.** The inherited valuation identity is
v_delta(g^t-1)=1+2v_p(t), until it vanishes in R_k. Take
f(x)=g^t x+b. For t nonzero mod p^k, set r=v_delta(g^t-1), and
write s=v_delta(b), with the zero residue assigned depth2k.

If s>=r, then b is in the image of g^t-1, so a translation conjugates
f to multiplication by g^t. Its order is p^m for
m=k-v_p(t), and its p^j-th power fixes
p^min(2k,r+2j) vertices. Successive fixed-point differences prove(4).
The number of such t is (p-1)p^(m-1); the number of image translations
is p^(2k-r)=p^(2m-1). Their product gives L_m. Summing(4) telescopes
to(5).

If s<r, the p^j-th power has linear part g^(t p^j) and translation
b(1+g^t+...+g^((p^j-1)t)). The geometric sum has delta-valuation2j,
by the same binomial valuation identity. Until s+2j>=2k, the translation
is outside the image of the linear part minus1, so this power has no
fixed points. At the first j with s+2j>=2k both parts become identity.
Thus every cycle has length p^ceil((2k-s)/2). The t=0 translation case
has the same conclusion directly.

There are exactly p^(3m) elements whose p^m-th power is identity:
t is divisible by p^(k-m), giving p^m choices, and b lies in
p^(k-m)O/(p^k), giving p^(2m). Hence exactly
p^(3m)-p^(3m-3) have order p^m. Subtract L_m to obtain U_m. QED.

### A short recursion for every finite named palette

Let Phi_k(Q)=sum_(f in A_k) Q^(number of cycles of f). Then the number
of Q-coloring orbits under THIS group is B_k(Q)=Phi_k(Q)/p^(3k), and

    Phi_0(Q)=Q,
    Phi_k(Q)=Phi_(k-1)(Q^(p^2))
             +L_k Q^(2p^k-1)+U_k Q^(p^k).             (6)

Indeed the old m<k type multiplicities do not change and each cycle
COUNT scales by p^2; the two order-p^k types supply the final terms.
For r fixed named binary ports use Q=2^r. This counts their complete
joint bit words, without introducing independent palette gauges.

The number of prime one-vertex attachment orbits under A_k is
B_k(2)-4 when k>=2: the two constant cuts and the two clone orbits
are precisely the forbidden types. This is a FRAME count. It becomes
the unmarked graph count only after a separate full transporter quotient.
S3b proves that A_k is always a proper subgroup here. Formula(6) is not A000568.

## S3. Symmetry left by a random attachment: two different regimes

For uniform Q colors on the vertices, S2 gives the exact expectation

    E|Stab_(A_k)| = 1 + sum_(m=1..k)
      [L_m Q^(-N(1-p^-m)^2)+U_m Q^(-N(1-p^-m))].      (7)

This uses E|Stab|=sum_f P(f fixes the coloring), which is Burnside
with the labelled measure retained. Every finite nonuniform IID
palette pi has the exact analogue: replace Q^(cycles-N) by
product_(cycle lengths j) c_j^(number of j-cycles), where
c_j=sum pi_i^j. This uses the FULL cycle inventory(4), not merely(5).

The least nonidentity permutation index is

    D=N(1-1/p)^2,                                     (8)

attained by precisely p(p-1) elements. All other nonidentity elements
have at most N/p cycles: for a fixed-point type with m>=2,
2/p^m-1/p^(2m)<1/p, and the remaining types are uniform cycles.
Consequently, with E=E|Stab_(A_k)|,

    0 <= E-1-p(p-1)Q^-D
         <= [p^(3k)-1-p(p-1)] Q^(-N(1-1/p)).          (9)

In particular even ONE uniformly random binary port destroys this
specified arithmetic subgroup with probability tending to1 exponentially
fast in N. More generally if max pi_i<=a_0<1, then

    P(Stab_(A_k)!=1) <= (p^(3k)-1) a_0^D,
    0 <= E|Stab_(A_k)|^s-1
       <= (p^(3ks)-1)(p^(3k)-1)a_0^D                 (10)

for every fixed s>0. The first bound uses c_j<=a_0^(j-1); the second
bounds the random variable on its nontrivial event. Thus every fixed
positive moment tends to1, without needing a growing palette.
Since p^(3k)=N^(3/2) and D>=16N/25, the bound is uniform over all
admissible p>=5,a,k as N grows, for fixed a_0<1 and s.

This differs sharply from the previous [recursive wreath law G3](../prime_power_plus_one_symmetry_20261030/PROOFS.md):
there independent bottom blocks leave many small-support symmetries,
while here the arithmetic group's minimum index is a fixed positive
fraction of N and its size is polynomial in N. These statements concern
specified permutation actions, not an inference from the spectral norm.

**Compulsory full-group hostile.** At p=5,k=1, A_1 has order125 and

    E|Stab_(A_1)|=1+20/2^16+104/2^20.

But the all-positive base bank is R5[R5], whose full group has order15625.
Its binary rigid-cut probability is1554*15625/2^25, about0.724.
The affine bound therefore cannot be advertised as the probability of
FULL graph rigidity. The next theorem makes this guard compulsory for
every bank in the present even-valuation-depth tower, including prime banks.

## S3b. Universal socle shifts and the actual full625 group

Let rho:R_k -> F_p be reduction modulo delta, and put

    eta=g^(p^(k-1))-1.

Then eta is nonzero of delta-valuation2k-1, eta^2=0, p eta=0,
and eta*x=eta*rho(x). For EVERY function h:F_p -> F_p define

    sigma_h(x)=x+eta*h(rho(x)).                        (10a)

These permutations form F_k isomorphic to the additive group C_p^p.
They preserve EVERY g-invariant Cayley tournament bank, regardless
of balance or the particular choices of orbit signs.

**Proof.** The residue is unchanged by sigma_h, so composition adds
the functions and inversion negates them. For two vertices with equal
residue their difference is unchanged. Otherwise their difference v
is a unit and its changed difference is v+eta*d for some d in F_p.
Choose s=d/rho(v). Square-zero binomial expansion gives

    g^(s p^(k-1))*v=(1+s eta)*v=v+eta*d.

Thus the ordered difference stays in its actual clock orbit and the
arc is unchanged. This proves preservation of the literal relation,
not just its spectrum or a quotient statistic. QED.

The affine group A_k normalizes F_k. Translations shift the argument
of h; the clock fixes both rho and eta and therefore commutes with F_k.
Moreover F_k intersect A_k consists precisely of h(z)=u+v*z, so

    H_k=F_k A_k,  |H_k|=p^(3k+p-2).                 (10b)

For the intersection claim, evaluate an affine map at0 to obtain
b in eta F_p. Its linear part minus1 must have image in eta F_p,
forcing t divisible by p^(k-1); the converse is immediate. Thus
H_k is a certified subgroup of every bank and A_k has index p^(p-2)
in it. These independent shifts do not require tournament modules:
they remain present in the prime balanced banks at k>=2.

There is an exact smaller-support symmetry law. Put

    D_0=N(p-1)/p^2.

A nonzero shift on one residue fibre moves N/p vertices in p-cycles,
so its permutation index is D_0. For uniform IID Q-colors,

    E|Stab_(F_k)|=(1+(p-1)Q^(-D_0))^p,
    P(Stab_(F_k)=1)=(1-Q^(-D_0))^p.                 (10c)

The p fibres have disjoint coordinates. On each fibre every nonidentity
shift has the same cycle partition because p is prime; its coloring
is fixed with probability Q^(-D_0). This proves both formulas.
In particular(10c) gives a lower bound on the mean residual FULL group
order and an upper bound on the probability of FULL rigidity.

In fact the least nonidentity index in H_k is exactly D_0, attained
by its p(p-1) single-fibre nonzero shifts. To check that no element
outside F_k competes, write it as g^t*x+b+eta*h(rho(x)). If t is
divisible by p^(k-1), it is in F_k when b is in eta F_p, and otherwise
has no fixed points. For other t, r=v_delta(g^t-1)<=2k-3. The fixed
equation, reduced modulo delta^(r+1), pins the residue of x; inside
that fibre it has at most p^r<=N/p^3 solutions. Since H_k is a p-group,
every moved cycle has length at least p. Its index is therefore at
least (N-N/p^3)(1-1/p)>D_0. The branch with other t is absent at k=1.
Two active fibre shifts already have index at least2D_0.
The IID event/moment bounds in(10) consequently hold for H_k after
replacing p^(3k) by p^(3k+p-2) and D by D_0. In particular all fixed
positive moments tend to1 at fixed p as k grows. This larger action
has a genuinely different leading exponent from A_k.

### A corrected all-depth cycle recurrence for the enlarged action

For k>=2 let Psi_k(Q)=sum_(f in H_k) Q^(number of cycles of f),
and let Phi_k be the S2 recurrence. Then exactly

    I_k(Q)=Q^N+p(p-1)Q^(N-(p-1)D_0)+(p-1)Q^(N-p D_0),
    F_k(Q)=[Q^(N/p)+(p-1)Q^(N/p^2)]^p,
    Psi_k(Q)=p^(p-2)[Phi_k(Q)-I_k(Q)]+F_k(Q).         (10e)

Here I_k is the cycle sum of A_k intersect F_k, and F_k(Q) the cycle
sum of the fibre group (the function and group are distinguished by
the argument Q). Thus S2 plus(10e) is an all-depth closed formula,
without enumerating H_k. Divide by p^(3k+p-2) for its Q-color orbits.

**Proof of the coset step.** Every element has form
f(x)=a(x)+eta*h(rho(x)), where a(x)=g^t*x+b. Put beta=rho(b).
Since the clock fixes eta and rho, iteration gives

    f^n(x)=a^n(x)+eta*sum_(j=0..n-1) h(rho(x)+j beta).

If beta=0, then f^(p^j)=a^(p^j) for every j>=1. For a outside F_k,
the fixed-point argument above also proves that f and a have the same
number of fixed points: when r<=2k-3, solve (g^t-1)y=-eta*h(z)
inside delta R_k in each residue z, giving a bijection of their fixed
sets; the remaining outside-F_k cases have no fixed points.

If beta is nonzero, a has all cycles of length p^k. For j>=2 the
extra sum vanishes in F_p, so f^(p^j)=a^(p^j). At n=p the extra
translation has delta-valuation at least2k-1>=3, whereas the
translation of a^p has valuation2 and its linear part minus1 has
valuation at least3. Hence f^p has no fixed points; f itself also
has none because it moves the residue. Thus f also has all cycles
of length p^k. This is where k>=2 is used.

The counts of fixed points of all p-power iterates determine the full
cycle inventory in a p-group. Consequently every nonidentity F_k
coset has a constant cycle type, and its A_k intersection has p^2
elements while the whole coset has p^p. Replace the identity coset's
I_k by its independent fibre sum F_k(Q) to obtain(10e). QED.

At k=1 the assertion is false: on the25-point ring the permutation
x -> x+1+eta*1_(rho(x)=0) is a single25-cycle, while the corresponding
translation has five5-cycles. The checker retains this minimal-depth
hostile. The socle-shift theorem(10a--c) itself still holds at k=1.

### FINITE-EXACT full group of the pinned625 source

For the particular balanced prime source from lifts.json, complete
individualization/refinement and literal arc checks give

    Aut(T_625)=H_2,  |Aut(T_625)|=5^9=1953125.        (10d)

Here vertices are (x,y) mod25, indexed25x+y. With ell=2x+y mod5
and d=(5,15), the full pointwise stabilizer of vertices0 and1 is

    (x,y) -> (x,y)+d*h(ell),  h(0)=h(1)=0,

an explicit C5^3 of order125. This d generates the same socle line
as eta. All125 maps are retained and separately checked against every
unordered pair in the actual tournament.

**Completeness certificate.** Rooted equitable refinement from0 has
five singleton cells,24 cells of size5 and20 of size25. The cell of
vertex1 is exactly its known25-element clock orbit. Thus the FULL
root stabilizer has orbit25 on1. After fixing1, canonical refinement
has250 singletons and75 size-five cells. Branching on every vertex
of the first smallest nonsingleton cell gives profiles

    (singletons, size-five cells; number of nodes):
    (250,75;1), (375,50;5), (500,25;25), (625,0;125).

Every leaf map relative to the first discrete leaf preserves all arcs;
no leaf is rejected. Conversely every automorphism fixing0 and1 maps
the first leaf's successive selections to one of these branches,
because the ordered refinement and cell-selection rules are invariant.
It therefore occurs among these125 distinct maps. Translation
transitivity and orbit-stabilizer give625*25*125=5^9. Formula(10b)
already supplies a subgroup of that size, so equality follows.
This completes ONE full graph group, not a census of all625 bank states.
For this source(10e) is therefore the FULL coloring-orbit formula:

    Psi_2(Q)=125[Phi_2(Q)-Q^625-20Q^225-4Q^125]
              +[Q^125+4Q^25]^5.

Its exact unmarked prime-child count is Psi_2(2)/1953125-4, by S4's
regular-parent theorem. This is a literal625-to626 branch of tournament
augmentation, not a new full A000568 term. The Q=2,3,4 orbit counts
and exact prime-child count are retained in parent_group.json.

At this prime source D_0=100 whereas A_2 has D=400. Even a highly
quasirandom prime source can therefore retain local shifts whose
color-stabilizer correction is exponentially larger than the affine
calculation predicts. The explicit626 child in S4 nevertheless destroys
this ENTIRE larger group. The literal socle operation, rather than a
coincidence of group orders, is the bridge to the earlier wreath theme.

## S4. A sharp universal budget for adjoining vertices

For ANY real skew matrix K and real vector b, put

    K'=[[K,b],[-b^T,0]],  M=||K||, B=||b||^2.

Then

    M^2 <= ||K'||^2 <= M^2+B.                          (11)

The lower bound is compression/interlacing. For the upper bound, consider
a real vector(x,t). Since x and Kx are orthogonal and ||Kx||<=M||x||,
Bessel's inequality gives

    (b^T x)^2+(b^T Kx)^2/M^2 <= B||x||^2.

Also 2t b^T Kx<=M^2 t^2+(b^T Kx)^2/M^2. Expanding
||K'(x,t)||^2 and applying these bounds gives
(M^2+B)(||x||^2+t^2). If M=0 the assertion is immediate. Equality
is possible, for example for the2-by2 skew block and b=(1,1).

For a tournament on N vertices, EVERY one-vertex cut has B=N. Iterating,
an arbitrary attachment of r vertices, with arbitrary mutual arcs, obeys

    ||K_extended||^2 <= ||K||^2+rN+binom(r,2).          (12)

This sharp per-step estimate complements the exact bordered characteristic
polynomial in [THM-1920, spectral insertion response](../../../01-canon/theorems/THM-1920-the-spectral-insertion-response.md).
It is uniform in the whole incidence pattern: selecting a cut to enforce
rigidity or a realization guard does not introduce a probabilistic
conditioning error into(11). The spectrum itself still depends on the
off-diagonal cofactor data, so a norm bound is not a complete observer.

### Rigid PRIME children without losing the spectral scale

Let T be a regular prime tournament whose every vertex stabilizer in
Aut(T) is nontrivial. There are at least two NONISOMORPHIC rigid prime
one-vertex extensions U of T, and EVERY such extension satisfies(11).

**Proof.** The CITED [Gluck Corollary1](https://doi.org/10.4153/CJM-1983-005-2)
gives a subset S with trivial stabilizer in any odd-order permutation
group; tournament automorphism groups have odd order because an
involution would reverse an arc. This known tournament application is
recovered in [THM-4743, converse-extended Gluck](../../../01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md).
The [regular-parent R1](../regular_parent_symmetry_20261029/PROOFS.md)
then makes the unmarked extension rigid.

The complete prime-parent forbidden list is inherited from
[THM-4169, Section1](../../../01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md):
only the two uniform cuts and the2N clone cuts fail to be prime.
A uniform cut retains Aut(T). A clone cut at v retains Aut(T)_v.
Thus none is a Gluck subset under the stated hypothesis, proving primeness.
S^c is also regular. No odd-order element sends S to S^c: squaring
such an element fixes S and regularity forces an involution. The two
cuts are in different Aut(T) orbits. Prime extensions have no clone
root collision, so the regular-parent classification makes their
unmarked tournaments nonisomorphic. QED.

Apply this to the ramified source at k>=2. Its point stabilizers contain
conjugates of <g>, of order p^k>1. Therefore there are rigid prime
tournaments on N+1=p^(2k)+1 with

    sqrt(pN) <= ||K_U|| <= sqrt((40000p+1)N).          (13)

The arithmetic lower scale persists by compression even though the
entire automorphism group has been destroyed. Thus that lower scale
cannot be interpreted as requiring a surviving symmetry of the final
object. Its arithmetic correlation is retained in the induced source.
For fixed p these rigid prime children still have ||K_U||/(N+1)->0.
For any disjoint vertex sets A,B their signed edge discrepancy is at
most ||K_U|| sqrt(|A||B|), by the matrix Cauchy--Schwarz inequality.
The cut may be chosen separately at each depth. Compatibility of the
distinguishing cuts under the ramified reduction is not asserted.

For any prime regular parent with t vertex orbits, its exact unmarked
prime-child count is (number of subset orbits)-2-2t. Indeed THM-4169
subtracts the source/sink and both clone orbits per vertex orbit, and
R1--R2 removes no further duplicates among prime children. The
transitive case is F-4. This combines the older prime-parent and newer
regular-parent results; it does not supply the unknown full group of
every ramified source.

### An orientation bit recovered from the root and parity

EVERY prime one-vertex extension of a regular tournament is nonself-converse,
whether or not it is rigid. Indeed its regular-deletion root is unique:
two such roots would be a clone pair by R1, contrary to primeness.
Any anti-isomorphism of the child to itself must preserve the set of
regular-deletion roots, and hence fix this root. But an anti-isomorphism
fixing a vertex equates its outdegree and indegree, impossible on an
even number of vertices. This contradiction proves the claim.

Thus no additional orientation selector is needed to distinguish these
children from their converses. It is the conjunction of an intrinsic
unique root and EVEN child order that recovers the sign. Dropping
primeness loses uniqueness: the clone children of C3,R5,P7 are explicit
self-converse hostiles. This is not a claim that every regular-parent
extension, or every re-marked/complemented target, has a free sign bit.

## S5. Actual Hamilton-path complements have bounded spectral cost

**INHERITED tariff, new composition/application.** The precise bound
below is already THM-4813(C),
[spectral exception budgets](../../../01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md),
with the earlier Session75 matrix/trace-norm work credited there.
[THM-4815(A1--A3), dipole principle](../../../01-canon/theorems/THM-4815-dipole-principle-native-move-distance-bounds-endpoint-relation-cut-or-colour-exactness-and-the-p7-parity-correction.md)
already organizes the same identity and potential method. The short
derivation is retained for the growth-word bound, not as a new tariff.

For a legal Hamilton path P of T let J_P be its skew adjacency matrix,
zero off the consecutive path edges. The ACTUAL complement retained in
the native port work is T^P=T^op xor E(P), so in the common named frame

    K_(T^P)=-K_T+2J_P.                                (14)

This is not merely global converse. After reordering along P and a
diagonal unitary change of basis, iJ_P is the symmetric path adjacency
matrix. Its sine eigenvectors give
||J_P||=2 cos(pi/(N+1)). Hence

    | ||K_(T^P)||-||K_T|| | <= 4 cos(pi/(N+1)) <4.     (15)

Re-marking is simultaneous permutation conjugation and costs zero.
After c legal complements at the same size, the norm changes by less
than4c. In particular o(N) complements preserve a starting norm o(N).
Primeness, regularity, group order and joint-port responses are NOT
asserted preserved by this spectral estimate.

For a growth word using r vertex insertions, c actual complements and
any re-markings, starting at size N and norm M, one has

    max(0,M-4c) <= M_final
       <= sqrt(M^2+rN+binom(r,2))+4c.                 (16)

To commute the scalar bounds, use that x->sqrt(x^2+n) is1-Lipschitz
on nonnegative x; each complement's additive cost may be moved to the
end of the bound. Insertion costs are the successive sizes N,...,N+r-1.
This is a compositional bound under the actual operations, not a
claim that these two scalar coordinates reconstruct their states.

The transitive tournament has norm cot(pi/(2N)), from its classical
explicit skew spectrum (also recovered in THM-1920's cited transitive
frame). Any fixed-size complement/re-marking route from a source of
norm M to a transitive tournament therefore needs at least

    [cot(pi/(2N))-M]_+ / [4 cos(pi/(N+1))]             (17)

complements, rounded up to an integer. With the source bound200 sqrt(pN)
and fixed p this is N/(2pi)-O(sqrt N). It is a lower bound on ANY
such route, not a claim that a route exists or a statement about a
different metagraph whose black moves change the literal tournament.

**Literal finite consequence.** The inherited625-vertex balanced prime
source has norm25 sqrt(5+2sqrt5), certified again by outward rational
Fourier bounds. Our explicit626-vertex rigid prime child therefore has
norm at most25(1+sqrt5), by(11). On the625-vertex source, the elementary
bound cot(pi/(2N))>=7N/11-11/(14N), using pi<22/7, together with the
certified parent upper norm, makes(17) exceed80 even with denominator4.
Thus ANY fixed-size native complement/re-marking route from this source
to the transitive tournament uses at least81 complements. The lower
bound does not presume connectivity or give a shortest path.

This81 is only a spectral control and is NOT an improvement on known
distance bounds. For any regular odd-order source,
[THM-4810, Theorem L](../../../01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md)
already gives the stronger lower bound(N^2-1)/8. Indeed sorted scores
define L=(1/2)sum_i |s_(i)-i|; converse preserves L and restoring a
Hamilton path changes only two scores by one, so each move changes L
by at most1. Regularity gives L=(N^2-1)/8 and a transitive target has
L=0. At625 the inherited bound is48828. The spectral coordinate is
useful for mixed insertions and irregular targets, where it composes
with(11), but not as a replacement for this stronger score obstruction.

## S6. What remains uncompressed

S2 retains whole permutation cycle types for joint named colors. S4
retains the actual cut when selecting a rigid child. S5 retains the
selected Hamilton path and its vertex transport. The spectral budget
cannot replace these data: the81/27 pair already separates core actions
behind equal cyclic triples; Session110 separates frame quotients from
full graph classes; and the previous weighted theorem separates an IID
palette law from an arithmetic realization law.

The all-size selected-complement integral recovery conjecture and the
smallest complement/re-marking carrier remain OPEN. The present new
operation bound supplies a genuine invariant constraint on their legal
paths. It does not turn the affine Burnside recurrence into a full
A000568 reconstruction.

The incoming [endpoint-shadow packet, G1--G5](../endpoint_shadows_path_corrections_20261030/TOURNAMENTS.md)
at de77c39fe5 supplies a decisive additional boundary: two legal moves
with the same endpoint ports can preserve all scores, have norm
correction below8, and still produce rank at least N-5sqrt(N)+6.
Its consistent path-program prime completions are not completions of
our fixed ramified source. Thus the small spectral budget cannot
license dropping full edge incidence, source marking or transporter
data. The two results compose as constraints, not interchangeable
representation theorems.
