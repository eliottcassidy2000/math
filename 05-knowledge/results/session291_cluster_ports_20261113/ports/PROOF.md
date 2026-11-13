# Paired ports: an exact arithmetic observation code, not the Golay code

**Status: PROVED / FINITE-EXACT where stated; independently reviewed.** All claims below concern
Session290's declared 24 even charts, with its third free port fixed at length
one. The underlying chart construction and its unit-defect identities are
inherited; the arithmetic observer and comparisons here are new applications.
No full Chan–Pak theorem or physical addition operation is asserted.

## 1. Recover the actual 12-to-24 map

Write U_j(m), 0≤j<24, for the selected-start slope in
[Session290's exact profiles](../../session290_core_descent_20261113/descent/PROOF.md).
Its actual tail has H=1+2^t U_j(m). At m=1, charts 2i and 2i+1 have literally
the same adjacency rows **and selected mask**, but different growing-port
pairs. Forgetting the pair gives 12 marked bases. There are only eight
literal adjacency-row lists, so “12 base graphs” must not be read as a claim
of 12 different unmarked tournament isomorphism classes.

The paired family labels are inherited, not chosen to fit a lattice. Choosing
one port pair in each marked base gives ε∈{0,1}^12. Let

    F_ε(m)=Σ_i U_(2i+ε_i)(m),
    Δ_i(m)=U_(2i+1)(m)−U_(2i)(m).

F is a formal aggregate of 12 separate observations. Every summand is a
physical source slope, but no tournament whose H equals this sum is supplied.
This is precisely the lost-operation boundary of the analogy.

## 2. Exact profile lattice and scalar observability

All 24 profiles lie in the exponential basis

    B=(1,2,3,4,6,8,9,12),  U_j(m)=Σ_b M_(b,j)b^m.

The integer 8×24 matrix M is included in receipt.json, extracted directly
from the inherited profiles. Its rank is eight and its nonzero Smith factors
are

    1,1,1,2,2,2,18,216.

The 8×12 paired difference matrix D has columns M_(·,2i+1)−M_(·,2i), rank
seven, and Smith factors

    1,1,2,2,6,72,216.

Thus seven rational response modes retain all 12 binary choices, even though
arbitrary rational combinations have a five-dimensional kernel. That last
qualification matters: real linear dependence does not imply a collision
between binary selections.

**Theorem 1 (all-length separating aggregate).** All 4096 binary choices have
distinct entire profiles F_ε. At m=1 every F_ε is equal. At m=2 there are 4014
values: 3932 singleton fibers and 82 double fibers. For every m≥3, all 4096
values F_ε(m) are distinct.

**Proof.** Exact enumeration of the finite cube gives 4096 distinct vectors
Dε and the stated first two fibers. It also tests every integer 3≤m≤89 with
exact integer arithmetic. This is a finite certificate for the small regime,
not extrapolation. For different ε,η, put c=D(ε−η). Its entries are integers,
not all zero, and Σ|c_b|≤L=Σ_(b,i)|D_(b,i)|=36288. Let B* be its largest
nonzero base. Every smaller positive base in B is at most (8/9)B*. For m≥90,

    |Σ_(b<B*) c_b(b/B*)^m| ≤36288(8/9)^m <1≤|c_(B*)|.

The exact integer check 36288·8^90<9^90 proves nonvanishing forever. Together
with the finite early check this establishes the all-order result. □

This theorem is an observer-selection tool, not a new path-count coverage
family. It gives a short reproducible separation certificate rather than
interpreting 24 as an ambient lattice dimension by itself.

## 3. Exact gcd and p-adic births of port blindness

Set x=2^m, y=3^m. The eight bases are the monomials

    (1,x,y,x²,xy,x³,y²,x²y).

Every difference polynomial vanishes at (x,y)=(2,3), because its two ports
share the same terminal graph. Let L0⊂Z^8 be the coefficient lattice of
polynomials in this span vanishing there. Evaluation is the primitive row
(1,2,3,4,6,8,9,12), so L0 is saturated of rank seven. The columns of D generate
a sublattice L⊂L0 with

    [L0:L]=373248=2^9·3^6.

The script checks the Smith form, and independently takes the gcd of all
792 maximal minors after deleting the constant row. This projection identifies
L0 with Z^7, proving the index independently of the Smith routine.

**Theorem 2 (exact scalar coupling divisor).** For every m≥2,

    gcd_(0≤i<12) Δ_i(m) = gcd(2^m−2, 3^m−3) =: G(m).

**Proof.** Every integer polynomial vanishing at (2,3) lies in the integral
ideal (x−2,y−3), so G divides every Δ_i. Conversely the finite index implies
2^9·3^6 L0⊂L. In particular 2^9·3^6(x−2) and 2^9·3^6(y−3) are integer
combinations of the difference polynomials. Thus their value-gcd divides
2^9·3^6G; only excess factors 2 or 3 are possible.

The i=7 difference is

    Δ_7=−485x²+336xy+531x−336y−130.

For m≥2, x=0 mod4 and y is odd, hence Δ_7=2 mod4. Therefore the common
2-adic valuation is exactly one, also v2(G).

The i=5 difference is

    Δ_5=−810x³+540x²y−150x²+926xy+389x−516y²−356y−22.

For even m, x=1 mod3 and y=0 mod3; Δ_5=1 mod3, so there is no factor3.
For odd m≥3, x∈{2,5,8} mod9 and y=0 mod9. The pairs
(Δ_5,Δ_7) mod9 for these three x values are respectively (3,0),(0,3),(6,6).
Therefore at least one is divisible by3 but not9. This equals v3(G)=1.
All prime valuations now agree. □

**Corollary (local blindness schedule).** For p≥5 let

    d_p=lcm(ord_p(2),ord_p(3)).

Then p divides all the paired differences precisely when d_p divides m−1.
When m−1=d_pq>0,

    v_p G(m)=min(v_p(2^d_p−1),v_p(3^d_p−1))+v_p(q).

This follows directly by lifting the exponent, since p divides neither base.
Together with v2G=1 and v3G equal to1 for odd m and0 for even m, it specifies
all prime-power losses of this observer. For example G(3)=6, G(5)=30,
G(7)=42. This is exact recurrence-controlled p-adic blindness, not a numerical
association between 12 ports and particular primes. Even when all differences
are divisible by a growing prime-power, Theorem1 distinguishes the actual
integers. Modular equality and integer equality play different roles.

A gcd of G does not mean binary subset sums fill any interval in GZ. Theorem1
is injectivity, not surjectivity, and no semigroup closure is granted.

There is also an exact finite-alphabet comparison: the coefficient vectors
Dε have respectively4,81,1480,3056,4096 images modulo2,3,5,7,13. Thus the
mod13 coefficient profile still identifies all twelve binary port choices,
whereas the mod2 profile loses all but two bits. Smith form gives full linear
image ranks2 at2,4 at3, and7 at every prime≥5; the displayed binary-image
cardinalities are separately finite-exact, not inferred from those ranks.
This is an arithmetic observation code with a genuine recovered-information
statement, not a binary[24,12,8] code or a statement about scalar residues.

## 4. The canonical paired code and its lattice

On the 24 presentation coordinates, one-choice-per-pair vectors form the
affine set

    A={ (ε_0,1−ε_0,...,ε_11,1−ε_11): ε∈F2^12 }.

Its difference code is the canonical paired repetition code

    C_pair={ (a_0,a_0,...,a_11,a_11): a∈F2^12 }.

It has dimension12, is self-dual for the binary dot product, and has minimum
distance2; its weight enumerator is (1+z²)^12. These facts hold for any twelve
pairs, so they are structural but not tournament-specific. Its usual scaled
Construction A lattice

    Λ(C_pair)={z/√2: z∈Z^24, z mod2∈C_pair}

is isometric to Z^24: each coordinate pair has orthonormal basis
(1,1)/√2, (1,−1)/√2. It is odd, minimum norm1, and has48 norm1 vectors.
It is not the even rootless Leech lattice, nor the extended binary Golay code.
The latter has parameters [24,12,8] and enumerator
1+759z^8+2576z^12+759z^16+z^24 (standard cited comparison, not derived here).

Adding parity checks to C_pair alone cannot repair this: a proper subcode has
dimension<12, so cannot remain self-dual of length24. A Golay construction
requires a new cross-pair coupling, not only deletion of the short words.

**Small-support obstruction.** A nonidentity coordinate permutation moving
fewer than8 coordinates cannot preserve a binary self-dual code of minimum
distance8. For every c, c+πc would be a codeword supported on fewer than8
positions, hence zero. Thus π fixes every word coordinatewise. If π moves i
to j, all codewords have c_i=c_j, putting the weight2 word e_i+e_j in the dual,
a contradiction. In particular a single port-pair swap cannot be a Golay
symmetry under any labeling. This does not prohibit simultaneous swaps of
four or more pairs, or prove that such swaps actually occur in our data.

The natural pairing symmetry is not automatically a symmetry of the
continuation observer. Theorem1 proves that for m≥3 each nontrivial binary
port-choice change changes the aggregate F. The full tensor of coefficients,
not the pair-forgetting map, is what retains that information.

## 5. An Eisenstein completion, and why it still is not Leech

A canonical order-three formal operator on each rank-two pair is

    J=[[0,−1],[1,−1]],  J²+J+I=0.

It preserves the A2 Gram form G=[[2,−1],[−1,2]]. On12 pairs this supplies an
Eisenstein-module completion with positive form A2^12, determinant3^12 and72
norm2 roots. It is not a Leech-lattice identification. Any lattice containing
this actual A2^12 as an isometric sublattice retains those roots; merely
adding glue cannot make it rootless. A different sublattice, scaling and/or
neighbor operation would have to be specified.

J is formal: it sends one positive basis vector to the negative sum of the
two. A genuine linear automorphism of the positive two-dimensional cone
permutes its extreme rays; an order-three such automorphism is the identity.
Thus a nontrivial J cannot describe a positive, pair-preserving port operation
without additional signed states. This assertion is per invariant pair, not
about permutations of three or more different pairs.

## 6. Historical connections and what is new

- [THM-481, Paley gauge/Gleason](../../../../01-canon/theorems/THM-481-paley-gauge-generates-gleason.md)
  and [THM-4857, GF2 tournament automaton](../../../../01-canon/theorems/THM-4857-the-gf2-automaton-of-a-tournament-classifies-the-merged-metagraph-exactly-through-order-five-circulants-are-additive-cellular-automata-and-paley-tournaments-give-quadratic-residue-codes.md)
  already give a genuine Paley23→Golay pathway. Our negative paired-code
  comparison does not deny that distinct, established map.
- [THM-480, skew-tower row code](../../../../01-canon/theorems/THM-480-skew-tower-row-code-is-dplus.md)
  explains paired doubling plus glue. Its [June reflection](../../../../07-reflections/border-remembers-doubling-forgets-the-gleason-ring-from-two-paley-tournaments-cbx0611s3.md)
  contrasts arithmetic retained by bordering with memory erased by doubling.
  The present rank7 aggregate keeps twelve binary choices despite rational
  dependence, a different exact example of representation-dependent memory.
- [THM-4752, score lattice](../../../../01-canon/theorems/THM-4752-tournament-scores-in-ii-25-1-backward-length-and-the-a24-deep-hole.md)
  supplies an exact 25-vertex score map into a Lorentzian lattice and its
  A24 cusp, while explicitly recording that score equality loses H. That
  invariant is not the present marked-port map.

The new exact data are the inherited24-profile lattice, all-length binary
separation, and the precise G(m) arithmetic schedule. Coding/lattice facts
are standard elementary comparisons. None establishes cluster mutation of
physical ports, an intrinsic Golay structure, or Chan–Pak completeness.

## Reproduction

From repository root:

    python3 -B 05-knowledge/results/session291_cluster_ports_20261113/ports/verify.py
    python3 -B -O 05-knowledge/results/session291_cluster_ports_20261113/ports/verify.py

The script independently checks the maximal-minor index, exact early scalar
separation, uniform dominance bound, literal pairing/selected masks, low-prime
valuation controls and256-length gcd probes. The latter probes check the
proved formula; they are not its proof. Source SHA and exact matrices are in
the deterministic receipt. SymPy is used for Smith form only; the lattice
index has a separate integer Bareiss/minor check.
