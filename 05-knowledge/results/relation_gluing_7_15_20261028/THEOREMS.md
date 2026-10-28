# Relation gluing, coefficient change, and arithmetic fusion

2026-10-28. **PROVED elementary statements**, with classical antecedents
credited below; explicit arithmetic realizations and finite enumerations are
identified separately. No claim of literature priority or universal solution
of the repository's open problems is made.

## T1 Exact gluing describes where dependencies are created

Let R be a commutative ring and let A:X->M and B:Y->M be R-linear maps.
Write K_A=ker A, K_B=ker B and K=ker[A B]. There is a canonical exact sequence

```text
0 -> K_A direct-sum K_B -> K -> im A intersect im B -> 0,
                                  (x,y) |-> A x.
```

**Proof.** If Ax+By=0, then Ax=-By is in both images. The kernel of this
map consists exactly of Ax=By=0. Every common image z has preimages x,y
with Ax=z and By=-z, proving surjectivity. All maps are linear. QED.

Over a field this gives

```text
dim K = dim K_A + dim K_B + dim(im A intersect im B).
```

Splittings into direct summands require choices; the sequence itself is
canonical. Over general rings a splitting need not exist. For example,
take A=B:Z->Z/2 to be reduction. Then K={(x,y):x+y even}, and the quotient
by 2Z direct-sum 2Z is Z/2. This sequence cannot split, since K is torsion-free.

For an inhomogeneous target q, the allowable intermediate values are
`im A intersect (q-im B)`. If nonempty, this is a torsor for
`im A intersect im B`; each intermediate value has a fibre which is
a torsor for K_A direct-sum K_B. Thus the same statement controls actual
solution gluing, not just counts of equations.

For a rooted binary tree whose leaves are columns of a matrix A, let E_v
be the leaves under v and K_v=ker A_(E_v). At an internal node with children
l,r, the quotient

```text
K_v / (K_l direct-sum K_r) ~= span(A_(E_l)) intersect span(A_(E_r)).
```

It records dependencies first joining those two branches. Choosing sections
produces a kernel basis assigned to lowest common ancestors. The individual
basis is not canonical, but each displayed quotient and its dimension are.
Every circuit meeting both children has a nonzero class there: otherwise
its two restrictions would already be dependencies on proper subsets.

Summing the dimensions over all internal nodes telescopes to
`|E|-rank A` when each leaf column is nonzero. Add the individual leaf
kernel dimensions if zero columns are present. This is the exact content
of the LCA decomposition; a tree shape alone does not determine the maps.

This is elementary module theory, organized for the present applications.
The related code-realization framework is classical:
[Forney, sections 5.5 and 6.3](https://arxiv.org/pdf/1306.6264).
The repo's [Session56 join theorem](../session56_rmq_lca_connections_20261025/correlated-gluing/THEOREM.md)
already treats actual compatible assignments. T1 specializes to linear
relations and exposes their kernels, image intersections, and extensions.
The incoming [Session85 linear-gluing work](../session85_arithmetic_codes_free_bits_20261028/LINEAR_GLUING_INFORMATION.md)
was read at commit dd908f4f88e7583a4391b6f1ad208f777848784c. Its row-space
LCA invoice counts independent equations; the kernel quotient here counts
dependencies among columns. For an n-column matrix their total dimensions
are rank(A) and n-rank(A). Neither is a count of all named constraints or
all circuits. Its exact projection-loss theorem is credited rather than
rederived as a new observer theory.

## T2 Projective binary dependencies have exact tree width

For columns indexed by a set E, write

```text
lambda(S)=rank A_S + rank A_(E minus S) - rank A.
```

This equals the dimension of the intersection of the two image spaces.
On a tree representing the zero-syndrome code, a cut can carry exactly this
intersection: the two partial syndromes must cancel. Leaf values and the
local rule that a parent syndrome is the sum of its child syndromes recover
every codeword, with a unique internal syndrome assignment.

There are `2^lambda(S)` possible live separator values. Fewer states cannot
realize this code on the same separating topology: two distinct syndromes
have valid completions, and merging them would admit a mismatched completion.
This is a statement about a separated constraint realization, not a lower
bound on every possible algorithm or numeric encoding.

**Theorem.** Take all nonzero columns of F2^r, r>=2, on
`n=2^r-1` labelled leaves. Over all binary decomposition trees,

```text
min_tree max_edge lambda(S_edge) = r-1.
```

The width convention here is lambda itself; conventions adding one give r.

**Upper bound.** Split the points into a linear hyperplane H minus zero
and its affine complement. The first has rank r-1; the second has rank r,
so their cut has lambda r-1. Recursively split the affine complement into
affine hyperplanes and continue by coordinates. Every proper descendant
piece of that affine complement has span dimension at most r-1. Every
piece inside H also has rank at most r-1. For any such piece S,
lambda(S)<=rank A_S<=r-1.

**Lower bound.** Suppress the degree-two root, giving a subcubic leaf tree.
A leaf-weight centroid has at most three branches, each with at most n/2
leaves; one has at least n/3. Thus some edge has both sides between n/3
and 2n/3. Each side has more than `2^(r-2)-1` points, so has rank at
least r-1. If both ranks were r-1, the whole projective point set would
lie in the union of two proper hyperplanes. Such a union has at most
`3*2^(r-2)-1 < 2^r-1` nonzero points. Therefore at least one side has
rank r, and lambda>=r-1. QED.

Thus the seven-point code has optimal two-bit separators, and the
fifteen-point code has optimal three-bit separators. The all-r statement
is about abstract binary configurations. The explicit arithmetic realization
below supplies r=3,4; it does not prove a prime realization for every r.

### T2a Exactly one fifteenth of Fano leaf trees are optimal

There are 11!!=10395 rooted, unordered binary trees on seven labelled leaves.
Exactly 693 have maximum lambda 2; the remaining 9702 have maximum lambda 3.

**Proof.** In the unrooted trivalent tree, only 3/4 splits can have lambda3;
they have lambda2 exactly when their three-point side is a Fano line.
There must be at least one 3/4 split. Two different such splits have disjoint
three-point sides, by compatibility of tree splits. Two Fano lines cannot
be disjoint, so an optimal tree has exactly one 3/4 split.

Choose its three-point Fano line in 7 ways. The three-point rooted side has
3 shapes (which pair is the cherry). The four-point side must split2+2,
since a 1+3 split would create another3/4 edge; its two cherries can be
paired in 3 ways. This gives7*3*3=63 unrooted trees. Each has 11 edges on
which to insert the root, giving 693. Since 10395=15*693, the proportion is
exactly 1/15 under the uniform distribution on these labelled rooted trees.
The independent recursion enumerates the entire10395-tree universe. QED.

This occurrence of 15 is a derived tree-count ratio, not a presumed
identification with the fifteen projective points.

## T3 A prime realization joins the seven and fifteen geometries

For distinct primes p_i, put `a_i=(v_l(p_i+1) mod 2)_l`. Unique
factorization and sigma's coprime product rule give

```text
sigma(product_(i in S) p_i) is square iff sum_(i in S) a_i=0.
```

The following explicit bank realizes every nonzero vector on coordinates
(2,3,5,7). The table gives d=sf(p+1) and p+1=d*a^2.

| binary vector as integer | d | a | p |
|---:|---:|---:|---:|
|1|2|2|7|
|2|3|1|2|
|3|6|1|5|
|4|5|2|19|
|5|10|3|89|
|6|15|2|59|
|7|30|1|29|
|8|7|6|251|
|9|14|1|13|
|10|21|2|83|
|11|42|1|41|
|12|35|2|139|
|13|70|9|5669|
|14|105|2|419|
|15|210|2|839|

Primality is checked by trial division. Each p is the least realizing prime
for its specified d: all smaller positive a were tested, and d*a^2-1 is
strictly increasing. These are finite witnesses, not prime-infinitude claims.
Notice that coordinate prime 7 is represented by input prime 251. Using atom4
instead would introduce the earlier2/4 support collision.

The first seven rows are exactly the prior arithmetic Fano bank. All fifteen
give the classical Hamming[15,11,3] code, hence 2048 square-total subsets.
There are35 circuits of size 3,105 of size 4,168 of size 5:308 total.
A circuit in rank 4 has at most5 elements. Every zero-sum support of size
3,4,5 is minimal here: a contained circuit of size 3 or 4 would leave a
zero-sum subset of size 1 or 2, impossible for distinct nonzero columns.
The enumerated weight counts therefore establish the displayed circuit counts.

More generally, character orthogonality gives the Hamming weight polynomial

```text
W_r(z)=2^(-r) [(1+z)^(2^r-1)
       +(2^r-1)(1+z)^(2^(r-1)-1)(1-z)^(2^(r-1))].
```

**Proof.** Replace the zero-syndrome indicator by the average of all binary
characters. The zero character contributes the first product. Each nonzero
character is positive on 2^(r-1)-1 nonzero columns and negative on 2^(r-1).
Factor the independent subset sum over columns. QED.

This coding identity is classical and already appears in related form in
[Session33, section 5](../session33_continuation_states_20261023/FANO_CONTINUATION_SPECTRUM.md).
The contribution here is the explicit fifteen-prime realization and its
integration with the following geometric and integral structures.

**Concurrent credit:** incoming [Session85 arithmetic codes](../session85_arithmetic_codes_free_bits_20261028/ARITHMETIC_FREE_BITS.md)
and [the free-bits packet](../free_bits_signed_incidence_20261028/README.md)
publish the identical least-prime bank and circuit counts. They were read
before close-out. This packet independently reproduced that shared result;
its additional work is the tree-width/count, integer-module distinction
and symmetry calculation, and the relation-to-object construction below.
Session85 also supplies an extended sixteen-prime code and an exact
five-of-seven parity-access guarantee; those are credited, not rederived here.

### T3a A single model explains several distinct roles of five and fifteen

Let

```text
V={even subsets of {0,1,2,3,4,5}} / <the full six-set>.
```

It is a four-dimensional binary space. Every nonzero class has a unique
two-element representative, a duad. Use basis [01],[02],[03],[04] to identify
it with the table. The alternating pairing is
`B([S],[T])=|S intersect T| mod 2`. It is well-defined, and nondegenerate
because the radical before quotienting consists of the empty/full sets.

Every projective line is a triple{x,y,x+y}. For two distinct duads:

* disjoint duads sum to the remaining pair, after complementation; the
  three form a perfect matching of six labels;
* intersecting duads sum to their two unshared endpoints; the three form
  a triangle on three labels.

Consequently the 35 arithmetic three-prime circuits split into 15 matching
triples and 20 triangle triples. The15 matching triples are the isotropic
lines of W(3,2), the classical doily. This is the repo's existing
[duad/Lagrangian carrier](../cardinality_bridges_20261018m.md#3-an-explicit-fifteen-bridge-including-its-action-boundary)
with an explicit arithmetic labelling. Its symplectic form is extra data;
the square-class matrix alone does not select it.

A spread partitions the 15 points into 5 disjoint three-point lines. The
checker finds all 56 projective spreads, of which6 are entirely isotropic.
One isotropic spread gives these five square-total arithmetic triples:

```text
{7,59,29}, {2,13,41}, {5,139,839},
{19,83,419}, {89,251,5669}.
```

The six isotropic spreads are the six pentads (one-factorizations of K6)
recovered in the [fifteen-kernel packet](../session52_saturation_residues_metagraph_20261024/fifteen-kernels/README.md).
Their5 is the number of rounds in such a factorization. An abstract
projective spread can also be obtained from the five one-dimensional F4
subspaces of F4^2: each has three nonzero vectors. This requires an F4
structure; it is not an identification of arithmetic prime 5 with a vector.

In the first seven rows, the duads are the six edges of K4 on{0,1,2,3}
plus the remaining pair45, represented by prime 29. Its Fano lines split
into four K4 triangles and three opposite-edge-pair lines through 45.
This recovers the earlier [tetrahedral Fano duality](../tetrahedral_fano_duality_20261021/session23b-duality/PROOF.md).
The distinguished point in this marking is29, not the input prime 7.

### T3b Literal cluster covers differ from algebraic decomposition width

By [Session56's hierarchy-cover theorem](../session56_rmq_lca_connections_20261025/hierarchy-covers/REPORT.md),
the number of same-label trees needed to display every requested subset as
a literal node cluster equals the chromatic number of their crossing graph.

The seven Fano lines pairwise intersect and none contains another: they need
seven trees as literal clusters, despite T2's single algebraic tree with
two-bit separators.

For the fifteen doily lines, at most five can share a tree because they must
be disjoint. Three trees would therefore require three disjoint spreads.
The exhaustive six-pentad certificate shows that every two isotropic spreads
share exactly one line, ruling this out. A four-colouring is provided in
normal.json. Thus the exact cluster-cover number is4.

For all 35 projective lines, at least7 trees are needed by the same
five-disjoint-lines bound. The receipt gives seven disjoint spreads covering
all 35 lines, proving that7 suffices. The spread and colouring certificates
are **FINITE-EXACT**, with complete finite universes and checkable witnesses.
These three tree counts and the separator dimensions answer different
questions; none can be substituted for another.

### T3c Relations at seven become objects at fifteen

Let P=(7,2,5,19,89,59,29), in binary-column order, and let C=ker H7.
For c in C, put N(c)=product of P_i over the support of c. Then

```text
N(c+d) = N(c) * N(d) / gcd(N(c),N(d))^2.
```

This is an exact group law on the sixteen square-total squarefree products
of the seven primes. Its identity is1. Choose the following basis of C,
writing supports by column numbers rather than arithmetic primes:

```text
g1={1,2,3}, g2={1,4,5}, g3={2,4,6}, g4={1,2,4,7}.
```

Their private coordinates 3,5,6,7 establish independence; all have zero
syndrome, so they span the four-dimensional kernel. The coordinate map
beta:C->F2^4 sends these generators to the four standard basis vectors.
For c!=0 label its point by the row-prime p_beta(c) from the fifteen table.
For distinct nonzero c,d,e we then have the exact equivalence

```text
c+d+e=0  iff  sigma(p_beta(c) p_beta(d) p_beta(e)) is square.
```

Thus the fifteen nonidentity square-total products of the first bank
become fifteen prime-labelled points of the next geometry. Its35 lines
are the three-term laws among those old dependencies. This is a literal
relation-to-object construction, not just the equality of two counts.
The construction needs the chosen code basis and prime realization; it
does not preserve numerical size, sigma value, or the original prime
labels. The basis is not canonical. All256 group-law pairs and all 455
three-point choices are independently checked against integer arithmetic.

This also separates two previously conflated fifteens: the 15 nonzero
words of Hamming[7,4,3] are inputs to this map; the 15 coordinates of
Hamming[15,11,3] are its outputs. The latter code has 2048 words.

## T4 Changing coefficients can create dependencies from torsion

Let A:Z^n->Z^r be an integer matrix and let ell be prime. Set C=coker A.
There is a canonical exact sequence

```text
0 -> ker(A)/ell ker(A) -> ker(A mod ell) -> C[ell] -> 0,
                                  x mod ell |-> [A x / ell].
```

Here C[ell]={c:ell*c=0}; it is Tor_1^Z(C,F_ell).

**Proof.** If A*x is divisible by ell, the displayed quotient is integral.
Changing the lift x by ell*y changes A*x/ell by A*y, so its class in C is
well-defined. Its class has order dividing ell. Conversely, if ell*z=A*x,
then x modulo ell lies in the modular kernel and maps to[z].
If A*x/ell=A*y, then x-ell*y is an integer kernel vector, proving exactness
in the middle. Finally the integer kernel is saturated in Z^n because
Z^r is torsion-free, proving the injection on the left. QED.

The full-rank special case was already used in
[THM-2703 / boundary-tree discriminant](../../../01-canon/theorems/THM-2703-c3-boundary-tree-arm-determinant-standard-plane-gate.md),
section 2. This is the elementary universal-coefficient mechanism, not a
new homological-algebra theorem.

For the arithmetic matrix retaining full valuations of p+1 on(2,3,5):

| source-prime subset | nonzero integer Smith factors | rational nullity | mod 2 nullity |
|---|---|---:|---:|
|{2,5,7}|1,1|1|1|
|{7,19,89}|1,1,6|0|1|

The first relation lifts from integers: coefficients(3,-3,1) on the
columns for(2,5,7) give zero. The second parity relation comes entirely
from the 2-torsion of Z/6. Its valuation matrix is

```text
       7 19 89
v_2    3  2  1
v_3    0  0  2
v_5    0  1  1
```

with determinant-6. Both subsets are Fano lines, but their integral
mechanisms differ. Gluing their integer image spaces and then reducing
mod 2 can therefore differ from reducing first and gluing: at the split
{7}|{19,89}, the integer intersection is zero while the modular
intersection has dimension 1.

For the signless triangle incidence matrix, the Smith factors are1,1,2.
Its mod 2 cycle likewise comes entirely from torsion. This is an exact
bridge between the square-class arithmetic and the odd-cycle integral
obstruction behind the earlier graceful/square-sum comparison.

Retaining only the full matrix's cokernel is insufficient: the complete
seven-column arithmetic matrix has trivial cokernel, while its displayed
three-column subset has Z/6. The appropriate classical enlargement assigns
every column subset S the module `Z^r/span_Z(A_S)`:
[Fink–Moci, definition 2.1 and equation 2.1](https://arxiv.org/pdf/1209.6571).
For reconstruction tasks retain the presented maps as well as abstract
module isomorphism types. The elementary calculations here do not import
every structural theorem about matroids over rings.

### T4a The arithmetic lift breaks Fano symmetry from 168 to 2

The binary Fano configuration has GL(3,2) as its full automorphism group:
every incidence permutation preserves x+y and therefore extends linearly.
Its order is7*6*4=168.

Decorate every subset by its full integer quotient-module type. An exhaustive
check of those168 permutations on all 128 subsets leaves exactly the identity
and the involution

```text
19 <-> 89,    59 <-> 29;    7,2,5 fixed.
```

This is **FINITE-EXACT**, with an explicit positive certificate. In the row
coordinates(2,3,5), the integer matrix

```text
G = [1 0 -3; 0 1 -2; 0 0 -1]
```

has determinant-1 and G^2=I. It fixes the first three columns and carries
the other four to negatives of their indicated partners. Hence it preserves
every integer-span quotient type. The exhaustive check proves that no other
binary automorphism does so. Preserving these modules is not preserving the
positive arithmetic labels: the signs and divisor-sum magnitudes change.

This parallels the mechanism of the inherited81/27 tournament example:
enriching a combinatorial object by its actual decorations restricts the
symmetry to their stabilizer. It does not identify the two group actions.

## T5 Every integer has a prime-power presentation, and overlap can be fused

Fix finitely many primes P and exponent bounds E_p. Choose exactly one
e_p in{0,...,E_p} for each p. The map

```text
(e_p)_p |-> n=product_p p^e_p
```

is a bijection onto that exponent box. Put S_p(e)=1+p+...+p^e.
The same assignment gives exact functions

```text
sigma(n)=product S_p(e_p), tau(n)=product(e_p+1),
Omega(n)=sum e_p, omega(n)=sum [e_p>0].
```

For any k>=2, form columns `v_l(S_p(e)) mod k` and one-hot variables
x_(p,e), with sum_e x_(p,e)=1. Then sigma(n) is a kth power iff the
selected column sum is zero in the module (Z/k)^L. For prime k this is
linear algebra over a field; for composite k it is a module problem, not
automatically a matroid. Every positive integer occurs in some such box.

The same presentation evaluates F,U,S exactly, for n>=2:

```text
F=tau-2, U=omega-[n prime], S=2^omega-1-[n squarefree].
```

Abundance is sigma(n)>2n. Liouville is(-1)^Omega; Möbius has that sign when
all exponents are<=1 and is zero otherwise. The inherited FUS shape theorem
and its abundance classification remain additional results about these
evaluations, not consequences of a rank count.

This presents actual numbers and all these predicates in one system.
It retains local choices, so it is not promised to be a compressed model.

### Exact fusion defect for unrestricted multiplication

Let v denote the full prime-valuation vector, with integer entries. Define

```text
d_p(a,b)=v(S_p(a+b))-v(S_p(a))-v(S_p(b)).
```

For all positive integers n,m,

```text
v(sigma(nm))=v(sigma(n))+v(sigma(m))
             +sum_p d_p(v_p(n),v_p(m)).
```

Only shared primes contribute. Moreover

```text
d_p(a,b)+d_p(a+b,c)=d_p(b,c)+d_p(a,b+c).
```

Both identities follow by cancellation of the displayed valuations.
The second is the cocycle identity ensuring that fusion is independent
of parenthesization. This is explicitly a coboundary of v(S_p(e)), so no
nontrivial cohomology class is claimed. Reduction modulo k preserves the
identities.

The earlier forbidden2/4 multiplication can therefore be performed with
its exact correction:

```text
sigma(8)/(sigma(2)*sigma(4))=15/21=5/7.
```

Modulo squares the defect is35. Thus the formal square1764 for atoms
{2,4,5,13} becomes actual sigma(520)=1260=36*35. The operation is no
longer discarded; its repeated-prime fusion supplies the missing factor.
Liouville has zero multiplicative defect everywhere. Möbius's zero states
require its separate repeated-prime cancellation, as the prior
[incidence-algebra account](../liouville_mobius_eleven_20261027/THEORY.md)
explains.

## T6 Exact weighted assembly

For finite local domains with weights in a commutative semiring, define
M_v(z) as the total weight of assignments below a tree node having boundary
value z. On disjoint leaf scopes, linear boundaries obey

```text
M_v(z)=sum_(u+w=z) M_l(u) M_r(w).
```

This follows by the bijection between a pair of child assignments and
their union; weights multiply once per original local factor. Boolean
weights test existence, natural-number weights count, and polynomial
weights retain the selected variables. Where scopes share variables, use
the actual equality-of-assignment fibre product from Session56 instead.
This is classical exact sum-product assembly, not an approximation by
independent marginals.

The extension beyond observation is concrete: specify generators and
relations, construct their solution spaces, glue them, retain coefficient
change and torsion, and only then evaluate the desired quantity. Nonlinear
restrictions such as positivity, all-different labels or primality remain
explicit factors. Their presence can make the resulting constraint
decomposition large even when the underlying unlabelled graph is a tree.

## T7 Beyond integer valuations: units, ideal classes, and Selmer groups

This is a classical number-field enlargement, not a new Selmer theorem.
Let K be a number field with ring of integers O, unit group U, fractional
ideal divisor group D, and ideal-class group Cl. For m>=2 put

```text
V_m = K^*/K^{*m},
E_m = {[alpha] in V_m : every finite-prime valuation of alpha is 0 mod m}.
```

There are exact sequences

```text
0 -> U/U^m -> E_m -> Cl[m] -> 0,
0 -> E_m -> V_m -> D/mD -> Cl/mCl -> 0.
```

The first is the empty-S case of the standard number-field Selmer
sequence; see the [Sage reference and implementation](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/selmer_group.html).
For m=2, [Dummit–Voight, section 3](https://arxiv.org/abs/1702.00092)
describes the same ideal-to-class map. Here the sequences are written
additively as abelian groups, even though field elements multiply.

**Proof.** If all valuations of alpha are multiples of m, unique ideal
factorization gives (alpha)=a^m. Send its power class to [a]. Multiplying
alpha by beta^m changes a by the principal ideal(beta), so this is
well-defined. Its image is exactly Cl[m]: an ideal class killed by m has
a representative a with a^m principal. Its kernel consists of principal
a=(beta), whence alpha/beta^m is a unit. A unit that is an mth power in
K is an mth power of a unit, since its root has zero valuation everywhere.
This proves the first sequence.

For the second, send a field element to its principal ideal divisor
modulo m. Its kernel is E_m by definition. A divisor D represents a
class in mCl exactly when D-mA is principal for some fractional ideal
divisor A; this is exactly membership in the image of V_m. Surjectivity
onto Cl/mCl follows from the definition of Cl. QED.

These formulas explain two distinct failures of naive parity coordinates:
some formal ideal-parity vectors cannot be principal; and zero ideal
parity need not mean a field element is a square. Units and class-group
torsion account for the latter. They do not assert that every associated
field extension is unramified at primes dividing m. That stronger local
condition needs further data.

For positive rationals, U contributes no positive nonsquare class and
Cl(Z)=0, recovering the original exact prime-parity decoder. For
Eisenstein integers the ideal class group is trivial, but the unit -1
is still not a square in Q(sqrt(-3)); its finite valuations are all zero.
Indeed (a+b sqrt(-3))^2=-1 with rational a,b forces ab=0 and then either
a^2=-1 or b^2=1/3, both impossible. Thus even this familiar grid needs a
unit coordinate if arbitrary nonzero field elements are allowed.

For the earlier two-grid fields K=Q(sqrt(-d)), d in{5,13,37}, put
a=(2,1+sqrt(-d)). The inherited [THM-4825 grid account](../two_grids_ring_or_ideal_class_13_37_20261027.md)
proves a^2=(2), with a nonprincipal. There is also a short direct check:
the quotient O/a has two elements; products of its generators are in(2),
and their generated ideal contains2, so equality holds. Explicitly
2(1+t)-(1+t)^2=1+d is2 modulo4, and 4 is another product generator.
Nonprincipality follows because x^2+d*y^2=2 has no integer solution.
Thus the square class of 2 has zero ideal parity but maps to the nonzero
class[a] in Cl[2]. Also2 is not a square in K: squaring a+b sqrt(-d)
again forces ab=0, giving a^2=2 or -d*b^2=2.

This identifies the user's centred-grid class with a specific missing
coordinate in the extended square-class theory. The obstruction already
exists for all squarefree d>1 with d=1 mod 4; the special role of 5,13,37
is the inherited class-number-two classification, not a new consequence
of the displayed parity calculation.

For finite arithmetic banks over K, form columns in the full subgroup of
V_2 generated by the chosen values. Its kernel again tests squarehood.
Computing it through ideals alone requires the unit and class data above.
The repo's [Dedekind conductor theorem](../fixed_twists_conductors_reverse_potential_20261025/conductors/DEDEKIND_CONDUCTORS.md)
already supplies a related warning: individual character components can
be nonfree ideals even when a total rank looks innocuous. Its assumption
that2 is invertible is retained; this paragraph does not remove it.
