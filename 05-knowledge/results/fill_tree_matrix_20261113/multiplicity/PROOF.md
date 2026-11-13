# Exact spatial fill inside a Hamilton-count fibre

**PROVED IN NOTES / FINITE-EXACT controls.** The neutral growth geometry,
its weighted uniform-fill process, and its primitive-cell difference law
are proved below. The underlying SCC multiplicity formula is inherited,
not new. No global odd-spectrum coverage is assumed; Chan-Pak remains
**OPEN**.

## 1. Inheritance and the four numbers

The closest proved mechanism is the actual-count SCC-word theorem in
[Session24's multiplicity proof](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md),
independently reviewed in
[its review](../../actual_hamilton_multiplicities_20261021/session24-multiplicities-review/REVIEW.md).
The corrected near miss is
[THM-1900, insertion-response calculus](../../../../01-canon/theorems/THM-1900-the-insertion-response-calculus.md):
deleting a newly inserted vertex from a Hamilton path need not leave a
Hamilton path. Its repaired all-order proof retains the additional
orphan paths and proves neutrality exactly at cuts between SCCs.
The hostile here is H9: forgetting the number of nontrivial components
makes a finite-difference cell counter negative. The least-used sidecar
is the number of literal insertion cuts mapping to the same quotient
successor.

Throughout H(T) counts actual directed Hamiltonian paths. N_n(h) counts
ordinary isomorphism classes, with converse-distinct classes kept distinct.
It does not count H/|Aut|, path-normalized presentations or merged classes.

The complete strong inventories for our four numbers are:

| h | Strong atoms, all orders | First order | N_n(h), all n>=1 |
|---|---|---|---|
| 5 | one rigid, self-converse atom of order4 | 4 | max(n-3,0) |
| 17 | two rigid, self-converse atoms of order6 | 6 | 2 max(n-5,0) |
| 7 | none | none | 0 |
| 21 | none | none | 0 |

These facts are inherited and independently recounted here. The script
enumerates every labelled tournament through order6, canonically tests
the indicated strong classes, and separately counts paths of the three
positive atoms by permutations. Busch's bound

    H(S)^3 >= 5^(|S|-1)

for a strong tournament S makes the census complete for these four
targets: 5^6 > 21^3 excludes every order>=7. The primary source is
[Busch, A Note on the Number of Hamiltonian Paths in Strong Tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
The exact lower-bound statement is used, not the loose wording that an
irrational expression is an attained minimum at every order.

For prime5 or17, any realization has exactly one nontrivial SCC, whose
count is that prime. Its other SCCs are singleton vertices. For21, any
nontrivial SCC count is a divisor in {3,7,21}; neither7 nor21 exists, and
a product of3s cannot be21. This proves the four all-order multiplicity
formulas from the strong inventories.

Thus5 and17 represent one and two populated geometric components,
respectively. Seven and21 represent empty fibres. A shared arithmetic
pattern does not make those roles interchangeable. The previous
unit-ready minimum theorem's appearance of5 and17 is a more restricted
rooted problem; its equality types are not this unrooted inventory.

## 2. The exact neutral-growth orthants

Fix an ordered word W=(C_1,...,C_r) of nontrivial strong tournament
isomorphism types. Put M=sum |C_i| and h=product H(C_i). For
g=(g_0,...,g_r) in N^(r+1), form the ordinal sum

    T(W,g)=TT_(g_0) -> C_1 -> TT_(g_1) -> ... -> C_r -> TT_(g_r),

where every arc between blocks points from left to right, and TT_0 is
empty. Singleton SCCs in the same gap form the displayed transitive run.

**Orthant theorem.** Ordinary isomorphism classes with fixed W are in
bijection with N^(r+1). Their order is M+sum g_i and their Hamilton count
is h. An H-preserving one-vertex birth changes exactly one coordinate
g_i to g_i+1. Conversely every such coordinate change is realizable.
There are exactly g_i+1 distinct literal attachment patterns producing
that same quotient successor.

Proof. The SCC condensation of a tournament is a linear order. Every
Hamilton path traverses the components in that order, so H is the product
of component counts. Its ordered SCC decomposition is canonical; this
recovers W and every g_i, proving injectivity. THM-1900's repaired theorem
says that H-preserving births are exactly terminal SCC cuts. In gap i
there are g_i+1 such cuts, including both ends. They insert a new singleton
between two consecutive positions in that run and all produce the same
unlabelled run of length g_i+1. No other H-neutral birth is possible. QED.

The underlying quotient graph has a coordinate edge g->g+e_i. For r>=1
it has commuting diamonds, so it is not itself a tree. A spanning tree
is obtained by giving each nonzero g the parent g-e_j, where j is the
largest index with g_j>0. The sum of coordinates strictly decreases,
so every cell has a unique finite route to the core at g=0.

This is the full H-neutral birth/deletion graph on ordinary types. It
does not classify arc surgery, marked-root repairs, or general metagraph
edges. Different W lie in different components even when their h agrees.
The H17 fibre therefore has two orthants, not one connected neutral tree.

## 3. A literal perfect-fill process, at every size

At each stage choose uniformly among every H-preserving attachment to
the current named tournament. At gap state g with d=sum g_i, the exact
quotient transition probabilities are

    Pr(g -> g+e_i)=(g_i+1)/(d+r+1).                            (1)

**Uniform-fill theorem.** Starting from W with no singleton padding,
after d births every weak-composition cell of size d has probability

    1/binom(d+r,r).                                          (2)

There are d! literal histories, with vertex names fixed by birth time,
over each such cell. Thus both the support and its probability measure
fill the entire integer simplex at every stage.

Proof by weighted histories. The denominator of every length-d history
is (r+1)(r+2)...(r+d). For a fixed gap vector g, the number of coordinate
orders is d!/product g_i!. Along each coordinate order the successive
numbers of literal cuts multiply to product g_i!. Thus there are exactly
d! literal histories above g, independent of g. Multiplying by the common
history probability yields d! r!/(d+r)!=1/binom(d+r,r). Distinct histories
produce distinct final named graphs: restriction to the first M+j names
recovers every preceding birth. QED.

Equivalently, (1) is the classical Pólya urn with r+1 colours, each
starting at weight1. This is a derived exact physical implementation,
not a claim of novelty for the urn identity. Upon scaling by d, the
uniform weak-composition law has the uniform-simplex limit; the finite
uniformity above is stronger than merely an asymptotic distribution.

**Hostile.** Uniform choice among the r+1 *quotient successor types* is
a different process. With r=1 and d=2 it gives gap probabilities
(1/4,1/2,1/4), whereas literal uniform attachment gives(1/3,1/3,1/3).
Discarding the g_i+1 multiplicity destroys perfect uniform fill.

This perfect fill never manufactures a new strong core W. It fills each
already-existing fibre component. It explains how occurrences of5 and17
proliferate, but it cannot fill an empty target fibre. Neutral growth
from a finite bank preserves that bank's finite set of Hamilton counts.

## 4. Counting and recovering newly born core cells

Let B_(h,r,m) count ordered nontrivial strong words W of product h,
length r and total order m. Let N_(n,r)(h) count ordinary tournaments of
order n, Hamilton count h and exactly r nontrivial SCCs. The inherited
auxiliary order-zero convention is H(empty)=1 and N_(0,0)(1)=1, with
all other order-zero slices zero. The inherited
SCC-word formula is the integer-simplex count

    N_(n,r)(h)=sum_m B_(h,r,m) binom(n-m+r,r),                 (3)

with each term zero for n<m. No quotient by permutations of components
is allowed: their order is intrinsic. For fixed h there are only finitely
many W, since every nontrivial component count is at least3, and Busch's
bound bounds the order of every component whose count divides h.

Define the backward difference nabla f(n)=f(n)-f(n-1), with f(n)=0 for
n<0. Then the exact primitive-cell recovery law is

    nabla^(r+1) N_(n,r)(h)=B_(h,r,n),                         (4)

for every integer n>=0, not only eventually. Indeed the ordinary series
of the r-th slice is B_(h,r)(z)/(1-z)^(r+1), and taking a backward
difference multiplies by1-z. This is also the iterated boundary count
for the lattice simplex in(3).

For total N_n(h), the eventual degree is the maximum feasible r, but
unmarked higher differences need not count positive primitive cells.
The explicit hostile is

    sum N_n(9) z^n = 2z^5/(1-z)^2 + z^6/(1-z)^3.

The two order-five H9 atoms have r=1; C3->C3 has r=2. Multiplying the
total series by(1-z)^3 gives 2z^5-z^6. The negative coefficient at n6
does not mean a negative number of births. It records lost dimension
information. Keeping r restores nonnegative primitive coefficients in(4).

The existing all-order formulas now have a direct geometric reading:
H5 is one length-one core at order4; H17 is two length-one cores at
order6; H7 and H21 have no core cells at any rank. The earlier Session24
degree-defect and repunit laws remain inherited, including
[THM-4749, the repunit degree law](../../../../01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md).

### A filled target can still have no incoming extension

The known all-order obstruction that no H27 tournament contains an H13
tournament as an induced subtournament was already proved in
[the November2 induced-obstruction packet](../../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md).
The coordinator's current all-target test independently recovers it.
Here is a smaller independent explanation using the inherited complete
atom inventory, with48 new five-subset checks and the neutral-cone
interpretation. This is not a new exclusion theorem.

H13 has one strong atom, of order5. Every H13 tournament contains it.
An H27 tournament has one of exactly ten nontrivial SCC words:

* five single strong H27 atoms, two of order6 and three of order7;
* four ordered words using C3 and one of the two strong H9 atoms;
* the single word C3->C3->C3.

Busch's inequality5^7>27^3 excludes strong H27 atoms of order>=8, so
the inventory is globally complete. A strongly connected induced
subtournament lies wholly in one SCC of its host. A component with
H3 or H9 cannot contain the H13 core, since H is monotone under induced
extension by the repaired insertion theorem. The remaining five strong
H27 atoms have three converse representatives, with rows

    (30,52,56,50,32,1),
    (62,124,120,112,32,64,17),
    (118,124,40,113,100,64,4).

All6+21+21 five-vertex subsets of these representatives have H different
from13, independently counted in the receipt. Converse preserves those
counts. Hence no H27 host, at any order or with any singleton padding,
contains any H13 source. Pure births cannot realize the arithmetic edge
13->27, even after arbitrary same-count source preparation.

Nevertheless the exact ordinary multiplicities are

    N_n(13)=n-4                  (n>=5),
    N_n(27)=n(n-5)(n-4)/6        (n>=6).

There are eventually quadratically more H27 objects than H13 objects,
yet none is an induced extension of any H13 object. This directly
separates populated fibres, spatial fill within fibres, and physical
edges between fibres. It does not exclude repairs of old arcs, deletion
of old vertices, or a different arithmetic parent of27.

## 5. What a global perfect fill would have to prove

Let S_n be the set of Hamilton counts attained by n-vertex tournaments.
Source padding gives S_n subset S_(n+1). Let b_n be the least positive
odd number outside S_n other than7 and21. It exists because S_n is finite.
Then the following are equivalent:

1. Every positive odd h other than7 and21 is attained at some order.
2. For each finite bound L, some S_n contains every allowed odd h<=L.
3. b_n tends to infinity.

Proof. For1=>2 take the maximum of the finitely many witness orders and
pad each to that common order. Statements2 and3 are the same first-gap
condition. Statement2 immediately supplies a witness for any fixed h.
The separate exclusion theorem supplies the assertion that the final
complement is *exactly*{7,21}. QED.

The complete inherited census gives the following first gaps:

| n | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| b_n | 5 | 9 | 17 | 35 | 63 | 611 |

The first four columns also agree with this session's independent
labelled census. The last two are inherited exact data, not rerun here.
They do not prove b_n->infinity. Repeated occurrences of earlier values,
or density zero of an omitted set, do not imply the first-gap condition.
For example a nested toy support can fill every odd number except one
fixed additional value while having missing density zero.

In geometric language the real missing theorem must supply at least one
core cell over every allowed target label. The neutral orthants and
their perfect fill then propagate those cells to all larger sizes.
Construction, component-level enumeration, and spatial propagation are
three distinct obligations; equations(1)-(4) settle only the latter two
once the actual strong core words have been supplied.

There is a stronger fixed-bank obstruction. If only finitely many strong
cores are allowed, neutral growth and ordinal sums produce counts in the
multiplicative monoid generated by their finitely many H values. If P is
the set of primes dividing those generators, every produced count is
P-smooth. Up to X there are at most

    product_(p in P) (1+floor(log_p X)) = O((log X)^|P|)

such counts. Infinitely many odd primes lie outside P and are missed.
This does not apply to a graph substitution that creates genuinely new
strong cores; it pinpoints why neutral fill plus a finite SCC alphabet
cannot supply the desired global coverage, however many histories occur.

## 6. Controls and reproducibility

Run `python -B 05-knowledge/results/fill_tree_matrix_20261113/multiplicity/check.py`
and repeat with `-O`. The receipts must agree as JSON.

The declared universe is all33,867 labelled tournaments of orders1..6;
all1,098 one-vertex attachments to labelled parents of orders1..4;
the literal neutral histories over one C3 through five births and over
C3->C3 through three births; and5,004 rational urn cells over r0..5,
d0..8. Literal atom paths are separately enumerated by permutations.
Integer Busch inequalities validate the all-order four-target cutoff.
Rankwise differences and the H9 negative unmarked coefficient are exact
controls. Old census and proof hashes are pinned in the receipt.
