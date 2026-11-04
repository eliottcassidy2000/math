# Hamilton frame counts and the exact first symmetric seven prime fibres

**PROVED reductions and marked-module inequalities / FINITE-EXACT complete
bases and inventories. The assertion that7 and21 are the only missing
positive odd Hamiltonian-path counts remains OPEN.**

The principal new classification here is all-order, although its proof uses
a finite complete quotient census. Among strong tournaments with
`H=7p<19467`, p prime and p>=11, the symmetric ones have exactly the following
counts and multiplicities. Here symmetry means a nonidentity vertex
automorphism, and classes are ordinary isomorphism classes, with converse
classes separate.

| H | p | Vertex order | Number of symmetric strong classes |
|---:|---:|---:|---:|
|4151|593|10|2|
|5789|827|10|2|
|14791|2113|12|2|
|16079|2297|12|4|
|17717|2531|12|2|

There are no other symmetric strong7p witnesses in this count range, at any
order. Every listed class has automorphism group C7. The result classifies
the correction needed to turn Hamilton-frame counts into class counts; it
does not claim new attainable values or compute the much larger asymmetric
parts of these fibres.

## 1 Inheritance and the connection to A000568

The closest proved mechanism is the Hamilton-frame fibre of
[THM-796, three-sorted recursive tiling line node incidence](../../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md):
a tournament T contributes `H(T)/|Aut(T)|` different presentations with the
fixed Hamilton spine `0->1->...->n-1`. The action on full vertex words is
free. The canonical hostile here is the symmetric H4151 class: dividing
its presentation count by H gives1/7, not one isomorphism class.

The corrected near miss is the claimed complete spectrum. The all-order
exclusion is inherited from
[THM-1370, H-spectrum omissions](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md);
[THM-4094, two-prime-lane reduction](../../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
leaves the prime and7p existence lanes open. Fixed-H rational series and
core-order membership bounds are also inherited, not new:
[Session24 multiplicities](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
and [joint spectrum](../../tournament_joint_parity_20261102/spectrum.md).

Most importantly,
[Session142, prime-lane asymmetry](../../session142_tournament_odd_spectrum_20261102/minimal_gap/PROOF.md)
already proves all of the following:

* Every prime-H witness with H>3 is asymmetric.
* If H=7p, p>=11 prime, its automorphism group has order1 or7.
* The first possible symmetric7p value is4151=7*593, and an actual witness
  attains it.
* An order7 automorphism with a nontrivial orbits forces `H>=175^a`, by
  the exact7-circulant minimum and disjoint induced path-product injection.

Those facts and their proof credit are retained. The new work is the exact
counting specialization, the complete marked-kernel comparison through
quotient order7, and the resulting all-order five-fibre classification.
The least-used coordinate is the **number of visits to a named substitution
module**, retained before taking a scalar H count.

The concept board is Hamilton frames / automorphism corrections / marked
substitution kernels / strong-deletion thresholds / actual target support.
Its anchor is the prime and7p inventory, its niche the weighted marked
strong floor, and its wildcard the possibility that one sharp two-visit
inequality controls the full weighted floor.

## 2 A class-count formula adapted to the two existence lanes

Let `N_n(h)` count ordinary n-vertex classes with H=h, and let `B_n(h)` count
Boolean tournaments with fixed Hamilton spine and H=h. No strongness filter
is imposed unless explicitly added. Then

    B_n(h)/h = sum_[T],H(T)=h 1/|Aut(T)|.                  (1)

Let `Fix_(n,h)(g)` count labelled n-vertex tournaments with H=h fixed by the
vertex permutation g. Burnside gives the exact class count. Its identity
term is `(1)`; a double count of pairs (labelled T, Hamilton path) also gives
`B_n(h)/h = (# labelled H=h tournaments)/n!`.

For a prime p>3 the inherited asymmetry theorem therefore gives

    N_n(p) = B_n(p)/p.                                   (2)

For h=7p, p>=11 prime, the only possible nonidentity automorphisms have order7.
Every such permutation has cycle type `7^a 1^b`, where b=n-7a. Consequently

    N_n(7p) = B_n(7p)/(7p)
              + sum_(1<=a<=n/7)
                  Fix_(n,7p)(g_(7^a1^b)) / (7^a a! b!). (3)

Terms with `175^a>7p` vanish. Equations(1)--(3) remain valid after restricting
both sides to strong tournaments, because strongness is isomorphism-invariant.
They are specialized applications of inherited orbit counting, not a new
Burnside theorem. At fixed h, they can avoid canonical labelling entirely;
their computational cost is still exponential in the number of free arcs.
Summing `N_n(h)` over h yields A000568(n), so this is an exact refinement of
that sequence rather than a conjectural association of two unrelated counts.

A useful alternate form is

    N_n(7p) = B_n(7p)/(7p) + (6/7) S_n(7p),             (4)

where S counts its symmetric classes. Our five-fibre classification computes
this correction completely for7p<19467. For a row of the table with core
order m and multiplicity c, ordered-SCC uniqueness gives

    S_n(7p) = c*max(n-m+1,0),
    sum_n S_n(7p) z^n = c*z^m/(1-z)^2.                 (5)

There is one nontrivial strong core, since a factorization into two
nontrivial SCC path counts would require the forbidden factor7. The extra
vertices are uniquely split into a source chain and a sink chain. This is
an application of the inherited SCC series, now with the new exact symmetric
atom inventory as input.

The loss ledger is explicit: forgetting the Hamilton frame retains H but
changes weights by |Aut|; forgetting the entire graph in favour of H loses
both the automorphism type and future attachment responses. Equations(3)
and(4) restore the first loss only. They do not construct a target witness.

## 3 The single module kernel and a coefficientwise insertion law

Let Q be a tournament with a named vertex v, and M any nonempty tournament.
Write `Q[M at v]` for uniform substitution. Define C_k(M) as the number of
ordered k-tuples of nonempty directed paths partitioning V(M), and define
A_k(Q,v) as the number of directed vertex words that visit v exactly k times
and every other quotient vertex exactly once. Consecutive repeated v symbols
are forbidden. The first coefficient is `A_1=H(Q)`.

The inherited maximal-run substitution law,
[THM-3121, path-cover walk-content kernel](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md),
specializes to

    H(Q[M at v]) = sum_(k>=1) A_k(Q,v) C_k(M).           (6)

Indeed, split a Hamilton path at changes between M and exterior vertices.
The M-runs form an ordered cover, and their contracted word is counted by
A_k. Conversely those objects concatenate uniquely. The implementation
computes A by subset DP on the exterior vertices with the previous vertex
and the visit number retained. No enumeration of the full substituted graph
is needed to evaluate(6).

**Coefficientwise insertion lemma.** If R is an induced subtournament of Q
containing v, then `A_k(R,v)<=A_k(Q,v)` for every k>=0 (with A_0 denoting the
Hamilton count on the exterior vertices). To insert a new exterior vertex x
into an old word, prepend x if it beats the first symbol, append it if the
last symbol beats x, or choose the first consecutive symbols a,b with
`a->x->b`. This always exists by the binary transition argument, including
when v occurs repeatedly. Deleting x recovers the old word and does not
change the visit number. Iterate over missing vertices. This is an injection,
not merely a numerical monotonicity test.

The map retains the actual named module and all visit coefficients. It loses
the quotient's isomorphism class and internal attachment geometry. The
positive coefficient law is sufficient for the lower bounds below; it does
not make A a complete observer for arbitrary native resets or port changes.

## 4 A finite complete basis gives universal weighted strong floors

Every strong tournament Q of order q>=4 has a strong induced subquotient on
q-1 vertices containing any specified v. Here is an elementary adaptation
of the inherited strong-deletion proof. First v lies in a directed triangle:
there is an edge from its out-neighbour set to its in-neighbour set, because
otherwise the out-neighbour set could not leave. Choose a largest proper
strong set S containing v. If an exterior vertex sees S in both directions,
it could be added. Otherwise all exterior vertices dominate S or lose to S.
Strongness supplies a pair a,b with `a->S->b->a`. There cannot be three or
more exterior vertices, since adding a,b would give a larger proper strong
set. If exactly two remain, delete any s in S different from v; the remaining
vertices lie on triangles `a->x->b->a` and form a strong(q-1)-set containing v.
Thus repeated deletion reaches any prescribed order r between3 and q.

Every rooted strong class also has a Hamilton path starting at v. One may
use the classical Hamilton-cycle theorem; the elementary proof is short.
Take a longest directed cycle. A mixed exterior vertex inserts at a cyclic
sign transition. If all exterior vertices are uniform to the cycle, strongness
supplies an exterior pair `cycle->a->b->cycle`, which also extends it.
Hence the cycle is spanning, and cutting immediately before v gives the path.
Therefore enumerating the `2^((q-1)(q-2)/2)` fixed-spine graphs with v=0 covers
**every rooted strong class**, without any assumption about its symmetry.

The complete scans atq=3..7 yield the following componentwise minimal vectors
`(A_1,A_2,...)`; all omitted coordinates are zero. The receipt gives a literal
fixed-spine source for each row. Each other kernel in the universe dominates
at least one displayed vector coordinate by coordinate.

| Quotient order | Minimal kernel vectors | Fixed-spine presentations tested |
|---:|---|---:|
|3|(3,1)|2|
|4|(5,2)|8|
|5|(9,3)|64|
|6|(15,5)|1024|
|7|(25,10), (27,9)|32768|

Together with(6), this proves exact weighted minima over all rooted strong
quotients of each listed order for **every** nonempty M, of arbitrary order.
In particular atq7 the exact minimum is

    min(25 H(M)+10 C_2(M), 27 H(M)+9 C_2(M)).            (7)

For |M|>=2, `C_2(M)>=2H(M)`. At order2 this is an immediate check. At order
m>=3, cut each Hamilton path at any of its m-1 internal positions; the
resulting ordered covers are distinct, so `C_2>=(m-1)H>=2H`.
Hence(7) simplifies, for every nontrivial M, to

    9(3H(M)+C_2(M)) = 9(H(M)+F_M(2)),                 (8)

where `F_M(2)=2H(M)+C_2(M)` allows empty slots. This bound remains valid for
**every strong quotient of order at least7**, by marked strong deletion and
the coefficientwise insertion lemma. It is attained at quotient order7.
This is a proved all-module/all-larger-quotient statement with a finite exact
basis; it does not assert a closed formula for the new exact minima atq>=8.

## 5 Complete symmetric seven prime inventory below19467

Let T be strong, symmetric, and H(T)=7p<19467 with p>=11 prime. Its group is
C7 by the inherited theorem. Two nontrivial orbits would force
`H>=175^2=30625`, so there is exactly one7-orbit M, and all other vertices
are fixed. They see M uniformly, making it a module. Contract it to v in a
strong rooted quotient Q.

The complete eight-presentation circulant table has six presentations with
`H=175, |Aut|=7` and two with `H=189, |Aut|=21`. The six are one ordinary
class, represented by forward steps{1,2,3} on Z/7. The latter option cannot
be M: every automorphism of a uniform module extends while fixing its
exterior, so Aut(T) would contain a group of order21. Thus M is the fixed
175-path class. Its exact ordered cover vector is

    (C_1,...,C_7)=(175,1638,6678,15120,20160,15120,5040).

Independent permutation-and-cut enumeration verifies every entry. By(8),
`|Q|>=7` would give

    H(T)>=9(3*175+1638)=19467,

contrary to the chosen strict bound. The quotient cannot have order2 if T
is strong; order1 gives175, whose quotient25 is not prime; order3 gives2163,
whose quotient309 is composite. The complete rooted quotient enumeration
at orders4,5,6 then gives precisely the five rows in the opening table.

The class multiplicities have two independent implementations. First,
rooted-first Hamilton frames contribute weight `1/s(Q,v)`, where s is the
number of Hamilton paths starting at v. In general this totals the reciprocal
rooted automorphism order. In this filtered family the rooted automorphism
group is trivial: it would otherwise extend to Aut(T) beyond the already
present C7. Thus the sum counts rooted quotient classes exactly. Second,
enumerate all labelled quotients on q vertices with v=0 and divide by(q-1)!,
which is legitimate in this same rigid rooted family. Successful sources are
also checked directly under all root-preserving permutations.

This quotient classification equals the T-class classification. The unique
nontrivial C7 orbit is characteristic under every isomorphism, so T determines
its marked quotient. Conversely a rooted quotient isomorphism extends over
the single fixed module class. The checker records actual full tournaments
for every listed H and independently recounts their paths with two vertex
DP algorithms.

As a hostile control for using only the identity term in(3), the two strong
H4151 classes contribute2/7 to that term. A fixed7-cycle with three fixed
vertices leaves exactly72 labelled H4151 tournaments invariant. Its Burnside
correction is `72/(7*3!)=12/7`, restoring the exact symmetric contribution2.
A uniform division by H would lose this contribution.

## 6 The sharper open lemma and proof frontier

The finite kernel basis suggests the explicit statement

    A_2(Q,v) >= f(q-2)  for every strong rooted Q, q>=5, (OPEN)

where f is the exact Busch minimum for strong tournaments. It holds for the
complete bases throughq7. If it holds in general, apply Busch to the actual
transitive-two-clone substitution Q[TT2], giving
`A_1+2A_2>=f(q+1)=5f(q-2)`. For every |M|>=2,

    H(Q[M]) >= H(M)(A_1+2A_2)+(C_2(M)-2H(M))A_2
             >= f(q-2)(3H(M)+C_2(M)).                (CONDITIONAL)

This is a precise positive-response obligation, not a proof by fitting the
five minima. No all-order assertion about this new two-visit bound is used
in the five-fibre classification. Nor would the bound itself establish
arbitrary target realization: it constrains how cheaply a marked module can
be placed inside a strong graph. Target selection and reusable attachment
states remain separate obligations.

## 7 Relation to the concurrent cyclic module and clone results

The incoming [Session175 arbitrary cyclic-module analysis](../../session175_tournament_spectrum_20261104/obstructions/PROOF.md)
uses the same inherited THM-3121 positive kernel and completely classifies
its constructor image through45. Its modules vary inside a fixed C3
quotient; here the rooted strong quotient varies around one fixed module.
Thus its fixed-size transitive minimum and our marked-quotient minimum
address different parameters of the same substitution law. The exact
low-spectrum exclusion supplies a hostile to treating uniform cyclic
composition as a universal constructor.

The [finite clone-library obstruction](../../session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md)
shows that finitely many fixed quotient templates, inflated only by
transitive modules, supply only polylogarithmically many prime counts below
a growing bound. Our formula(3) is an exact enumeration refinement, not a
finite-template realization strategy, so there is no conflict. A successful
prime-lane construction still needs growing interface structure or modules
with additional arithmetic information.

## 8 Reproduction and scope

From the repository root:

```sh
python3 -B 05-knowledge/results/orbit_spectrum_session_20261103/bridge/check.py > /tmp/bridge-check.json
python3 -B -O 05-knowledge/results/orbit_spectrum_session_20261103/bridge/check.py > /tmp/bridge-check-optimized.json
cmp /tmp/bridge-check.json /tmp/bridge-check-optimized.json
```

The program is standalone standard-library Python. It exhausts all Hamilton
frames through quotient order7, all labelled rooted quotients through6, all
eight7-circulants, and all permutations/cuts of the selected7-module. It checks
72 quotient kernels by an independent permutation/gap algorithm,4732 induced
coefficientwise comparisons, the five full target witnesses and five minimum
witnesses by independent Hamilton DP, the rooted class multiplicities by
literal canonicalization, and the prime-H frame formula against labelled
counts through6. There are no symmetry or target filters before the stated
universe generation; filters are applied explicitly afterward. Checks raise
exceptions and remain active under optimization.

The frozen normal and optimized receipts agree byte for byte. The source
script hash is inside the receipt; `sources.json` records inherited source
paths and Git/blob hashes. No all-order external novelty claim is made.
