# Coupled residual moves and shared-frame compatibility

23 October 2026. **PROVED, typed all-size criteria; FINITE-EXACT controls.**

This packet extends Session38's transport-cycle uniqueness criterion to a
finite integer attachment model with shared-owner constraints. It does not
claim that every anonymous metagraph wiring satisfying counts is native.
It also gives the correct frame-comparison problem for odd-cell blocks:
nonexistence of a parallel owner section does not obstruct the existence or
classification of the tournament.

## Inheritance and contribution

- Session38 `PROJECTED_TRANSPORT_UNIQUENESS.md`: exact residual-SCC
  potentials for a bounded integral network transport fiber. Its atomic
  moves are single directed residual cycles. Session38's independent audit
  explicitly retains reflection equalities rather than dropping them.
- Session39 `AFFINE_CONNECTIONS.md`: actual cyclic-triple block decoding,
  jointly based affine holonomy, and distinction between sections and
  frames. THM-2622 and THM-2089 are earlier general mechanisms.
- Incoming d45865, `metagraph_owned_boundary_frontier_20261023.md` and
  `hopf_sierpinski_port_recovery_20261023/REPORT.md`: all 420 reached
  order-31 eight-triple systems have no parallel section, but rotational
  centralizers distinguish the 252 C3 cases from the 168 S3 cases.
- Session40 `PATH_TRANSPORT_AND_FOREST_RULE.md`: complete actual path
  transport and the at-most-six-unmarked positive rule; target class is not
  literal port identity.
- Canon THM-4009 and THM-4014 already use conformal/Graver decomposition in
  arithmetic relation compression. We transfer that standard mechanism to
  the typed attachment fiber, not claim a new general integer-programming
  theorem.

**New increment:** a finite bounded atomic-move iff criterion including all
shared-port rows; a two-switch hostile showing why individually admissible
cycles are insufficient; an explicit finite-domain compilation retaining
one common cell frame; and a correct arbitrary-block frame-comparison
certificate. These distinguish integer move coupling, additive potentials,
nonabelian holonomy, and actual-owner recovery.

## 1. Exact bounded coupled-move theorem

Let A be any integer matrix, l,u integer vectors with finite l<=u, and

    F = {x in Z^d : Ax=b, l<=x<=u}.

A contains **all** linearized obligations: transport margins, actual shared
owner identifications, reflection matching, one-hot choices, and any other
constraints used in this declared model. Let P be the integer observation
matrix whose uniqueness is asked. Auxiliary coordinates are allowed; P
simply ignores them when appropriate. Let D=u-l.

For vectors g,h write h <=c g if every coordinate has the same weak sign
and |h_i|<=|g_i|. Define the bounded atomic bank

    G_D(A) = nonzero g in ker_Z(A), |g|<=D,
             with no nonzero proper h in ker_Z(A) satisfying h <=c g.

This bank is finite by the displayed bounds. It is the part of the usual
conformally primitive kernel bank that fits the declared coordinate widths.
No unproved finite generation assertion is needed for this definition.

**Theorem 1.** For any feasible integer x, the following are equivalent:

1. P is constant on F.
2. Pg=0 for every g in G_D(A) for which l<=x+g<=u.

**Proof.** If condition 1 holds, every feasible atomic augmentation leaves
P unchanged. Conversely let y be any feasible point and d=y-x. Repeatedly
choose a nonzero conformally minimal kernel subvector of the remaining d.
Such a vector exists because the finite conformal interval contains d; it
belongs to G_D(A). Subtracting it leaves a kernel vector in the same orthant,
and decreases the integer l1 norm. Thus d is a finite conformal sum of
atoms. Every atom fits coordinatewise between 0 and d, so x+g lies between
x and y and is feasible. Condition 2 annihilates every atom, hence P(y-x)=0.
This also proves that a failed uniqueness certificate has a **single
feasible coupled atom** as an exact witness. QED.

A purported list of atoms is not a certificate of completeness without a
coverage proof. A slow independent verifier can enumerate the finite box
[-D,D], check every kernel vector, and identify its conformal minimality.
A faster generator must prove the same coverage. The theorem is all-size;
it does not promise a uniformly small bank or efficient enumeration.

### Why Session38 is the network specialization

When A is a signed network incidence matrix, conformal primitive moves are
signed simple cycles (with distinct parallel channels retained). Bounds
select precisely the positive-capacity residual cycles. Hence Theorem 1
recovers the Session38 result. Once shared-port rows are added, primitives
can combine several former cycles. Deleting those rows silently changes the
question. Keeping them but testing only individual cycles is also wrong.

## 2. Small exact hostile: two cells must switch together

Take two 2-by-2 permutation tables X,Y. Every row and column sums to one;
all eight entries are binary. Add the shared-port equation

    X_00 = Y_00.

The feasible set is exactly (I,I) and (S,S), where S swaps the two labels.
At x=(I,I), neither the cycle switching X alone nor that switching Y alone
is allowed by the shared equation. Nevertheless switching both is allowed
and changes P=X_00 from one to zero. The coupled atom is

    (-1,+1,+1,-1 ; -1,+1,+1,-1).

Its l1 norm is eight, twice that of each uncoupled elementary switch. The
bank has exactly this vector and its negative. Thus "every individually
admissible ordinary residual cycle has zero observed cost" is vacuously
true and yet projected uniqueness is false. This is an abstract attachment
hostile, not an asserted new native tournament counterfeit.

The opposite relaxation error is equally concrete. Restrict one 2-by-2
permutation table by X_00=1. The declared fiber is a singleton. Its relaxed
transport graph still contains a nonzero-cost switch. Therefore failure of
an uncoupled potential test does not imply failure in the coupled fiber.
An uncoupled zero-cost potential **does** remain sufficient for any subset.

### Integer widths matter, not only rational cycle directions

The bounded fiber 2a+b=1 with 0<=a,b<=1 has the unique integer point (0,1),
but also rational points (t,1-2t) for 0<t<=1/2. The primitive kernel move
(1,-2) cannot fit the b width. Real linear dependence or a rational-flow
ambiguity is therefore not an integer attachment witness. This example is
not itself a tournament encoding; it isolates the forbidden inference.

## 3. A lossless finite-domain shared-frame compiler

For each cell i let D_i be a finite set of permissible **whole cell frames
or whole local attachment states**. Every actual marked owner must be fixed
by a permissible frame. For each interacting pair ij retain an explicit
compatibility relation R_ij subset D_i x D_j. Keep the same variable for a
cell wherever that cell occurs; independent copies would recycle ownership.

Introduce binary z_(i,a), one for each a in D_i, and binary w_(ij,a,b) for
each (a,b) in R_ij. Impose

    sum_a z_(i,a)=1;
    sum_b w_(ij,a,b)=z_(i,a) for every a;
    sum_a w_(ij,a,b)=z_(j,b) for every b.

An empty sum is zero. At an integer solution each cell selects one state,
and every edge selects exactly that pair; hence the pair must be allowed.
Conversely any globally compatible state assignment gives exactly one
integer solution. This is a bijection, not just a necessary margin test.
Every linear port observation can now be tested by Theorem 1. More general
finite observations can be compiled as retained one-hot output states.
For a fixed linear electrical guard Cx=d (mod 8), add integer variables z
and exact rows Cx-8z=d. For each row, interval arithmetic on l<=x<=u gives
finite integer bounds ceil((min Cx-d)/8)<=z<=floor((max Cx-d)/8). These
auxiliary quotients are retained throughout atom generation. Actual vertex
labels remain fixed. A nonlinear linking invariant is not linear merely
because its output is reduced modulo eight.

Higher-arity constraints require factor variables indexed by **whole
allowed tuples**, linked to the same z variables; pair projections alone
are not generally sufficient. Congruence guards can likewise be compiled
by exact finite relations, rather than being confused with real equalities.

**Forest corollary.** On a forest of pair constraints, every nonnegative
real solution of these equations is a convex combination of compatible
integer assignments. To prove it, choose a root according to its z
marginal and, along each oriented edge, sample the child's state with
conditional probability w_(ij,a,b)/z_(i,a). Zero-probability parent states
are never used. This reproduces all node and edge marginals. Different
components can be sampled independently. Thus real feasibility and linear
extrema are exact on this carrier. This does not assert that a tournament's
actual interaction graph is a forest.

**Cycle hostile.** On a triangle, let every D_i={0,1} and require unequal
states on every edge. There is no integer assignment. Nevertheless each
z value 1/2 and each allowed edge-pair w value 1/2 satisfies every displayed
real equation. All local relations are nonempty and every one-cell marginal
agrees. The obstruction is the common owner on the closed cycle. Exhaustive
controls cover all 16^3 binary relations on a triangle.

## 4. Frame comparison, not mandatory parallel ownership

Let A and B be two fully specified tournaments with the same named cell
partition. Retain the internal tournaments and first choose a reference
internal isomorphism per cell that sends every retained source mark to its
prescribed target mark (reject if no such mark-preserving reference exists).
After using these
references to identify the cell sets, put D_i=Aut(A[C_i]), restricted to
fix actual retained marks. Define R_ij to consist of pairs (g_i,g_j) with

    B(g_i(u),g_j(v)) = A(u,v) for every u in C_i, v in C_j.

**Theorem 2.** Globally compatible assignments in this finite-domain model
are in bijection with cell-preserving isomorphisms A->B fixing the retained
marks (after the reference identifications).

The proof is direct: internal arcs are preserved by D_i, cross arcs by
R_ij, and every arc is one of these two kinds. Conversely every such
isomorphism restricts to exactly these frames. The carrier is valid for
arbitrary odd cells, not just permutation blocks, and also for even cells
when their internal tournaments are supplied. No small-observer or fast
algorithm claim follows: supplying all blocks already supplies the arcs.
This is a rigorous baseline against which a proposed compressed menu must
prove losslessness.

This model chooses a coordinate identification, not a vertex which must be
transported unchanged around every cycle. It always has the identity
assignment when A=B, including valid tournaments with no parallel owner
section. Applying the section-existence test instead would incorrectly
reject all 420 incoming order-31 eight-triple examples.

### Bijection blocks: simultaneous conjugacy certificate

Suppose each edge block really selects a bijection alpha_ij for A and
beta_ij for B, in one common group G acting on the cell labels. On a
connected support pick a spanning tree and root r. Let a_i,b_i be the
respective root-to-i tree transports. Any comparison must satisfy

    g_i = b_i g_r a_i^(-1).

For each chord ij form H_A=a_j^(-1) alpha_ij a_i and
H_B=b_j^(-1) beta_ij b_i. The necessary and sufficient remaining conditions
are, **with one common g_r**,

    H_B = g_r H_A g_r^(-1),
    g_i belongs to the actual permitted D_i for every i.

This proves an all-size constructive certificate: one root frame and all
chord equations, with membership checked at every cell. Edgewise or
cyclewise unrelated conjugators do not suffice. Uniform/polarity/internal
arc data remain separately retained. If A=B and all D_i=G, solutions are
exactly the centralizer in G of the subgroup generated by all H_A. With
restricted D_i they are the subgroup of that centralizer passing all
transported membership tests.

### Exact order-31 lesson in the triple carrier

For G=Aff(F3)=S3, allowed root charts are rotations C3. Write t(x)=x+1 and
s(x)=-x. A connection with holonomy subgroup <t>=C3 has no parallel owner
section, but all three root rotations commute with its holonomy. A
connection with holonomy subgroup <t,s>=S3 also has no parallel owner
section, but only the identity rotation commutes with its holonomy.
Thus these two no-section situations have respectively three and one
compatible self-comparison frames. This is precisely the distinction used
by the incoming 252/168 split, without independently re-enumerating that
census here. A single marked actual owner kills a nontrivial root rotation.
For seven-cell Paley carriers the allowed group and matching relation must
be changed; a triple-carrier conclusion is not automatically a seven-cell
conclusion.

## 5. Additive potentials and nonabelian holonomy are different

An additive cost potential h gives c_ij=h_j-h_i; its sum around every closed
walk is zero. A group-valued frame gives transports g_j g_i^(-1); their
ordered product around every closed walk is identity. The first is an
abelian linear statement. Taking an abelianized shadow of the second can
lose exactly the obstruction one wants to detect.

Concrete S3 hostile: t=(012), s=(12), and sts=t^(-1). The commutator

    t s t^(-1) s^(-1) = t^2

is a nonidentity fixed-point-free rotation, while its image in the sign
quotient S3->C2 is identity. Both translation holonomy and identity have
positive sign, but different section sets. Ordinary scalar electrical
potentials cannot recover the discarded label action. In the special
translation subgroup C3, additive F3 potentials are appropriate, but
integer/mod-8 potentials are not a replacement for them.

No result here discards converse, the electrical mod-8 guard, integer
clocks, or literal path ownership. Those remain independent constraints
in the actual application. No implication to Collatz coverage follows.

## Scope and remaining obligation

The all-size conclusions are exact **for a supplied complete finite
attachment model**, and for cell-preserving frame comparison with supplied
blocks. To classify complete native tournament metagraph wiring one still
must prove that the chosen states, factors, and observations cover exactly
all native source/path/target records. Counts alone do not provide that
exhaustion. The practical next target is to compile surviving port choices
using actual shared owners, then seek a small structured atomic bank or a
real forest/tree decomposition, rather than hoping that uncoupled cycles
or empty parallel-section summaries settle the problem.

## Reproducible finite controls

Run `python3 check.py` and `python3 -O check.py`. The two receipts agree.
The checker exhausts 250 declared small integer matrices, 33,750 basepoint
observation tests, and all 4,096 binary pair-relation systems on a triangle.
It also checks the coupled-switch and integer-width hostiles, exact mod-8
auxiliary encoding, centralizers, and shared-root conjugacy. There are 54
ordered holonomy-pair comparisons in the tested S3 universe that pass
separate rotation conjugacy but fail simultaneous rotation conjugacy.
These finite controls support, rather than replace, the displayed proofs.

The executable `compile_attachments.py` builds the actual bounded integer
matrix, with `encode`/`decode` and a full-block frame-relation constructor.
`check_compiler.py` independently compares 32,768 assignment encodings and
324 actual six-vertex tournament frame comparisons. Frame domains must
already retain actual marked owners; the constructor does not guess marks.
