# Signed converse transport controls the lift of class symmetries

**Status: PROVED IN NOTES after independent proof/code audit; see [AUDIT.md](../AUDIT.md).**
The elementary statements S1–S4 are proved all-order operator statements;
S5 is **FINITE-EXACT** for orders3–7. No scarce theorem ID is reserved.
No preferred fixed-spine presentation lift is asserted for an arbitrary
class-matrix symmetry. Global converse has its inherited literal reflection
lift. Chan–Pak and LRC(14) remain **OPEN**.

## Inheritance and board

The closest proved mechanism is the [general converse-orbit lift](../../joint_metagraph_dyadic_20261113/PAIRING_WIRING.md):
if only the merged diagram is retained, its compatible ordinary-class
permutations form `C2^t semidirect G`, and the prior involution polynomial
counts them exactly. That theorem and its [independent review](../../joint_metagraph_dyadic_20261113/metagraph/PAIRING_LIFT_REVIEW.md)
are inherited, not new here. Its exact order-five result has two independent
chiral flips and polynomial `x^8(1+x^2)^2`.

The canonical hostile is the order-five full-Q pair, codes8 and10, with
identical path-cover profiles but different transitive-node incidence;
see [full-slot wiring](../../joint_metagraph_dyadic_20261113/SLOT_WIRING.md).
The corrected near miss is treating independent flips inside converse pairs
as already equipped with a coherent action on actual path presentations.
The least-used sidecar pursued here is the **difference** between two
ordinary-class transport channels that merge to the same edge.

Portfolio: anchor = exact incidence-compatible pairing; niche = signed
converse eigenspaces; wildcard = cycle-index consequences. The live board
contains merged wiring, ordinary wiring, the signed difference operator,
component-wise flips, and cycle/Burnside counts. All source conventions and
inherited-file fingerprints are in [source-pins.json](source-pins.json).

## Objects and normalization

Fix n>=3. Let X be ordinary tournament classes, with converse R. Let
`Y=X/R`, with s singleton fibres and t double fibres. Thus `|X|=s+2t` and
`|Y|=s+t`. A representative T and a directed Hamilton vertex order P define
the native move that **complements all nonpath arcs**, retaining P.

Let C(T,U) be the number of Aut(T)-orbits of Hamilton paths whose native
partner is of ordinary class U. Equivalently it is the raw path count
restricted to that partner class, divided by |Aut(T)|. The action on full
path orders is free. The resulting C is a symmetric nonnegative integer
matrix, commutes with R, and has row sum `d(T)=H(T)/|Aut(T)|`. Parallel
lines and doubled loop entries are retained. Define its blue part C_B by
the inherited reflection-fixed-path criterion and let C_K=C-C_B.

The merged matrix is the **sum over both endpoint fibres**:

```text
A(u,v)=sum_(x in u,y in v) C(x,y).                    (1)
```

It is not one selected row. If m(u)=|u| is1 or2, its row degree is
`d_Y(u)=m(u)d(x)`. Both C and A retain literal dart multiplicities.
Blue darts have self-converse source and target, so C_B is supported on
`Y_s x Y_s`; the converse action commutes with the native partner map.
Consequently all signed chiral transport below is black transport.

## S1. The missing signed matrix reconstructs ordinary wiring

Choose one representative x_u in each double fibre. Define

```text
D(u,v)=C(x_u,x_v)-C(x_u,Rx_v),      u,v in Y_c.        (2)
```

This is symmetric. Changing representative choices by signs epsilon_u
changes D to `diag(epsilon) D diag(epsilon)`; it leaves A unchanged.
In particular its off-diagonal support graph is independent of the choices.

The data `(A, fibre sizes, chosen D)` reconstruct C exactly. On a pair of
double fibres its block is

```text
[[a,b],[b,a]],       a+b=A(u,v)/2,       a-b=D(u,v),
a=A(u,v)/4+D(u,v)/2, b=A(u,v)/4-D(u,v)/2.             (3)
```

On a singleton–double block both entries are A(u,v)/2; on two singleton
fibres the entry is A(u,v). Thus double-double data must satisfy
`A(u,v)/2 = D(u,v) (mod2)` and `|D(u,v)|<=A(u,v)/2`; mixed blocks require
A even. Diagonal C entries also retain their actual loop convention.
These are reconstruction guards, not a characterization of all formally
possible pairs `(A,D)` as tournament incidence matrices.

*Proof.* R-equivariance forces the displayed block forms; sum and difference
solve for a,b. Symmetry and gauge covariance follow directly. ∎

The source is an ordinary-class incidence operator; the quotient map loses
a signed difference per pair of chiral fibres. This is a genuine information
loss, even when the entire merged matrix is known.

## S2. The exact kernel is constant on signed-support components

Let Sigma have vertex set Y_c and an edge u–v for distinct u,v exactly
when D(u,v) is nonzero. Include isolated chiral fibres as components. Let
its components be Sigma_1,...,Sigma_c, of sizes m_1,...,m_c.

A permutation fixing every merged vertex and commuting with R independently
chooses whether to exchange each double fibre. Such a permutation preserves
C **if and only if its choice is constant on every component of Sigma**.
The kernel of the action on merged vertices is therefore exactly `C2^c`,
not automatically `C2^t`. The same kernel preserves blue and black matrices.

*Proof.* Flipping fibre u multiplies its signed coordinate by epsilon_u.
By S1 the C-preservation condition is exactly
`epsilon_u epsilon_v D(u,v)=D(u,v)` for all u,v. For an off-diagonal nonzero
entry this says epsilon_u=epsilon_v. Diagonal entries impose no condition.
These equalities are precisely constancy on connected components. Mixed
and singleton blocks are unchanged, and blue support is singleton-only. ∎

**Conditional all-order consequence.** If the typed merged weighted graph
is vertex-rigid, its entire R-commuting ordinary-class matrix automorphism
group is `C2^c`. Every member is an involution, and its fixed-point polynomial is

```text
I_C(x)=x^s product_(i=1..c) (1+x^(2m_i)).             (4)
```

The minimum fixed count is s and its unique minimizer is R, including empty
chiral sector. If in addition Sigma is connected and nonempty, the group
is exactly `{identity,R}` and (4) is `x^s+x^(s+2t)`.

This is a statement about the class dart-adjacency matrix. A matrix symmetry
can be lifted to arbitrary bijections of parallel darts, but that does not
supply the particular fixed-spine transformations or coherent marked-path
operation demanded by a stronger presentation action.

## Which base symmetries lift

**S3.** Let g preserve the fibre sectors and A (and whichever extra labels
or color matrices are required). An ordinary lift has signs epsilon_u on
the chiral fibres. It preserves C exactly when

```text
D(g(u),g(v))=epsilon_u epsilon_v D(u,v) for all u,v.   (5)
```

Thus g must preserve zero positions, absolute values and the diagonal of D.
For each nonzero off-diagonal edge, the ratio
`r_uv=D(g(u),g(v))/D(u,v)` must be ±1. A lift exists exactly when the product
of these ratios around **every support cycle** is+1. If it exists there
are exactly2^c sign choices, one overall sign per component.

*Proof.* Equation(5) is S1 after permutation and gauge transport. Necessity
of the cycle test follows by multiplying vertex signs twice around a cycle.
For sufficiency, choose a sign at one root of each component and propagate
along a spanning tree; every nontree edge is consistent exactly when its
fundamental-cycle product is+1. Loops were already handled separately. ∎

The image group G_lift can be smaller than the merged graph group, and the
exact sequence with kernel C2^c is not asserted to split. Raw nonzero D
entries are integer currents, not invertible transports. Only the ±1
**ratios for a tested base permutation** form the coboundary problem.
The companion completion-transport lane uses the same equation on its own
carrier; no identification of those two graphs is assumed.

A synthetic hostile in the checker takes four double fibres, identical
merged off-diagonal weights6, and signed off-diagonal weights+1 except one
edge of weight−1. All24 base permutations preserve A, but only4 lift,
yielding8 ordinary symmetries. Changing the representative gauge supplies
positive cases requiring nonconstant signs. Brute force over all8! ordinary
permutations independently checks the complete answer.

## S4. Representation and cycle-index meaning

Over the reals the ordinary permutation space decomposes under converse as

```text
R^X = V_+ direct-sum V_-,
dim(V_+)=s+t,       dim(V_-)=t.                      (6)
```

V_+ consists of functions constant on each R-fibre. With the orthonormal
basis of singleton indicators and normalized double-fibre sums, C acts by
`M^(-1/2) A M^(-1/2)`, where M=diag(m(u)). On the normalized differences
`(e_xu-e_Rxu)/sqrt2`, C acts by D. Hence

```text
det(lambda I-C)
 = det(lambda I-M^(-1/2) A M^(-1/2)) det(lambda I-D). (7)
```

For random-walk matrices, the plus block in the natural fibre-constant
basis is the merged stochastic operator `P_Y=diag(d_Y)^(-1)A`.
In the orthonormal plus basis above it is instead
`M^(1/2) P_Y M^(-1/2)`. The minus block on normalized differences is
`diag(d(x_u))^(-1)D`. These weights and basis choices are essential:
A itself is not the plus block in an orthonormal basis.

For the kernel C2^c, V_+ is fixed pointwise. The minus space splits into
c real sign-character spaces of dimensions m_i, one for each signed-support
component. Therefore an independent sign cannot be changed inside one
nonzero-current component. This is a literal representation of commuting
involutions. It is not an anticommuting Clifford-generator construction,
and implies no Bott or eightfold periodicity.

When the base graph is rigid, the exact cycle index of this ordinary-class
action is

```text
Z_C(x_1,x_2,...)=x_1^s / 2^c
                 product_i (x_1^(2m_i)+x_2^(m_i)).  (8)
```

A component left unflipped contributes2m_i fixed elements; a flipped
component contributes m_i transpositions. This proves both(8) and(4).
Burnside substitution x_j=q counts q-colorings of **ordinary classes**
modulo this matrix action, giving
`q^s product_i(q^(2m_i)+q^(m_i))/2^c`.

For a nontrivial base group the same cycle-index construction is available
only after solving(5): sum over liftable g and their valid epsilon solutions.
A singleton g-cycle of length l contributes x_l; a chiral g-cycle contributes
x_l^2 when the product of its transport signs is+1, and x_(2l) otherwise.
Normalize by the actual number of lifts. This also counts only the valid
involutions after imposing g^2=1 and epsilon_u epsilon_g(u)=1.

These formulas extend the inherited lift bookkeeping to the extra signed
constraint. They do not derive A000568 from prior scalar counts: s,t and
especially the component sizes or liftable group remain inputs from actual
incidence. The [pointed-cycle Burnside recurrence](../../session315_dyadic_pairing_20261113/SPECTRUM_AND_SIEVE.md)
and the prior gcd-cycle inventory remain the existing count recurrences.
Replacing their retained cycle variables by class-action dimensions would
lose their necessary input, not improve that recurrence.

## S5. Exact orders3–7 and the first kernel collapse

The following results use the full weighted merged matrix, initially
without H, type or color guards. A multiset refinement starts from weighted
degree and loop weight and records the actual edge weight into each current
vertex color. A discrete final partition proves vertex rigidity. Exhaustive
backtracking in its remaining cells checks every possible permutation.

| n | ordinary | merged | SC | chiral pairs t | uncolored refinement cell counts | component sizes of Sigma |
|---:|---:|---:|---:|---:|---|---|
|3|2|2|2|0|1→1|empty|
|4|4|3|2|1|3→3|1|
|5|12|10|8|2|8→10→10|1,1|
|6|56|34|12|22|26→33→34→34|22|
|7|456|272|88|184|145→271→272→272|184|

At order3 the uncolored graph has its two-vertex swap, which changes H.
At orders4–7 the merged weighted graph is vertex-rigid without any guard.
Thus the complete **H-preserving, R-commuting ordinary C-matrix** groups
have orders `1,2,4,2,2` at n=3,4,5,6,7. Their involution polynomials are

```text
n=3: x^2;
n=4: x^2+x^4;
n=5: x^8+2x^10+x^12;
n=6: x^12+x^56;
n=7: x^88+x^456.                                    (9)
```

The same results hold with the blue/black matrix guard. The new n6,n7
minimum fixed counts are12 and88. In contrast, retaining only merged
wiring gives the inherited polynomials `x^12(1+x^2)^22` and
`x^88(1+x^2)^184`. Their minima agree, but almost all their independently
chosen pair flips fail the ordinary incidence guard.

The first failure of arbitrary independent chiral flips occurs at order6:
there are no nonzero signed off-diagonal entries through5. At6, chiral
representatives2 and18 satisfy `C(2,18)=1`, `C(2,R18)=0`, so flipping only
one changes a matrix entry. Its complete signed-support graph is connected.
A spanning tree with22 vertices and21 nonzero entries certifies that all
flip choices are forced equal. The row's
`signed_spanning_forest_C_and_D` field records every tree edge with both
ordinary C entries and their difference. At7 the corresponding certificate
has184 vertices and183 tree edges; the full signed support has3377 edges.

## Verification and retained scope

[check.py](check.py) is independent of the producer code. At n=3..6 it
partitions every labelled orientation into complete S_n orbits, then lists
every literal Hamilton permutation of every ordinary class. It computes
native partners and blue tests directly, with source automorphism weights.
Every ordinary and merged matrix, both colors, and the class metadata agree
entrywise with [closure/data.json](../closure/data.json).

At n=7 the full matrix is the producer's complete finite census, whose two
independent path generators agree. This checker independently recomputes
six entire low-H ordinary rows by literal paths and unrestricted permutation
canonicalization, including H, Aut and converse. It independently audits
the complete exported matrix's reconstruction, all weighted graph groups,
the signed support, all component-flip generators, and every actual
ordinary lift in the small resulting groups. It does not call those six
rows a second full n7 census.

The normal and optimized [receipts](normal.json) pin the producer and checker
hashes and retain all consequences, including exact spanning forests.
Reproduce from the repository root:

```sh
python3 -I -B 05-knowledge/results/incidence_closure_20261113/pairings/check.py
python3 -I -B -O 05-knowledge/results/incidence_closure_20261113/pairings/check.py
```

The conditional all-order test is now precise: establish merged vertex
rigidity and connected nonempty signed support. Neither property is proved
at arbitrary order here. The new sidecar classifies a specific matrix-lift
problem; it does not make graph connectedness, equal Q, or scalar class
counts a complete tournament-spectrum theorem.
