# A uniform strict Gold certificate in a two-interface, seven-thin family

Date: 2026-10-19. Status: **PROVED**, by the kernel and finite hand sums below; **FINITE-EXACT** independent controls. The general Gold Partition and 1/3–2/3 conjectures remain **OPEN**.

## Result

An explicit nine-element core B has two fixed comparisons that satisfy the stronger budget

    t_0 - t_1 - t_2 >= (109/1174) t_0 > 0                    (1)

in every allowed two-interface context. Concatenating any number of copies of B through those interfaces produces an unbounded family of width-three, seven-thin, ordinal-indecomposable posets. Every block, including an interior block, supplies the same strict Gold certificate, irrespective of the length or structure of the surrounding permitted context.

The construction does not prove Gold for every seven-thin or width-three poset. It supplies a context-stable certificate family beyond simply repeating a fixed finite example or applying equal-fibre module substitution.

The Gold inequality is Peczarski's original t_0>=t_1+t_2, not a strict conjecture:
https://cdm.ucalgary.ca/article/download/62349/46785/178291
The present strict margin is a property of this family. The stronger companion `gold_unbounded_flag_context_20261019.md` allows arbitrary terminal-upset and initial-ideal collars, using a different fixed strategy that also preserves 2/5 balance.

The underlying two-interface kernel and balance example are also developed in the companion `poset_four_collar_boundary_20261019.md`. This note strengthens its one-comparison balance certificate to one fixed, strict two-comparison Gold strategy.

## 1. Two-interface kernel and four-corner Gold theorem

Let B be a finite poset, a a minimal element, and b a distinct maximal element. Let A,C be arbitrary finite posets; when A is nonempty choose a maximal u in A, and when C is nonempty choose a minimal v in C.

Start with the ordinal sum A<B<C. Remove just the cross-relations u<a and b<v, making those pairs incomparable; retain all other cross-relations. Minimality/maximality ensures the resulting relation remains transitive. Empty context on either side is allowed.

For sigma in L(B), define

    I(sigma)=1{a is first in sigma},
    J(sigma)=1{b is last in sigma},
    r=Pr_A[u last], s=Pr_C[v first],

with r=0 or s=0 for the respective empty context. Then the exact restriction fibre, divided by e(A)e(C), is

    w_{r,s}(sigma)=(1+r I(sigma))(1+s J(sigma)).              (2)

**Proof.** After fixing orders on A,B,C, the only possible nonordinal interleavings swap u with a when u is last and a first, and swap b with v when b is last and v first. Each eligible swap is optional and the swaps are independent, because a and b are distinct. All extensions arise uniquely in this way. Summing over the exterior orders gives (2). ∎

Fix two core comparisons, and let X_0,X_a,X_b,X_{ab} be their four-cell count matrices on all extensions, extensions starting a, extensions ending b, and extensions having both flags. The context matrix is

    X(r,s)=X_0+r X_a+s X_b+rs X_{ab}.                         (3)

For a table X=(x_{ij}), the four Gold deficits are

    D_{ij}(X)=sum_j x_{1-i,j}-x_{ij}.                         (4)

They are exactly t_0-t_1-t_2 for the corresponding path.

**Four-corner theorem.** A fixed strategy satisfies D_{ij}>=delta sum(x_{ij}) for every r,s in [0,1] iff it does at the four corners (0,0),(1,0),(0,1),(1,1). In particular delta=0 gives exact Gold transfer.

**Proof.** Every entry of X(r,s) is the bilinear interpolation of the four corner entries, with nonnegative weights

    (1-r)(1-s), r(1-s), (1-r)s, rs.

Both the deficits and the total are linear in the entries, so the claimed inequalities interpolate. Necessity follows by taking the corners, realized by empty and singleton contexts. Equivalently, normalized path deficits are convex combinations of their corner values, with weights additionally multiplied by the corner totals. ∎

This theorem uses one common strategy at all corners. Four different corner strategies would not establish it. It also uses the specific interface geometry: arbitrary induced-core fibres need not have form (2).

## 2. The core and the exact joint count

The core B has labels 0,...,8 and predecessor masks

    (0,0,1,0,8,24,3,5,189).

Equivalently its covers are

    0<2<7<8,
    3<4<5<8,
    0<6, 1<6.

The marked ports are a=0 and b=6. Use these two comparisons, in this order:

    first: 2 versus 4;
    second: 1 versus 2.                                     (5)

The natural increasing answer is the first row/column. The exact flag tables are

    X_0  = [[170,174],[236,84]],
    X_a  = [[67,129],[58,26]],
    X_b  = [[26,54],[44,36]],
    X_ab = [[10,39],[10,11]].                                (6)

Their totals are 664,280,160,70. Thus

    Z(r,s)=664+280r+160s+70rs.

### Short hand derivation, independent of poset enumeration

Temporarily delete vertices 1 and 6. The remaining extensions are the 20 shuffles of the three-element chains (0,2,7) and (3,4,5), followed by 8. Let t be the position of 0 and q the position of 2 in the six-letter shuffle. Then

    2<=q<=5, 1<=t<q,
    multiplicity at (t,q) = 6-q,
    2 precedes 4 iff q<=3.

Insert vertex 1 into a gap j in {0,...,7} of this seven-letter backbone. Then 1 precedes 2 iff j<q. The number of allowed slots for 6, which must follow both 0 and 1, is

    W(t,j)=8-t if j<t, and 8-j if j>=t.                      (7)

The flag a-first imposes t=1,j>=1. The flag b-last replaces W by 1. For each of the four cells, sum (6-q)W(t,j) under the two event tests; impose the flags or replace W as appropriate. These short integer sums give (6). For example the upper-left unflagged cell is

    4*14 + 3*20 + 3*18 = 170.

No large census or floating-point step is involved in these counts.

## 3. Exact context budget and sharp uniform margin

The four corner tables are

    (0,0): [[170,174],[236,84]], total 664;
    (0,1): [[196,228],[280,120]], total 824;
    (1,0): [[237,303],[294,110]], total 944;
    (1,1): [[273,396],[348,157]], total 1174.

Their four deficits, in row-major order, are respectively

    (150,146,108,260),
    (204,172,144,304),
    (167,101,246,430),
    (232,109,321,512).                                      (8)

Every entry is positive. The smallest normalized deficit over all corners is 109/1174. By the four-corner theorem, this proves (1) in every context.

The constant is sharp for strategy (5) over this context class: at two singleton collars, the first increasing answer leaves 669 extensions and the second decreasing answer leaves 396, so

    1174-669-396=109.

For direct use without a table, the four deficit polynomials are

    D_00=150+17r+54s+11rs,
    D_01=146-45r+26s-18rs,
    D_10=108+138r+36s+39rs,
    D_11=260+170r+44s+38rs.                                 (9)

The two negative coefficients in D_01 explain why checking only coefficientwise positivity would miss a valid certificate. Corner interpolation is the relevant positivity test.

The first comparison alone has the previously obtained balance bound 101/236; that did not automatically imply this two-step result. The new information is the joint table for the second comparison.

## 4. An unbounded seven-thin family

Take k>=1 copies B_1,...,B_k. Order every earlier block below every later block, except make the adjacent pair b_i and a_{i+1} incomparable for i=1,...,k-1. All other interblock relations remain.

Fix any block B_i. Its prefix is an allowed A-context with u=b_{i-1}, and its suffix is an allowed C-context with v=a_{i+1}; empty contexts cover the ends. Therefore the two comparisons in (5), translated into B_i, satisfy (1). This proves Gold for every k, with the same strict margin in every chosen block.

The structural parameters are explicit:

- Width exactly three. The core antichain {0,1,3} witnesses the lower bound. The core chain cover (0,2,7,8), (3,4,5), (1,6) extends blockwise, since no missing bridge joins consecutive members of one of these three chains.
- Maximum incomparability degree exactly seven. Core degrees are (4,7,5,5,5,5,6,5,2). An interior a-port gains one incomparable predecessor, and an interior b-port gains one incomparable successor; their new degrees are 5 and 7. Vertex 1 already has degree 7. No other cross-incomparabilities occur.
- Ordinal indecomposability. The core incomparability graph is connected, and each missing bridge joins two consecutive core incomparability graphs. Hence the whole incomparability graph is connected, precluding a nontrivial ordinal-sum cut.

This construction is not simple lexicographic module substitution: a block's a and b ports have different exterior relations from its other vertices. The core itself has only chain-valued proper modules (checked separately), so a nonchain proper module does not explain its certificate.

Peczarski's published six-thin theorem does not directly apply to a seven-thin source. Nevertheless no claim is made that this particular infinite family has never been covered by another known sufficient criterion. The proved new result is the explicit arbitrary-context, fixed-strategy bound (1). No majority-oriented pair with at most one slave exists at any chain length, so the elementary low-slave screen does not supply these certificates. The stronger companion `gold_unbounded_flag_context_20261019.md` now upgrades this exclusion to every chain length by a sixteen-row core certificate and an exact cross-bridge argument; the proof of (1) does not depend on that exclusion.

## 5. Validation and remaining proof obligation

`gold_four_collar.py` computes the hand sums, independently enumerates all 664 core extensions, and independently counts the actual four collar posets by ideal dynamic programming. It checks 81 nontrivial exterior-context/port combinations and every block in bridge chains of 1 through 12 copies (up to 108 vertices). Normal and `python3 -O` outputs agree. All comparisons and margins use exact integers or rational numbers.

Reproduce from the repository root; helper scripts are together in `05-knowledge/experiments`:

    python3 05-knowledge/experiments/gold_four_collar.py
    python3 -O 05-knowledge/experiments/gold_four_collar.py

The all-length proof is kernel (2), hand count (6), and corner inequality (8); the 108-vertex controls are checks of that proof's consequence rather than its logical foundation.

The remaining general challenge is to find an appropriate finite-boundary kernel and one common two-comparison strategy at every boundary extremum. These examples show that the right criterion is simultaneous positivity of the same path deficits. They do not show that every nonchain poset contains a core/interface system of this kind, and they give no bound on the necessary interface size in general.
