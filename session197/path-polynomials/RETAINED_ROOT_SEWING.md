# A positive four-state retained-root sewing grammar and its exact finite inverse

Session197. **PROVED all-order physical recurrence and target-bounded inverse; FINITE-EXACT bank spectrum.** The retained mark stays in the old left graph. This is different from moving the mark into the right gadget. With every right block of order at most5, the exact arbitrary-depth spectrum through5000 contains every odd number except7 and21. No coverage beyond the stated bound, universal digit return, or full tournament spectrum is claimed.

## 1. Source attribution and the physical operation

The scalar single-edge join law is inherited from the pinned S195 rooted-interface proof:
https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/rooted-interfaces/PROOF.md . That source proves a three-coordinate grammar when the mark moves into the new right gadget; it also gives stationary retention for a singleton gadget. The pinned S196 star synthesis is a separate physical construction:
https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/spectral-moves/PROOF.md . Neither named result is being re-labelled as a new scalar law here.

For disjoint nonempty rooted tournaments (A,a),(B,b), form A->B and reverse exactly the edge a->b to b->a. Retain the mark a. Write S(A,a;B,b) for this actual graph. All old vertices and arcs remain, but A is generally not a module in the output. A word starts at a singleton and repeatedly applies this operation with a right block from a specified finite bank. The root never moves.

The additions here are the arbitrary-right-block stationary four-state closure, a positive coordinate system, its minimal linear dimension, a finite exact inverse with an explicit word/order bound, and the bounded bank certificate. The proof includes the scalar decomposition to fix all endpoint conventions.

## 2. Four physical coordinates and the exact transfer

For a rooted graph (A,a), let h=H(A), e be its number of Hamilton paths ending at a, and d=H(A-a), with H(empty)=1. Let R count ordered two-bin path covers (P,Q), where P may be empty and Q is a nonempty path starting at a. Define

    x=R-d, y=d-e, state(A,a)=(h,x,y,d).                  (1)

These are nonnegative integers, and h,d>=1. Indeed the root-singleton Q and any Hamilton path on A-a give d covers, so R>=d. Deleting the last a from a Hamilton word ending there injects into words on A-a, so e<=d. Also d<=h by deterministic tournament path insertion.

For a fixed right block (B,b), abbreviate

    H=H(B), E=end_b(B), S=start_b(B),
    Q=R_b(B), L=L_b(B), C=B_B(2),

where L counts ordered covers with the first nonempty path ending at b and the second path optional, and C counts all two labelled path bins, both allowed empty. Then the exact transfer is

    h' = H h + L x + S y + (L-S)d,
    x' = C x + (C-Q-H)d,
    y' = (H-E)d,
    d' = H d.                                          (2)

Every coefficient is a nonnegative integer. This is a physical four-coordinate closure for any sequence of actual right blocks, not a claim that every nonnegative state or matrix is realizable.

### Literal fragment proof

In the original coordinates (h,e,R,d), the four formulas are

    h'=Hh-Se+LR, e'=Ed, R'=CR-Qd, d'=Hd.                (3)

A Hamilton word avoiding b->a is one A-path followed by one B-path. The forbidden boundary a->b removes exactly eS words. A word using b->a has the unique run order

    A_before, B_before, A_after, B_after,

where A_after starts at a and B_before ends at b; the first and last runs may be empty. Choosing an R-cover of A and an L-cover of B gives these words bijectively. This proves the inherited scalar row.

A Hamilton word ending at a must use b->a as its final edge. Everything in A-a precedes B, whose path ends at b. This gives dE, including A singleton. Deleting a removes the only backward cross edge and leaves (A-a)->B, giving dH.

For a two-bin cover whose designated path starts at a, the edge b->a cannot be used in either bin: a is first in its own bin and absent from the other. Each bin therefore traverses its A part before its B part. An R-cover of A and an arbitrary labelled two-bin cover of B initially give RC choices. A forbidden boundary a->b occurs precisely when the designated A path is the singleton a and its B part starts at b. The other A path is a Hamilton word on A-a, and the B parts form a Q-cover. There are exactly dQ invalid covers. This proves the R row. Substituting R=x+d and e=d-y yields (2).

### Positivity and strict activation

Splitting each Hamilton word of B immediately before b injects into its Q-covers; splitting immediately after b injects into its L-covers. Hence Q,L>=H, while S,E<=H.

The two labelled bins treat the bin containing b symmetrically. Therefore C is twice the number of covers with b in one specified bin, and Q is a subset of those covers. Thus C>=2Q>=2H. It follows that C-Q-H>=0 and L-S>=0.

For |B|>=2, the coefficient C-Q-H is strictly positive. Cutting every Hamilton word at any of its |B|+1 gaps and putting its prefix/suffix in the fixed bin order gives an injection, so C>2H. Together with Q<=C/2, this implies C-Q-H>0. For B singleton, (H,E,S,Q,L,C)=(1,1,1,1,1,2), and this coefficient is0. Thus precisely the singleton letters fail to activate x from x=0,d>0.

The familiar retained-root singleton transfer is recovered exactly:

    (h,x,y,d)->(h+x+y,2x,0,d).

## 3. Minimal rational linear dimension

Four coordinates cannot be replaced by a three-dimensional rational linear realization for this fixed retained-root grammar, even when the right bank contains only rooted singletons and rooted transitive edges.

Let F be the edge gadget marked at its source and K the edge gadget marked at its sink. Starting from the singleton, take the four physical inputs given by words empty,F,K,FF. Their states are

    (1,0,0,1), (1,1,1,1), (3,2,0,1), (3,5,1,1).

Use the four physical continuations empty, one singleton, two singletons, and K. Their actual Hamilton-count matrix, with contexts as rows, is

    [1 1 3 3
     1 3 5 9
     1 5 9 19
     3 5 9 15],    determinant16.                       (4)

Every rank-three rational linear state/continuation representation would factor this matrix through a three-dimensional space, impossible. The initial encoding need not be linear for this rank argument. Equation(2) supplies dimension four, so the dimension is exact under these assumptions. This does not concern a nonlinear or unrestricted representation.

## 4. An exact finite inverse, without a guessed depth

For any physical input and any bank letter, (2) gives

    h'>=h+x, x'>=2x, d'>=d.                             (5)

Fix a target bound M. Start at (1,0,0,1). For each retained state and each bank letter, compute (2). If h'<=M, record the count and a physical predecessor recipe. Expand the new state only if h'+x'<=M, and only once for each exact four-tuple. A terminal state with h'<=M is still recorded even if h'+x'>M.

This algorithm terminates. Every expanded state has 1<=d<=h<=M, 0<=y<=d and 0<=x<=M-h; hence only finitely many integer states are possible. Identical states have identical scalar counts after every future bank word by the full returning recurrence, so merging them is exact.

It is complete for arbitrary word depth. Counts never decrease. If a word has final count at most M, each state before its last step satisfies h+x<=M by (5), so induction visits it. The final value is recorded even when its x is too large for another step. Thus the finite search decides membership in this bank language and provides an actual adjacency-matrix construction whenever successful.

### Exact initial-padding transport and graph-order bound

At the singleton seed, a singleton letter leaves the entire four-state equal to (1,0,0,1). Any initial string of singleton letters can therefore be removed without changing the state seen by any later right-block continuation. This is a proved transport lemma, not an assumption that internal inactive padding can be deleted. No internal letter is stripped.

If the word has any nonsingleton letter, take it as the first letter after this normalization. Strict activation gives x_1>=1. For a normalized word of m>=1 letters, (5) implies x_j>=2^(j-1) and

    h_m>=1+sum_(j=1)^(m-1)2^(j-1)=2^(m-1).

Hence a target N has a normalized witness, if any, with at most1+floor(log2 N) letters. If the bank's largest right block has k vertices, a canonical witness has at most

    1+k(1+floor(log2 N))                                (6)

vertices. Arbitrarily padded witnesses may be larger; the bound is on a count-equivalent canonical recipe. The all-singleton language only gives count1, represented by the singleton itself.

## 5. Complete small bank inventory and bounded spectrum

The verifier enumerates all1099 labelled tournaments through order5 and all5405 roots. It computes each letter's six constants by endpoint subset DP and independently by literal vertex permutations with every possible path cut. Deduplicating only equal six-tuples gives65 distinct retained-root transfer matrices, each with a literal physical representative. Restricting to order at most4 gives19 signatures. These counts happen to agree with the moving-root source's bank sizes, but the transfer matrices and retained-root graph words are different.

The exact arbitrary-depth inverse gives:

- Order-at-most4 right blocks, M=2000:20389 distinct reachable four-states, maximum breadth-first representative depth10,997 attained odd counts. The missing odds are exactly7,21,173.
- Order-at-most5 right blocks, M=5000:152044 distinct reachable four-states, maximum breadth-first representative depth12,2498 attained odd counts. The missing odds are exactly7,21.

Here reachable states include recorded terminal states that are not expanded. Counts and depths are reproducible implementation receipts, not minimal witness-order claims. Because the inverse is complete, the omissions within each stated bank and target range hold at every word depth. In particular173 is excluded by the entire order-at-most4 bank but is realized after allowing an order5 block.

One actual H173 witness has outgoing rows

    [188,225,194,198,206,220,1,64], root0.

It has8 vertices and is independently counted by the source's endpoint DP. The receipt also supplies complete words, graph rows, four-states and direct counts for67,163,179,181,247,511,539,541,821 and4999. These are independent constructions of those values. They do not assert that every existing H179 parent has a541 child or that the same return can be repeated.

The observed exact interval through4999 is not an induction theorem for larger targets. The unrestricted spectrum and even the all-order scalar coverage of this fixed bank remain open here.

## 6. Reproduction and verification scope

Run `python3 verify_retained_root_sewing.py` and its `-O` form. The script uses exception gates rather than removable assertions. It regenerates the complete labelled bank, checks both counting routes, verifies the full four-state against actual graphs for every rooted pair through order3 and further larger controls, checks (4) literally, runs the target-bounded inverses and constructs every reported witness graph.

The bank inventory and target recipes are written as reproducible JSON artifacts. An optional larger bound changes only the finite census and cannot certify unbounded coverage. The sole imported helper is the frozen `verify.py`, SHA25610c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb. The scalar source, its moving-root scope and its singleton stationary row are attributed in Section1. The all-order proofs above are independent of a finite numerical extrapolation.
