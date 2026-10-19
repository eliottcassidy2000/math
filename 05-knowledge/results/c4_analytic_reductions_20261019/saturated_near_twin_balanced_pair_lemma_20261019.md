# Saturated common neighbors give exact balanced pairs

Date: 2026-10-19. Status: proved structural reductions, with a separate finite audit. No completeness claim for mixed (1,2), c=4 follows from this note.

Write D(v) and U(v) for the strict downset and upset of v in a finite poset. A pair is balanced when its ordering probability in a uniformly chosen linear extension lies in [1/3,2/3].

## 1. Equal lower sets and one extra upper element

**Lemma.** Suppose x and z are incomparable and

    D(x)=D(z),       U(x)=U(z) union {a},       a not in U(z).

Then

    1/2 <= p(x<z) <= 2/3.

**Proof.** Necessarily x<a and z is incomparable with a. The alternative a<z would imply x<z. Partition the extensions into three classes according to the relative order of x,z,a:

    E0: z<x<a,       E1: x<z<a,       E2: x<a<z.

There are no other possibilities. Swapping the labels x and z is a bijection E0 -> E1. Both have the same predecessors, and every successor of z is also a successor of x. In E0 the extra successor a follows x, and in E1 it follows z; consequently swapping preserves every required relation. In particular, |E0|=|E1|.

There is also an injection E2 -> E0. Replace the labels in their three occupied positions by z,x,a, respectively. Thus z moves to the old position of x, x to the old position of a, and a to the old position of z. All other positions are unchanged. This fixed permutation of three labels is injective. To verify that it gives a linear extension:

- Every predecessor of z is a predecessor of x, so it was before the old position of x.
- Every predecessor of x was also before the old position of x, and hence before its new, later position.
- Every successor of z is a successor of x and was after the old position of z, so it remains after both new positions.
- Every successor of x other than a is a successor of z and was after the old position of z. The new a follows the new x.
- Every predecessor of a other than x was before the old position of a and hence before its new position; x now occupies the old position of a. Neither z nor a is a predecessor of the other.
- Every successor of a is also a successor of x and is not a itself, hence belongs to U(z) and was after the old position of z.

Therefore |E2|<=|E0|. Since extensions with x<z are E1 union E2 and those with z<x are E0,

    |E0| <= |E1|+|E2| <= 2|E0|.

Division by the total proves the bounds. The upper endpoint is attained by the three-element poset having only x<a. The pair is balanced in either event. The order-dual statement follows by reversing every relation. ∎

**Twin case.** If instead D(x)=D(z) and U(x)=U(z), swapping the two labels is an involution exchanging the two ordering events, so p(x<z)=1/2.

## 2. Application to a seven-thin mixed (1,2), c=4 core

Use root labels x,y, lower slave b, upper slaves a1,a2, and four common incomparable neighbors z1,...,z4. Thus b<y, x<a1,a2, and every zi is incomparable with both roots. The core Q is I(x) union I(y), including the roots. The exterior is A below both roots and B above both roots.

Suppose the ambient poset has incomparability degree at most seven. If a common vertex z has exactly one comparable core neighbor, its seven other core neighbors exhaust its incomparability allowance. The saturation lemma therefore forces every A element below z and every B element above z.

1. If the sole core neighbor is b, its direction must be b<z: z<b<y would contradict z incomparable with y. Hence

       D(z)=A union {b}=D(y),       U(z)=B=U(y).

   The pair z,y consists of true twins, with probability exactly 1/2.

2. If the sole core neighbor is an upper slave aj, its direction must be z<aj: x<aj<z would contradict x incomparable with z. With ak the other upper slave,

       D(z)=A=D(x),       U(z)=B union {aj},
       U(x)=B union {aj,ak}=U(z) union {ak}.

   The lemma gives 1/2<=p(x<z)<=2/3.

These are ambient identities, not inequalities only on a relaxed generator cone. They rule out the corresponding core in a poset without any balanced pair, without requiring that x<y be a majority edge or have maximum index span. No assertion is made here about a saturated common vertex whose sole neighbor is another common vertex.

## 3. Separate exact finite census

The naturally labeled c=4 morphology family contains 96,428 cores. Exactly 79,121 satisfy the seven-thin degree bound at every core vertex. Among the latter:

- 13,885 have a saturated common vertex whose sole core neighbor is b;
- 21,893 have one whose sole core neighbor is an upper slave;
- 4,210 meet both tests;
- 31,568 meet at least one test;
- 47,553 remain after both exact reductions.

These are finite census counts, not the proof of either lemma. Remaining cores are not claimed to be counterexamples or to admit any particular certificate. The separate full c=4 cone scan retains its own scope and completion status.

## Reproduction

Run `python -O verify_saturated_near_twin_lemma_20261019.py` in this directory. The self-contained standard-library verifier checks the extension maps on every naturally labeled poset through six vertices and independently regenerates the c=4 census. It writes `saturated_near_twin_lemma_controls_20261019.json`.
