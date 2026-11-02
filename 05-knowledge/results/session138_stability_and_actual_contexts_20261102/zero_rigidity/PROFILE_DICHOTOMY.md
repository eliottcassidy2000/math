# Ordered-profile dichotomy and restricted strong zero rigidity

**Status: PROVED, all orders, with a finite 15-core certificate explicitly
listed and independently counted.** This does not settle arbitrary outside
cuts in the strong-quotient zero-rigidity problem.

Let `T=TT_k[C3,...,C3]`, delete v in block r, and let U be **any** completion
of that card by x. Let f_W(j) count ordered covers of all vertices of W by
j nonempty directed paths. Use the block indicator conventions of
`ORDINAL_ZERO_CLASSIFICATION.md`: `(a,b)` on the surviving pair u->w,
and t_i on intact triangles. Put `f_W=(f_W(1),...,f_W(3k))`.

## Theorem A: every completion is profile-comparable

Exactly one of the following occurs:

1. U is the canonical completion, and `f_U=f_T`.
2. `f_U<=f_T` entrywise, with `H(U)<H(T)`.
3. `f_U>=f_T` entrywise, and either `H(U)>H(T)`, or
   `H(U)=H(T)` and `f_U(2)>f_T(2)`.

In particular no noncanonical completion has the same full ordered profile.
The local H-zero cases are precisely those already classified in
`ORDINAL_ZERO_CLASSIFICATION.md`.

## Theorem B: cyclic strongification is zero-rigid for internal edits

Let A,B be any two nonempty tournaments, and let

    G = C3[T,A,B].

Delete v inside T. Suppose the new vertex x has its original orientation to
every vertex of A and B, but arbitrary orientation to the surviving vertices
of T. Then

    H(G-v+x)=H(G)  if and only if  x is the canonical completion of T-v.

Thus every nontrivial **internal** one-vertex replacement of the ordinal
triangle module changes H. The ambient quotient may be strong; this is a
restricted zero-rigidity theorem, not a statement about cuts also changing
arcs to A or B.

## 1. Colored odd-cycle polynomial and profile positivity

For integer m>=0, color all vertices of W with m available colors, and in
each color class choose one directed Hamiltonian path (empty classes count1).
The count is the polynomial

    P_W(m) = sum_j binom(m,j) f_W(j).                    (1)

Apply the inherited odd-cycle-family formula for H independently within
each color class. A family C of vertex-disjoint odd directed cycles, all
of length at least three, contributes

    2^|C| m^(|V(W)|-sum_(cycle in C)(|cycle|-1)).         (2)

Indeed each cycle must be monochromatic, and every unused vertex has its
own arbitrary color. Thus P_W is a sum of nonnegative monomials indexed by
actual odd-cycle families. This is derived directly from the Session137
OCF identity, rather than imported as a new external theorem.

A polynomial with nonnegative monomial coefficients has nonnegative
binomial-basis coefficients, since

    m^e = sum_j j! S(e,j) binom(m,j),                    (3)

where S(e,j) counts set partitions. Therefore inclusion of odd-cycle-family
collections gives entrywise ordered-profile domination.

Write `p(m)=m^3+2m=P_C3(m)`. An ordinal sum has multiplicative P, because
all directed cycles lie in individual summands. More generally, if a
vertex partition of W has parts W_1,...,W_s, then

    P_W - product_i P_(W_i)

has nonnegative monomial coefficients: retain only the cycle families
lying within the individual parts. No assumption on the cross arcs is
needed. Multiplication by p preserves profile inequalities, either by
counting colorings or by the positive ordinal-join law (6) in
`ORDINAL_ZERO_CLASSIFICATION.md`.

## 2. The low-profile states, completely explicit

Call the intact profile **monotone** if it is a sequence of zeros, at most
one mixed indicator1 or2, and then threes (omitting the broken block r).

The low states are exactly:

* `(a,b)=(0,1)` with the canonical external profile (all earlier t_i=0,
  all later t_i=3);
* `(a,b)=(0,0)`, with no selected vertex in an earlier block and with a
  monotone intact profile;
* the dual condition for `(a,b)=(1,1)`: no unselected vertex in a later
  block, and a monotone intact profile.

If no intact block is mixed, U is an ordinal sum of the k-1 intact
triangles and three transitive singleton vertices (possibly grouped into
TT2 or TT3). Hence

    P_U-p^k = -2m p^(k-1).                             (4)

If one intact block is mixed, its union with x is the unique strong
four-vertex tournament: it has two triangles and no longer odd cycle, so
its P is `m^4+4m^2`. The surviving root pair is a separate transitive
summand. Consequently

    P_U-p^k
      = p^(k-2)[m^2(m^4+4m^2)-p^2]
      = -4m^2 p^(k-2).                                (5)

The right sides have strictly negative value at m=1 and nonpositive
monomial coefficients. Equations (1),(3) prove alternative2 of Theorem A.
The mixed case requires k>=2; k=1 is covered by (4).

## 3. Every other state has a bounded positive core

If `(a,b)=(1,0)`, all original triangle cycles and their disjoint families
survive. A wrong incidence with an earlier block creates the new triangle
`x->z->w->x`; a wrong incidence with a later block creates
`x->u->z->x`. Thus either the cut is canonical or P_U-p^k has nonnegative
monomial coefficients and is strictly positive at1. This proves the claim
for the cyclic pair state.

For the other three pair states, use reversal to reduce to these cases:

* **Middle pair, wrong external incidence:** `(a,b)=(0,1)` and x fails to
  beat all vertices of some later intact block (or the reversed earlier
  case). The root pair, x, and that intact triangle form a six-vertex core.
* **Pair00, a selected earlier vertex:** root pair, x, and that earlier
  intact triangle form another six-vertex core.
* **Pair00, no selected earlier vertex, but not monotone:** there are two
  later intact blocks i<j with t_i>0 and t_j<3. They, the root pair, and x
  form a nine-vertex core.

These cases exhaust the non-low, non-cyclic states. In each core, delete
all intervening blocks: ordinal orientation between the retained blocks
is unchanged. Triangle rotation permits the canonical t_i representatives.
Thus only the following fifteen finite cores occur.

### Exact local profile differences

For the six-vertex cores, list the first four coordinates of
`f_core-f_(C3=>C3)`; the last two are zero. Pair state2 means `(a,b)=(0,1)`;
state0 means `(0,0)`.

| root position | pair state | intact indicator | difference coordinates1–4 |
|---:|---:|---:|---|
| earlier |2|0| (6,60,144,96) |
| earlier |2|1| (4,56,144,96) |
| earlier |2|2| (0,24,72,48) |
| later |0|1| (6,60,144,96) |
| later |0|2| (16,128,288,192) |
| later |0|3| (24,168,360,240) |

For the nine-vertex cores the root pair is first, its state is00, and the
following intact indicators are s in{1,2,3}, t in{0,1,2}. Their differences
from `f_(C3=>C3=>C3)` are entrywise at least

    (4,648,8088,35040,67920,60480,20160,0,0).            (6)

Equality in all coordinates of this lower envelope is attained at(s,t)=(1,2).
Every one of the nine full vectors and every literal tournament mask is
retained in `profile_core_certificate.out`; this is not an inferred bound.

For any selected core C using s original blocks, retain the actual cycles
within C and within each of the k-s other intact triangles. Section1 gives

    f_U >= f_(C => C3 => ... => C3)
         >= f_(TT_s[C3] => C3 => ... => C3) = f_T.       (7)

The second inequality is the positive ordinal-join profile law applied to
the displayed finite core differences. If the core has positive first
coordinate, (7) is strict in H. If it is the unique zero-first-coordinate
six-vertex type, its second coordinate increases by24 and every additional
ordinal C3 multiplies that first nonzero difference by12. Thus (7) either
is strict in H or is strict in coordinate2. This completes Theorem A.

## 4. Proof of restricted strong rigidity

For a fixed substitution context, the Hamiltonian-path count is a
nonnegative linear functional of the substituted module's f profile:
contract each maximal module run. For the outer C3 context, the word

    module, A, B

gives a positive coefficient of f_module(1), and

    module, A, B, module

gives a positive coefficient of f_module(2), since H(A),H(B)>0.
Alternative2 of Theorem A therefore strictly decreases the ambient H;
alternative3 strictly increases it. Equality holds only for the canonical
module completion. This proves Theorem B.

## 5. Independent leaf verification and boundaries

`profile_core_certificate.py` uses a direct subset DP for ordered covers:
a valid next adjacency multiplies the path-cover polynomial by1+t and an
invalid adjacency by t. It outputs all six six-vertex and nine nine-vertex
core differences. The six-vertex cases are also checked by literal720-word
enumeration inside that script.

`profile_literal_audit.cpp` is an independent literal-permutation verifier:
it counts descents of every labelled ordering and chooses separator cuts.
Its seventeen inputs in `profile_leafs.txt` are the two ordinal controls and
the fifteen cores. It checked all3,633,840 permutations, reproducing every
whole profile—not just the sign of a selected coefficient. The compilation
uses `-std=c++17 -O2 -Wall -Wextra -Werror`. Input/output receipts are retained.
No large-order finite extrapolation is used: all-order coverage follows from
the explicit case split and odd-cycle-family inclusion.

Outside cuts that alter x's orientations to A or B destroy the module
context and are NOT covered. Strong quotients with arbitrary structure and
arbitrary attachment changes remain OPEN. The result proves that the full
ordinal-module strongification route is closed, rather than merely showing
that its first small controls fail.

An additional whole-theorem control, `profile_dichotomy_audit.py`, checks all
228 block profiles at all roots for k=1,2,3. It verifies the exact low
polynomials (4),(5), every full-profile inequality, and then directly counts
each outer `C3[U,K1,K1]` by ordinary vertex-subset Hamiltonian DP. Every
context has only the canonical zero, matching the independently derived
formula `3f_U(1)+f_U(2)`. Ordinary and `-O` runs agree.
