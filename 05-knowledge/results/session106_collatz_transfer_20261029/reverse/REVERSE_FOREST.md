# A sharp width-three reverse forest and native-clock uniqueness

**PROVED, scoped to the five growing returns on positive multiples of 81 in shifted coordinates m=(n+1)/2.
Thus S consists of m=81t; the native odd source is n=162t−1, not a
multiple of81.**
The incoming macro suffix decoder has a stronger single-return consequence:
**at any fixed completed-return depth an endpoint has at most three legal
ancestors, and their native clock pairs are distinct.** Endpoint plus native
clocks therefore determines at most one legal source and word in this layer.
This is not a forward termination or reset-payment theorem.

## Provenance and independent checks

The five guarded maps and incoming adjusted-unit decoder are from
[`seven_return_native_recursions_20261029/COLLATZ.md` at
63ef3ca1172665f7100d06d405a3aa34c42efa57](https://github.com/eliottcassidy2000/math/blob/63ef3ca1172665f7100d06d405a3aa34c42efa57/05-knowledge/results/seven_return_native_recursions_20261029/COLLATZ.md).
Original-source marking and the native clocks are inherited from
[`atomic_tail_reset_20261029/ATOMIC.md` at
201fa96604d51171909318cc482220e00df54c38](https://github.com/eliottcassidy2000/math/blob/201fa96604d51171909318cc482220e00df54c38/05-knowledge/results/atomic_tail_reset_20261029/ATOMIC.md).
These identify the supplied incoming paths/commits; the local snapshots are
`../history/SEVEN_INCOMING.md` and `../history/ATOMIC_INCOMING.md`.

`check_reverse.py` imports no producer code, uses integer inverse equations,
and verifies each accepted edge by literal odd Collatz valuations. Its
`checks.json` records 10,000 endpoints in the 81-layer to reverse depth 12,
a sharp width-three witness, all five branches, the 27-layer collision,
the marked-source hostile example, and runs of up to 257 zero letters.
Both receipts include the checker's SHA-256. Reproduce from the workspace root:

```sh
python session106_collatz_transfer/reverse/check_reverse.py
python -O session106_collatz_transfer/reverse/check_reverse.py
```

The second command writes `checks_optimized.json`, separately from the normal
`checks.json`. All invariants use explicit exceptions, not removable asserts.
Finite checks support implementation correctness; the following proof
covers unbounded valuations and arbitrary finite depths.

## 1. Exact local predecessor rule

Write S=81 Z_{>0}, v=v_3(y), and
`chi(y)=2^v (y/3^v) mod81`. The incoming guarded maps have the form
`B_i(x)=3^e(a x+b)/2^d`:

|i|e|a|b|d|x mod 2^d|native J|
|--:|--:|--:|--:|--:|--:|--:|
|0|1|1|0|1|0|1|
|1|4|1|1|6|63|4|
|2|6|3|11|11|679|7|
|3|5|9|25|11|1135|7|
|4|4|27|59|11|1439|7|

For each i form

    x_i=(2^d y−3^e b)/(3^e a).

Retain it exactly when it is a positive integer in S satisfying its source
guard. These tests are necessary and sufficient; integrality alone is not
a replacement for the guard. In particular, the zero predecessor is
`2y/3`, and it exists in S exactly when v>=5.

For a nonzero predecessor, `v=e_i` and `chi(y)=c_i`, where

    (e_1,e_2,e_3,e_4)=(4,6,5,4),
    (c_1,c_2,c_3,c_4)=(61,13,70,34).

This follows because a source in S makes `a x+b` a 3-adic unit and reduces
it to b modulo81. Consequently at most one nonzero predecessor exists.
The test uses exact v, not a bounded residue approximation.

## 2. Comb theorem: a continuing path with short terminal side paths

At each reverse vertex choose its nonzero predecessor if one exists;
otherwise choose its zero predecessor if one exists. Call the resulting
path the **spine**. It is finite for a fixed positive endpoint since every
legal inverse strictly decreases the positive integer. This observation
alone says nothing about infinite forward paths.

The only branching occurs when a nonzero predecessor and a zero predecessor
both exist. It has exactly two possible types:

* Nonzero type 3: v=5, chi=70. Taking the zero alternative leaves valuation4
  and chi70. It has no further zero predecessor in S and no nonzero
  predecessor, since the valuation4 nonzero types have chi61 or chi34.
* Nonzero type 2: v=6, chi=13. Taking zero leaves valuation5 and chi13, where
  type3 is excluded. One more zero leaves valuation4 and chi13, where both
  valuation4 types are excluded. The side path therefore has length2.

Types1 and4 have v=4, so they cannot coexist with a zero predecessor in S.
The zero inverse preserves chi. Thus every off-spine component is exactly
a zero-only terminal side path of length1 or2, never a recursively
branching subtree. Arbitrarily long forced zero segments may occur on the
spine and have not been discarded.

**Width bound.** At reverse depth k there is at most one spine vertex, one
side vertex from a split at depth k−1, and one from a split at depth k−2.
Hence the width is at most3, independently of k and endpoint height.

**Sharpness.** The endpoint 97551864 has precisely these depth2 ancestors:

|source|forward word|native J,D|
|--:|:--|:--|
|43356384|00|2,2|
|60901146|02|8,12|
|85545639|22|14,22|

The last path is `85545639 -> 91351719 -> 97551864`, both type2.
Every displayed source belongs to S and every edge passes its actual guard.

## 3. Native clocks remove the remaining ambiguity

Let the **clock defect** be H=D−J. A zero edge contributes0; type1
contributes2; types2,3,4 contribute4. At a fixed reverse depth, different
ancestors have distinct H: each side path retains the spine's defect at
its splitting vertex, whereas continuing along the spine immediately
adds a positive defect. If two different side paths reach the same depth,
the later split comes after that strictly positive addition. The spine
vertex, when present, likewise has larger defect than every earlier side
path reaching its depth.

Every branch also obeys `5J_i−3D_i=2`. Therefore a k-return word satisfies

    k=(5J−3D)/2.

**Clock uniqueness theorem.** For an endpoint y in S and specified native
clocks J,D, there is at most one legal growing-return word, and at most one
source in S, having those endpoint and clocks. Indeed the clocks fix k,
and fixed-depth ancestors have distinct H. There may be no accepted path.

Equivalently, with return count k and defect H, the native clocks are
`J=k+3H/2`, `D=k+5H/2`. H must be a nonnegative even integer. These are
necessary arithmetic checks, not sufficient legal source guards.

This gives an exact finite-width reverse controller: enumerate inverse
edges with the local rule while retaining remaining clocks, and stop
exactly at k. There are at most three candidates at each layer, and the
clock filter leaves at most one final candidate. At most three vertices and
fifteen candidate branch equations are examined per layer, so a breadth-first
implementation uses O(k) inverse-node/edge checks and constant-width active
storage (apart from counters and word reconstruction). This is an arithmetic
operation count, not a bit-complexity bound: integer sizes, valuation
extraction, and division costs remain relevant. Test that candidate
against the supplied marked source. An alternative greedy controller can
use the fact that a zero side branch cannot pay any positive remaining
clock defect. Either version retains the unbounded valuation counter.

### Complete predecessor enumeration has logarithmic arithmetic depth

Every forward map has nonnegative additive constant and slope at least
rho=2187/2048>1 (the minimum is shared by types2,3,4). A k-edge predecessor
chain in S ending at y therefore satisfies

    y >= 81 rho^k, equivalently y*2048^k >= 81*2187^k.

Let L be the largest nonnegative integer satisfying this integer inequality.
There are no reverse vertices beyond depth L. By the width bound the entire
reverse tree has at most 1+3L vertices, and hence O(log(y/81)+1) candidate
inverse-edge checks suffice to enumerate all ancestors, not just one
requested clock layer. Every such word has native clocks J<=7L and D<=11L.
The integer inequality is the exact cutoff; floating-point logarithms are
unnecessary. This is still only an arithmetic-node count, not a claim of
logarithmic bit complexity. It bounds histories of a supplied finite
endpoint, not the size of a future endpoint of a supplied source.

## 4. Marked boundary, hostile tests, and actual Collatz relevance

The source mark cannot be replaced by the spine's beginning. In

    6399 --1--> 8100 --0--> 12150,

the clock pair (1,1) recovers 8100, whereas (5,7) recovers 6399. A query
marked at 8100 with clocks (5,7) is rejected, not silently re-marked.
Thus the new theorem makes the incoming overshoot warning an enforceable
clocked boundary condition, including when the mark is inside a long run
of zeros.

The theorem genuinely needs the stable layer: on 27Z,
`B_1(4671)=B_4(5535)=5913`. Both inverses are legal, so the local comb proof
fails there. Distinct chi constants require modulus81; their collisions
at smaller precision are not merely numerical noise. Our theorem does
not assert the same width or clock uniqueness outside S.

This yields a concrete leverage point for tournament-style recursive
certificates: **a clock-preserving smaller-source proposal can be checked
without exponential reverse branching once both sources lie in S.** In
fact, two different sources in S cannot have the same endpoint and the
same native clock pair through this five-return language. A proposed
collision switch with those properties must leave the language/layer or
change clocks. This is a useful obstruction as well as a decoder.

It still gives no lower bound on the least forward survivor s_k. The
spine may keep acquiring larger endpoints in forward time; clocks and
sources cannot be reset when a local return fails. A descent theorem
requires an additional source-preserving rank or an admissible clock-
changing switch whose growing-return count and final source debt are
proved. Unbounded survivor escape, uniform seven-return growth, and
original-source reset payment remain open.
