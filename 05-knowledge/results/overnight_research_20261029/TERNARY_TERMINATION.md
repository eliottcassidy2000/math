# Pointwise termination of the cheaper ternary parser

Status: **PROVED elementary**, with exact finite controls. This proves
termination of a specified auxiliary parser on every positive odd
integer. It does not prove Collatz convergence or a decreasing return.

**Incoming integration, 2026-10-29.** While this alternative potential
was being checked, the concurrent packets
`session97_overnight_arithmetic_carriers_20261029/controllers/INTEGER_TERMINATION.md`
and `ternary_integer_returns_20261029/INTEGER_RETURNS.md` reached the same
termination theorem with stronger bounds. The former uses the sharp
branch factor 27/28; the latter proves J<=69K, D<=139K for
K=max(1,ceil(log2 n)), and classifies every expanding first return into
five exact source cylinders. Those results supersede a search for this
termination statement or its finite non-descent classification. The
potential here remains an independently derived alternative proof and
the general (d,r) criterion remains useful; no priority claim is made.

## T1 The retained source and queue

At a positive odd integer n take the actual accelerated Collatz step

    a = v_2(3n+1),        n' = (3n+1)/2^a.

Start with any positive integer queue Q. The cheaper ternary rule of
[MOMENTS](../affine_tree_recurrences_20261028/MOMENTS.md) updates

    Q' = Q+2 if 3 divides a, and Q'=Q-1 otherwise.

Stop at the first Q=0. The queue counts open ordered tree ports. The
rooted-tree case starts at Q=1; arbitrary Q represents a supplied forest.
No chronological letters are reordered or redrawn independently.

**Theorem.** This stopping time is finite for every positive odd source n
and every positive starting queue Q. Thus the exceptional set of
nonterminating odd 2-adic sources has empty intersection with the positive
integers, for this particular grammar.

## T2 A height and port potential

Set `V(n,Q)=n(8/5)^Q`. While n>=9, every legal step satisfies

    V(n',Q') <= (224/225) V(n,Q).                           (T1)

For a leaf, a>=1 and Q'=Q-1, so

    V'/V <= (5/8)(3+1/n)/2 <= 35/36 < 224/225.

For a branch, a>=3 and Q'=Q+2, so

    V'/V <= (64/25)(3+1/n)/8 <= 224/225.

If the queue never empties and the orbit never enters {1,3,5,7}, (T1)
would force V to zero, whereas V>=9*(8/5). This is impossible. Once
the orbit enters that four-element core, its future is explicit:

    1 -> 1,
    3 -> 5 -> 1,
    5 -> 1,
    7 -> 11 -> 17 -> 13 -> 5 -> 1.

The last path has valuations (1,1,2,3,4), queue increments
(-1,-1,-1,+2,-1), and net change -2. If Q reaches zero earlier we stop;
otherwise it reaches 1, whose valuation is always 2 and hence consumes
one queue unit per step. Every case therefore terminates. QED.

This proof uses positivity and an Archimedean height. It does not infer a
pointwise conclusion from a Haar-null set. The unbounded, history-dependent
queue also lies outside the finite-valuation obstruction in **THM-4507**,
`01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`.

An explicit deadline is available. Put V0=n(8/5)^Q0, rho=224/225.
If n>=9, let h be the least nonnegative integer with
`rho^h V0 < 72/5`. Within h steps the parser either stops or enters the
core. At first core entry the queue is at most
`floor(log_(8/5) V0)` because n>=1. Core passage then takes at most
`5+floor(log_(8/5) V0)` more steps. Both bounds can be computed by exact
rational comparisons, without floating logarithms. For an initial core
source, `Q0+5` is an immediate bound.

## T3 Why 27 versus 32 is the decisive boundary

For a grammar with d children at valuations divisible by r, try
`V=n c^Q`. At large n the worst leaf multiplier is `3/(2c)` and the
worst branch multiplier is `3 c^(d-1)/2^r`. A strict contracting choice
of c exists exactly when

    3/2 < c < (2^r/3)^(1/(d-1)),
    equivalently 3^d < 2^(r+d-1),                         (T2)

for d>1. This equivalence concerns this potential ansatz, not every
possible proof of termination. In the (3,3) case it is 27<32, and c=8/5
works above the exact core. The old critical (3,2) parser has 27>16;
this argument does not extend to it. For other grammars satisfying (T2),
one still owes the finite-core orbit certificate before concluding
pointwise termination.

The factor 27/32 also occurs in the preceding source-conditioned moment
bound in [REFINEMENT](../affine_tree_recurrences_20261028/REFINEMENT.md).
Here it arises from balancing two leaf costs against one branch cost;
the earlier derivation used a three-branches-per-bit tail budget. The
common exponential balance is real, but the pointwise proof additionally
needs the queue and the small-integer certificate.

## T4 Source preserving return controller

The controller keeps `(n_initial,n_current,Q,J,D,C)` with

    n_current = (3^J n_initial+C)/2^D.

After a letter a, update `C <- 3C+2^D`, `D <- D+a`, `J <- J+1`, and Q
as above. A complete return word has exact source guard

    n_initial = (2^D-C) 3^(-J) mod 2^(D+1).

This is the inherited faithful affine carrier of **THM-4604**,
`01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md`.
The return is certified decreasing precisely when

    (2^D-3^J) n_initial > C.                              (T3)

For concatenation uv, keep
`C_uv=3^(J_v) C_u+2^(D_u) C_v`, add the two clocks, and pull the second
guard back to the original source. Resetting Q for another return costs
a factor 8/5 in V; that reset is not free. Repeated parser termination
therefore supplies no proof that (T3) eventually holds.

**Hostile:** source 3 has a one-letter completed tree and returns to 5.
The exact census through 199,999 finds 51,710 nondecreasing returns among
100,000 positive odd sources. This is a finite statistic, not a density
claim. The same finite certificate checks arbitrary queues 1 through40
at the nine displayed test sources.

For `3n-1`, the large-height contraction is even easier, but the small
core changes: 5 -> 7 -> 5 is a genuine positive cycle. Both letters are
leaves for this parser. Parser termination can coexist with a nontrivial
arithmetic cycle. This is another decisive reason not to equate the
auxiliary termination statement with a Collatz conjecture.

## Recovery and the next obligation

The least-used useful predecessor is **THM-4654**,
`01-canon/theorems/THM-4654-safe-mahler-codes-contract-collatz-logarithmic-memory-and-complexity-approximation.md`:
it combines real contraction, a bounded rational window, and 2-adic
memory. The present proof keeps an unbounded port queue instead of its
safe-code constraint. The map is from an actual integer orbit to the
augmented `(n,Q)` state, preserving every step and valuation; projecting
away Q destroys the potential. No map between the Mahler problem and
this parser's nontermination set is asserted.

The remaining arithmetic problem is a source-preserving controller that
pays every restart and forces (T3), or a precise obstruction to such a
controller. That is stronger than the now-proved parser termination.
