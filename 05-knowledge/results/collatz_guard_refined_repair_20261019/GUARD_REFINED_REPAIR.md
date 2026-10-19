# Source-aware guard refinement beyond the 120-bit obstruction

Date: 19 October 2026.

## Main result

The 25-episode first-crossing theorem holds with a repair allowance of 183
Terras steps. Its new ingredient is a small symbolic partition of the original
source prefix into two longer lawful cylinders. This supplies genuine new
parity information and repairs the offset 12852 that the earlier 120-bit guard
could not handle.

A separate computer-assisted corollary extends the episode limit to 1000,
with allowance 273. It uses the same two cylinders and 15 small complete
refinement trees. The important mechanism is the all-lift guard refinement,
not the numerical episode limit.

## 1. Definitions and inherited symbolic facts

Use T(x)=x/2 for even x and T(x)=(3x+1)/2 for odd x. Read parity words in
chronological order. Use exactly the guarded banks

    K = 110110110110110101011011101010,
    F = 10,       H = 110101000,       V = K^4 F.

K and V are expansion blocks; F and H are resets. A V is one 122-step
expansion macroblock, including its internal F. A source is fully guarded
when all indicated parities are actual parities of its orbit.

Take a finite bank word beginning in K or V. Stop at the first endpoint of an
explicit F or H macroblock where the cumulative affine multiplier is less
than one. Call that Terras time N. Assume it exists, and let m count maximal
nonempty expansion runs before N. Every earlier explicit reset boundary has
multiplier at least one. The exact rewrite K^4F=V is permitted before counting
m; it preserves the first-crossing time and can reduce m.

The preceding
[interleaving proof](../collatz_interleaved_height_repair_20261019/INTERLEAVED_HEIGHT_REPAIR.md)
establishes, for every m and unrestricted individual word lengths:

- The full branch is Sx+e, with S<1 and e<mE, where

      E = 83986881009841395726881741282576539890
          /157489106280555890461927692448373459.

- If d=T^N(x)-x>=0, then d is an integer with d<mE.
- A zero displacement gives a lower cycle with least period N and minimum x,
  and a separately lawful adjacent companion cycle containing (x+1)/4,
  strictly between 2 and x. The cycles are distinct and nontrivial.

For clarity, the carry bound uses the maximal expansion center C=E-1.
If episode i has multiplier si and carry ei, then ei<C si+1. With cumulative
multipliers Si>=1 before the final S=Sm<1, summing ei S/Si gives e<m(C+1).
Every nonempty prefix of K and V expands; every proper F/H prefix at input
>1 exceeds its reset endpoint. Those facts prove the zero-displacement cycle
minimum, including times inside V. The companion identities are exact guarded
integer identities, not a density statement.

**Theorem.** Under these hypotheses, if m<=25, exactly one of the following
occurs:

1. T^t(x)<x for some 1<=t<=N+183.
2. T^N(x)=x, the lower cycle has least period N and minimum x, and a different
   nontrivial positive cycle contains (x+1)/4<x.

A least positive nonconverging integer cannot start such a word. This is a
statement about the original source of this word. It supplies no arbitrary
pointwise entry theorem and no claim about descent below a source that
preceded a later, higher bank entry.

## 2. A symbolic partition with only thirteen short exclusions

Any such first-crossing word either already descends at its endpoint, or its
literal source parity prefix is one of

    A = K^5       (150 bits),
    B = V K^3    (212 bits).

Here is a complete proof, including the short branches. Beginning with K,
read consecutive K blocks. If a V appears after one, two, three, or four K
blocks, its own initial K^4 makes the source begin K^5. Five consecutive K
blocks do the same. Before that, the possible resets are

    K^ell F or K^ell H,  ell=1,2,3;
    K^4 H;
    K^4 F = V.

The first seven words are already contracting and descend on every positive
lawful lift. The final word is expanding and is exactly V as a literal parity
word, whether it was originally parsed as a V or as four K blocks and a reset.
Thus the remaining branch begins V.

After V, either a new V appears before three further K blocks, or three K
blocks appear. In both cases the source begins V K^3, because a V starts K^4.
Before that, the only resets are

    V K^i F or V K^i H,  i=0,1,2.

All six words contract and descend on every positive lawful lift. Together
these are thirteen short exclusions. A first crossing cannot have occurred
earlier on the surviving path: every earlier explicit reset, if any, was the
expanding K^4 F. A complete word cannot end at a surviving expansion block,
because N is a reset boundary.

The all-lift check for each short word is particularly small. Reconstruct its
branch (a x+b)/M and unique lawful residue r modulo M. Verify

    a<M,      r>0,      (M-a)r>b.

Then the endpoint is below x for every x=r+uM, u>=0. The exact checker performs
these three tests on all thirteen words, and lists their residues, lengths,
and fixed points in the receipt. No small positive source in these cylinders
is skipped.

This partition is independent of the episode count and of the eventual
endpoint displacement.

## 3. The two-cylinder all-lift offset lemma

The two source cylinders are

    rA = 380358739127104371168868555681161388769345531 mod 2^150,
    rB = 2467150623986679092662520573921733544503296291794223140097949691
         mod 2^212.

For every d=1,...,13332, the actual prefix of rA+d has a contracting affine
multiplier within 135 steps, and that of rB+d within 183 steps. In each case,
if the first contracting branch is (a z+b)/2^j, the exact check gives

    (a d+b)/(2^j-a)
      <= 810692053755261919/1738366812781 < 466400.

The maximum threshold is the same as in the preceding 24-episode result.
Both source residues are much larger than 466400. The prefix length j is at
most the number of known source bits. Consequently every source lift
x=rA+u2^150, or x=rB+u2^212 respectively, satisfies the same continuation guard
at x+d and

    T^j(x+d) = (a(x+d)+b)/2^j < x.

This is a finite, exact, computer-assisted offset lemma. The verifier loops
over the full offset interval, reconstructs the affine branch, authenticates
the congruence and modulus containment, and checks the displayed inequality.
The checker also runs the larger interval of Section 5, groups identical
first-contracting cylinders, and
checks the worst offset in every class as an additional control. All offset
counts refer to displacement parameters; no conclusion is inferred from a
bounded set of Collatz starting integers.

Since

    13332 < 25E < 13333,

a failed endpoint descent with m<=25 has 0<=d<=13332. If d>0, apply the symbolic
partition and the offset lemma at time N. If d=0, use the cycle alternative
from Section 1. This proves the main theorem and its least-counterexample
consequence. The alternatives are exclusive because the lower cycle in the
second case has minimum x forever.

## 4. General finite guard-refinement criterion

The same idea can be formulated as a reusable certificate criterion for a
fixed positive displacement d.

Construct a finite forest of macroblock prefixes, with roots K and V. Every
nonterminal node has all four children obtained by appending K,V,F,H. Every
nonterminal node has cumulative multiplier at least one. Allow two types of
terminal node:

**Source descent.** Its multiplier is below one, and its own lawful residue r
satisfies the exact all-lift endpoint descent inequality (M-a)r>b. Because its
parent did not cross, this node is the first contracting reset boundary on
that path.

**Offset repair.** For a source prefix w of n bits and lawful residue r, some
prefix of r+d has length j<=n and branch (a z+b)/2^j with a<2^j, satisfying

    r>0,      (2^j-a)r>a d+b.

The branch guard must be checked modulo 2^j. It then applies to x+d for every
source lift x=r+u2^n and gives T^j(x+d)<x.

A finite complete forest with these terminal tests proves the following:
every fully guarded first-crossing word either already has endpoint descent,
or, if its endpoint displacement is d, admits the certified continuation
below its original source. The worst j over its offset leaves is a valid
uniform repair allowance.

To justify coverage, follow the given bank parse from its actual root. It
cannot end at a nonterminal node, since its endpoint must be a first crossing.
It cannot pass a source-descent leaf before its declared first crossing.
It therefore either ends at such a leaf with actual descent, or reaches an
offset-repair leaf whose source prefix occurs no later than N. Completeness of the
four children and finiteness leave no unexamined path. This uses the full
source prefix at a leaf, not guessed parity bits beyond an inherited guard.

## 5. Separate exact computational corollary: 1000 episodes

The same theorem holds for m<=1000 with N+273 in place of N+183.
The exact carry budget is

    533286 < 1000E < 533287.

On the two cylinders A,B, all offsets 1,...,533286 are repaired within 183
steps except for the following fifteen offsets, where at least one cylinder
needs refinement:

    14196, 15364, 18436, 49324, 128308,
    215276, 266092, 288932, 295660, 310916,
    390388, 394148, 480032, 522404, 529652.

For each exceptional offset the verifier constructs a complete forest of the
kind in Section 4. All terminal inequalities pass; no unresolved leaf is
accepted. The fifteen forests together have 838 nodes. Their largest offset
repair time is 273, and the longest recorded source guard at an offset leaf
has 394 bits. The full forests are stored in the receipt, including every branching node.
Offset-repair leaves record their source residue and exact threshold.
Source-descent leaves record their macroblock word and length; the verifier
reconstructs their residues and checks their exact all-lift inequalities.
This is a compact certificate rather than an omitted search.

The base two-cylinder coverage is authenticated over the entire declared
interval. There are 33,990 first-contracting offset classes in the combined
base calculation. Their worst exact threshold is

    19411997879661734889626095/420491770248316829,

less than either base source residue. The checker's separate per-offset test
and its congruence-class endpoint/count checks cover all offsets, including
the fifteen exceptions later handled by trees.

This corollary uses the same symbolic carry and zero-displacement argument.
It changes no claim about pointwise entry or unbounded episode counts. Its
finite arithmetic is completely reproducible by the single companion script,
without network access or importing an earlier verifier.

## 6. What an unbounded refinement theorem would need

A sufficient next theorem would say: for every positive integer displacement
d, a finite complete bank-prefix forest satisfying Section 4 exists. Combined
with e<mE for each finite m, this would remove the episode-count restriction
for the first-crossing bank language. It would still leave the separate
pointwise bank-entry problem.

Finite experiments do not prove this universal termination statement. There
is a concrete rational stopping condition already present on its simplest
infinite branch.

Write

    c = C0/D = 4912087879/88519643.

The infinite source itinerary K K K ... represents the odd-denominator
2-adic fixed point X=-c of the K branch. Its finite lawful residues rL modulo
2^(30L) are nondecreasing and unbounded: they converge 2-adically to the
negative rational -c, which is not a nonnegative integer. For each fixed d,
the parity prefixes of rL+d stabilize to those of the rational

    X+d = d-c.

Here T on an odd-denominator rational uses its numerator's parity after
putting it over an odd denominator; the same two affine formulas apply.

Along this infinite K path, no source-descent leaf is possible, because
all K-prefix multipliers exceed one. A finite offset-repair leaf exists on
that path if and only if d-c has a finite coefficient stopping time, meaning
that one of its parity prefixes has multiplier 3^h/2^j<1.

The forward implication follows from the matching guarded prefix. For the
reverse implication, fix such a contracting prefix. Its affine threshold
(a d+b)/(2^j-a) is finite. Eventually 30L>=j and rL exceeds that threshold,
so K^L supplies the required offset-repair certificate.

Thus universal finite-tree termination would in particular require finite
coefficient stopping for all positive rationals d-c with integers d>=56.
This is a precise necessary subproblem for this certificate method, not an
assumed consequence of a density-one theorem. It is not claimed to be
necessary for every possible proof about actual positive-integer bank words:
an infinite 2-adic branch can be an obstruction to finite cylinder covering
without itself being a positive integer source.

The distinction between coefficient contraction and actual source descent
is also essential. The latter needs the carry inequality. For background on
this distinction for integer sources, see Rozier–Terracol,
[Definition 1.2 and the discussion of the CST conjecture](https://arxiv.org/html/2502.00948v5).
No general CST assertion is used here.

## 7. Reproduction and boundaries

Run, from this directory:

    python verify_guard_refined_repair.py
    python -O verify_guard_refined_repair.py

Both outputs must match verification-receipt.json. Checks remain active with
Python optimization. The receipt distinguishes the clean 25-episode theorem
from the separate 1000-episode computational corollary and records complete
refinement trees. The written symbolic partition, carry argument, and finite
forest criterion supply the unbounded-word-length quantifiers.

The earlier 24-episode theorem remains intact. Its observation that offset
12852 is unresolved by the K^4 cylinder alone was correct; this note resolves
it using additional original-source guards. No Haar-almost-everywhere,
density-one, residue-equidistribution, or finite-word-frequency statement is
promoted to a pointwise integer entry theorem.
