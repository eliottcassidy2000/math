# Marked periodicity and the factorial leakage barrier

## Different meanings of a return

The October 7 [marked-periodicity session](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_bott_marked_periodicity_20261007d.md)
already identifies a clean distinction. For the negative-square real
Clifford convention, `Cl_(0,k+8) ≅ Cl_(0,k) ⊗ M16(R)`; Bott's stable
homotopy theorem gives a period of eight in the real case. These are
transport theorems in specified categories, not numerical repeats of
raw dimensions or integer states. See [Atiyah–Bott–Shapiro, §4 and §6](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf)
and [Bott, §1](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf).
At the elementary affine level, tensoring a word matrix with I16 is
lossless if the marked matrix corner is retained. Similarity without the
mark already fails: plus-Collatz valuation words (1,2) and (2,1) have
equal multipliers 9/8 and similar matrices but translations 5/8 and 7/8,
and their native source classes differ (11 and 9 mod 16 respectively).
The prior note also records why replacing Z8 by E8 does not preserve the
theta series, even though a discriminant form is preserved by adding E8.

This is the useful question for the tournament mod-8 clock: *in which
class of physical continuations is a proposed equivalence actually
invisible?* Below is an exact answer for unmarked cyclic joins and a
quantitative obstruction to turning an invisible difference into +2.

## Literal probes recover the retained interface

Let P_i(T) count spanning ordered covers of the nonempty tournament T
by i nonempty directed paths; allow singleton paths. P_1=H, and P_i=0
above |T|. Write C3[A,B,C] for the cyclic join with every cross arc
A→B→C→A. Counting via maximal visits gives

    H(C3[A,B,C]) = Σ_w P_{a(w)}(A) P_{b(w)}(B) P_{c(w)}(C),

where w runs over module words that use all labels, and successive labels
advance around the cycle. Within each module the visit order is already
specified. This is the inherited [run-word path-cover kernel](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
Define the literal probe F_1(T)=H(T) and, for r≥2,

    F_r(T)=H(C3[T, TT_{r−1}, TT_{r−1}]).

Then F_r is triangular in P_1…P_r. Between two visits to T there is a
visit to each filler. No T-visit count exceeds r. At exactly r visits
both fillers must appear exactly r−1 times, in the unique word starting
and ending in T. All vertices of a filler are used as singleton runs,
so the top coefficient at P_r is ((r−1)!)². This proves:

1. F_1,…,F_r separate precisely the same pairs as P_1,…,P_r over the
   integers. They are physically evaluable, not merely formal probes.
2. Equality of P_1,…,P_r is necessary and sufficient for all cyclic
   continuations with at most r possible visits to the source (necessity
   uses the displayed probes and the bare source). In particular
   equality of the full sequence P_i is equivalent to agreement under
   every such continuation.
3. For a positive-integer modulus m, equality of the probe values
   mod m need not invert this triangular matrix by division. The
   determinant in ordered-cover coordinates is the product of
   ((j−1)!)² up to r. In unordered-cover coordinates the j-th
   diagonal is j!((j−1)!)². If the relevant determinant is a unit
   mod m, equality in that coordinate system lifts; otherwise a
   modular inverse is not automatic.

If the source order n is known, fewer probe observations can suffice.
The polynomial
`q_T(t)=Σ_i binom(t,i)P_i(T)` is the number of covers in t labelled
optional slots. For each permutation π, require nondecreasing slot
numbers, with a strict increase at each of its b(π) backward
adjacencies. This gives

    q_T(t) = Σ_π binom(t+n−1−b(π),n).

It is monic of degree n. Generalized binomial reciprocity and reversal
b↔n−1−b give q_T(−t)=(−1)^n q_T(t); also q_T(0)=0. There are exactly
⌊(n−1)/2⌋ potentially free coefficients below the monic term. After
subtracting the monic term and dividing by t for odd n or t² for even n,
one interpolates in t² at that many distinct positive slot values.
In particular H and q_T(2) recover the full ordered interface at known
orders through six. This is a dimension bound, not an assertion that
all coefficient vectors occur. Related prior work: [THM-3134 endpoint
jet and reversal](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md)
and the [Rédei–Berge power-sum formula](https://arxiv.org/abs/2307.05569).

## Factorial leakage barrier: high-order visibility alone cannot give a unit repair

**Proposition.** Suppose two source tournaments A,A' agree on
P_1,…,P_r for some r≥1. Then for every pair of nonempty filler
tournaments B,C, without order restrictions,

    H(C3[A',B,C]) − H(C3[A,B,C]) ≡ 0 (mod (r!)² (r+1)!).

**Proof.** In the run-word formula the terms with at most r visits to
the source cancel exactly. In every remaining term, the source has
i≥r+1 visits and each filler has at least i−1≥r visits. Ordered path
covers come in free orbits under permuting their distinct nonempty
components. Hence P_j(X) is divisible by j! for any X, and so is the
difference of two such values. Every summand in the difference is
divisible by (r+1)! from the source and r! from each filler. Sum. ∎

Thus once a source is an exact match on **both** Hamilton counts and
ordered two-path covers (equivalently H and q_T(2)), *no* number of
unmarked cyclic-reader trials can uncover a hidden higher-order
difference with step +2. Already r=2 forces step divisible by 24
(and hence by 8); r=3 forces divisibility by 864. Increasing reader
size to reveal higher jets does not recover the lost first bit.
There is an elementary parity refinement. Dividing P_i(T) by i! and
summing over unordered partitions into i nonempty subsets, Rédei's
parity for each induced subtournament gives

    P_i(T)/i! ≡ S(|T|,i) (mod 2),

with S the Stirling number of the second kind. If this Stirling parity
agrees at i=r+1 between the two sources (in particular, at equal
order), the first surviving source coefficient gains a factor two.
Every i≥r+2 already contributes the additional even factor
(r+2)(r+1)² in the product of the three factorials. Thus in this case
the displayed divisor doubles.

The r=1 boundary is genuinely different: both orders ≥2 have odd
S(n,2)=2^(n−1)−1, so all unmarked cyclic readers forbid a step of ±2
even with only H matched. But source TT1 versus TT2 have equal H=1,
and singleton fillers give the step 3→5. Both are literal strong
tournaments; the singleton exception is real.

The mechanism generalizes. In a fixed tournament skeleton Γ with
unmarked uniform substitution, suppose every directed closed walk
from a marked module back to itself, with no internal occurrence of
that module, visits each member of a fixed set J at least once. If
source candidates agree through r visits, surviving terms force at
least r visits to each of those mandatory return modules. Therefore
every count difference is divisible by `(r+1)! (r!)^|J|`. If no module
is individually mandatory, this argument claims only (r+1)!; alternate
return routes do not imply the stronger factor. The proof is termwise
using the same run-word expansion valid for any tournament skeleton.

## Sharper theorem: exact matched moments make a factorial band

Use the principal specialization (each of t variables is 1) of
[Grinberg–Stanley, Theorem 1.39](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge.pdf).
Here a fundamental quasisymmetric term counts exactly the nondecreasing
optional slot assignments, so the specialization of the tournament
symmetric function is q_T(t). Their tournament power-sum sum is over
permutations with only directed odd cycles and fixed points; with
p_j=t, a term is t raised to the number of cycles, times two per
nontrivial cycle. In particular q_T(t) has integer coefficients, the
same parity as its order, leading coefficient 1, and every other
coefficient even. This coefficient divisibility is stronger than
merely knowing the parity from reversal.

**Theorem.** If |A|=|A'|=n and P_1,…,P_r agree, let ε=1 for odd n
and ε=2 for even n. All unmarked cyclic readers satisfy

    H(C3[A',B,C]) − H(C3[A,B,C])
         ≡ 0 (mod 2 (r+1)^(ε−1) (2r+1)! (r!)²).             (∗)

If A,A' have different orders of the same parity, the same assertion
holds without the leading 2.

**Proof.** With Δq=q_A'−q_A, the zero match on slot labels 0,1,…,r
and parity force roots 0,±1,…,±r, with a double zero root in the even
case. Same order cancels the monic terms, giving integral coefficient
content divisible by 2; different orders of the same parity keep
integer coefficients without that factor. All linear root factors are
monic, so in the same-order case exactly

    Δq(t) = 2 t^ε ∏_(j=1)^r (t²−j²) R(t),   R ∈ Z[t].

The product with only one factor t is (2r+1)! binom(t+r,2r+1).
The second factor t in the even case brings a fixed factor r+1 because

    t binom(t+r,2r+1)
      = (r+1)[binom(t+r+1,2r+2) + binom(t+r,2r+2)].

Thus every nonnegative-integer slot evaluation has the claimed source
divisor. The cover difference ΔP_i is the i-th forward difference
of Δq evaluated at 0 (the coefficient of binom(t,i)), hence retains
that divisor for all i. For every nonvanishing run-word term each
filler still contributes a factor r!, proving (∗). The argument
for differing orders of the same parity omits the content-2 step. ∎

For example, r=2 gives divisibility by 960 at equal odd orders,
and by 2880 at equal even orders. The modulus includes more power
of two than the elementary orbit argument alone. An explicit order-nine
pair shows the **2-adic** part of this statement is sharp at r=2 for odd
orders. Number vertices 0,…,8, list pairs i<j lexicographically from
least to most significant bit, bit 1 giving j→i. Codes 38858915863
and 15014846659 have common P_1=375, P_2=10242, so common
q(2)=10,992. Their P_3 are 95562 and 98202, respectively. The
TT2,TT2 probe separates by exactly 10,560=11·960, with v₂=6.
This certifies that even at fixed n the two-coordinate shortcut
fails by order nine, not that nine is minimal. Replay
`python bott_check.py` independently counts both thirteen-vertex
cyclic joins, not just the formula.

## Consequence for exact placement

The Bott analogy is **marked equivalence**, not the number eight:
the inverse/corner dictates which tests commute with the transport.
On a tournament cyclic join the factorial orbit explains why one
pattern works on the first hidden coordinate but breaks on the next:
all repairs in this unmarked class are forced onto a coarse lattice.
A unit-step construction that wants to keep H and P₂ fixed must break
the interchangeable run order, for example with marked endpoints and
nonuniform cross arcs, rather than simply add more unmarked visits.
The prior six-site mod-8 bank uses literal distinguished choices. The
factorial barrier explains what an alternative gadget must avoid; it
does not prove a marked gadget hits every absolute target, nor that
Bott periodicity or the Cayley–Dickson tower implies Chan–Pak.
