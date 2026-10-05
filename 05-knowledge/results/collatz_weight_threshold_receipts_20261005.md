# A finite compiler for a source-specific positive weight floor

**PROVED:** the exact finite superlevel compiler, uniform-witness compactness,
and floor-to-receipt interface. **FINITE-EXACT:** the declared controls.
**OPEN:** obtaining a positive floor for every supplied positive odd source.
The inverse closure and mixture weight are inherited; no new basin is claimed.

Artifacts: [script](../../04-computation/experiments/collatz_weight_threshold_receipts_20261005.py)
and [output](collatz_weight_threshold_receipts_20261005.out).

## 1. Four less-used mechanisms and the transfer

This pass starts outside the recent Collatz/modular thread:

| Canonical mechanism | What actually transfers | What does not |
|---|---|---|
| [THM-4026, Sun counterexample](../../01-canon/theorems/THM-4026-sun-two-four-six-eight-binomial-counterexample.md) and [THM-4027, universal modular solubility](../../01-canon/theorems/THM-4027-sun-two-four-six-eight-universal-modular-solubility.md) | The distinction between compatible modular witnesses and an actual integer witness with bounded height | Universal local solubility does not imply global representability |
| [THM-4210, Rule 30 lossless Cartier tree](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md) | Retain the actual admissibility data, not only a lossless ambient encoding | The current tree does not establish physical admissibility or termination merely by being lossless |
| [THM-4286, signature-response nonfactorization](../../01-canon/theorems/THM-4286-signature-response-nonfactorization-and-two-deck-surgeries.md) | A target predicate descends through an observer exactly when every observer fibre is pure | Matching observed signatures do not transport a certificate across a mixed fibre |
| [THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md) | Match the exact labelled boundary state before gluing two proofs | Two positive side masses can have disjoint supports; the Petersen example has zero pairing |

The quantitative repair here is to require every replacement witness to
carry the **same positive weight floor**. That confines witnesses to one
finite, computably enumerated set. Refinement can then isolate the actual
source. The live board is **source / exact receipt / weight floor / finite
candidate set / refinement / independently justified floor**.

For the Sun sum, canon proves that
$N=896315812331399$ has no exact representation but has representations
modulo every modulus. If $m>N$, every such nonnegative integer sum $S$
satisfies $S\ge N+m$: the only smaller nonnegative candidate in that
congruence class is the forbidden $S=N$. Thus these local witnesses must
escape in ordinary size. This is the concrete hostile motivating the
bounded-witness condition, not a new computation of the counterexample.

## 2. Exact weight and the finite bounds

Write $U(n)=(3n+1)/2^{v_2(3n+1)}$, and stop a receipt at its first visit
to ROOT $1$. The inherited [adaptive mixture](collatz_adaptive_mixture_flow_20261005.md)
defines $W(1)=1$ and, for a nonroot first-hit word
$a=(a_1,\ldots,a_\tau)$,

\[
L=\tau-1,\qquad K=\sum_{i=1}^{\tau}\left\lfloor\frac{a_i-1}{2}\right\rfloor,
\qquad
W(n)=w(L,K)=\frac{2K!(L+1)!}{(L+K+2)!}.
\tag{1}
\]

Outside the rooted component $W(n)=0$. This definition does not assume that
every positive odd integer is rooted. The final valuation of any nonroot
receipt is even and at least four: $3x+1=2^a$ gives
$x=(2^a-1)/3$, and $a=2$ would be the forbidden ROOT self-edge. Therefore
$K\ge1$.

For $0<\epsilon\le1$, let

\[
B=\max\{b\ge0:\ \epsilon(b+1)(b+2)\le2\}.
\tag{2}
\]

The sharper total-counter bound is shared with the concurrent
[floor deadlines](collatz_floor_transport_deadlines_20261005.md). Its short
proof is included so the compiler is self-contained. Put $N=L+K$.
For a nonroot receipt, $1\le K\le N$, so

\[
w(L,K)=\frac{2}{(N+2)\binom{N+1}{K}}
\le\frac{2}{(N+1)(N+2)}.
\tag{3}
\]

Consequently $W(n)\ge\epsilon$ forces

\[
N\le B,\quad 0\le L\le B-1,\quad 1\le K\le B,\quad
\tau=L+1\le B,\quad A=\sum_i a_i\le2(N+1)\le2(B+1).
\tag{4}
\]

ROOT is treated separately. The looser deadline $\tau\le B+1$ is also
valid, but the concurrent forward checker now uses the sharper (4), which
retains the compulsory final $K\ge1$.
Each individual exponent is at most $2B+2$.

Clearing the complete actual word gives
$3^\tau n+C=2^A$, with positive integer carry $C\ge1$. Hence

\[
n\le H_\epsilon:=\frac{4^{B+1}-1}{3}.
\tag{5}
\]

Equality is attained by the ROOT sibling $S^B(1)$, where $S(x)=4x+1$.
For $B\ge1$ its word is $(2B+2)$ and its weight is
$2/((B+1)(B+2))\ge\epsilon$; for $B=0$ it is ROOT with the empty word.
Thus (5) is the exact largest source in the superlevel set, not only a
rough size bound.

## 3. The compiler and its completeness

Start with the only seed $1$. For each certified parent $p$ and positive
valuation $a$, form

\[
n=\frac{2^ap-1}{3}.
\tag{6}
\]

Retain it only when it is an integer greater than one. It is automatically
odd, and $3n+1=2^ap$ has exact valuation $a$. Prepend $a$ to the parent's
first-hit word. The new source cannot be ROOT; the parent's word first hits
ROOT at its end, so the new word also does.

For parent ROOT, the child counters are
$(0,\lfloor(a-1)/2\rfloor)$. For every other parent they are

\[
(L+1,\ K+\lfloor(a-1)/2\rfloor).
\tag{7}
\]

Use (4) to bound the candidate exponents, and retain a child precisely when
its weight is at least $\epsilon$. Every nonroot inverse extension strictly
decreases weight: increasing $L$ gives factor
$(L+2)/(L+K+3)<1$, and increasing $K$ also decreases it.
ROOT's children have weight at most $1/3<1$.

**Exact superlevel theorem.** This finite algorithm returns exactly

\[
\mathcal H_\epsilon=\{n>0\text{ odd}:W(n)\ge\epsilon\},
\tag{8}
\]

together with the unique actual first-hit word for every member.

Soundness follows from (6), strict ROOT handling and exact rational
weight comparison. Finiteness follows from (4): only bounded-length words
over a bounded exponent alphabet can be visited. For completeness, every
member of (8) has a first-hit word; deleting its first edge moves toward
ROOT and increases weight. All suffixes therefore survive pruning, so
induction backwards from ROOT generates that source. Duplicate inverse
addresses cannot occur: a positive integer has one actual $U$-parent and
one exact valuation. A duplicate along its own ancestral chain would
contradict the already verified first-hit word.

The program neither searches an unknown forward orbit nor imports any
nonroot seed. Its public **threshold_receipts(epsilon)** returns the exact
dictionary. **certificate_weight(source,word)** independently replays every
claimed exponent and rejects ROOT padding. **source_receipt(n,epsilon)**
returns the word or None; the latter means exactly

\[
W(n)<\epsilon,
\]

which includes both small positive weights and zero weights. It is not a
nonconvergence decision. The threshold must be an exact Fraction in
$(0,1]$; source and word types are checked.

The inherited bound $\sum W\le16/3$ gives
$|\mathcal H_\epsilon|\le16/(3\epsilon)$. This is a useful size estimate,
not a polynomial-time claim in the bit length of $\epsilon$. Tiny floors
can require enormous finite banks. For a single supplied source, the
concurrent bounded forward test is often cheaper than enumerating this
whole reusable bank.

## 4. A weight-qualified refinement really can lift a witness

Fix a positive odd integer $n$ and $\epsilon>0$ in the compiler's domain.
The following are equivalent:

1. $W(n)\ge\epsilon$.
2. For every $j\ge1$ there is a rooted integer $n_j$ satisfying
   $n_j\equiv n\pmod{5^j}$ and $W(n_j)\ge\epsilon$.

The forward direction uses the same source at every scale. For the reverse
direction, all $n_j$ lie in the one finite set (8), whose members are at
most $H_\epsilon$. Choose $5^j>\max(n,H_\epsilon)$. Both positive integers
then lie below the modulus, and congruence forces $n_j=n$.
In fact this one sufficiently fine scale suffices.

This is a constructive bounded-witness compactness statement. It asks for
one **individual witness weight** at least $\epsilon$, not merely a
positive sum of masses in a residue class. The latter can be positive at
every scale while the selected atom is zero, as proved by the hostile in
[representation positivity](collatz_representation_positivity_20261005.md).
Nor does a fixed modulus globally decide membership in (8): any residue
containing a heavy source also contains arbitrarily large sources outside
the finite heavy set. The supplied source and the refining modulus remain
essential.

For $\epsilon=1/1000$ the exact qualified witnesses for supplied $27$ change:

| Modulus | Largest qualifying atom in its class | Weight |
|---:|---:|---:|
| 5 | 17 | $1/30$ |
| 25 | 227 | $2/105$ |
| 125 | 277 | $1/140$ |
| 625 | none | all atoms in this class have weight below $1/1000$ |

The lack of a qualifying witness does not refute ROOT: an independently
supplied and replayed word for $27$ has weight $1/11274451650$.
For supplied $7$, refinement instead stabilizes at its own word and weight
$1/84$ already modulo 25. These are finite controls, not a new coverage
result.

## 5. Connect an independent atom lower bound to an actual receipt

The incoming injection measure is
$\lambda(n)=W(n)$ on positive odd multiples of three, and zero elsewhere.
For such a supplied source $n$, an independently proved rational inequality

\[
\lambda(n)\ge\epsilon>0
\tag{9}
\]

forces the compiler to return its ROOT receipt. This direction is
unconditional once the premise (9) is proved; it does not assume universal
Collatz.

The [polynomial atom dual](collatz_atom_polynomial_dual_20261005.md)
supplies a precise possible frontend. Its signed polynomial is bounded
above by the indicator of a selected support point, and truthful moment
intervals yield a rational lower readout. A positive readout for source
$6m+3$ is a premise of type (9). The intervals must be bounds for the
actual measure; arbitrary numerically plausible intervals are not proof.

The companion tests this interface at sources $3$ and $9$, using the exact
positive readouts from that package's independently grounded sixteen-atom
head. The returned words are $(1,4)$ and $(2,1,1,2,3,4)$, with actual
weights $1/6$ and $1/126$. This is an already-grounded interface test,
not new sources or a way to manufacture unknown moment premises.

## 6. Reproduction and the honest cost comparison

From the repository root:

    python 04-computation/experiments/collatz_weight_threshold_receipts_20261005.py
    python -O 04-computation/experiments/collatz_weight_threshold_receipts_20261005.py

The finite thresholds are $1/3,1/6,1/10,1/20,1/100,1/1000$.
Their bank sizes are $2,4,5,7,29,188$, respectively.
At the four smallest comparison searches, three independent routes agree:
the weight-pruned inverse tree, the complete rectangular counter bank
filtered afterwards, and direct whole-word enumeration followed by solving
the affine source equation and literal replay.

| Threshold | Pruned sources / inverse candidates | Rectangle sources / inverse candidates |
|---:|---:|---:|
| $1/3$ | 2 / 4 | 2 / 6 |
| $1/6$ | 4 / 8 | 5 / 18 |
| $1/10$ | 5 / 16 | 12 / 44 |
| $1/20$ | 7 / 28 | 28 / 108 |

An inverse candidate is one tested pair $(p,a)$, including failed
integrality/ROOT checks. Whole-word tests have a different cost and are
reported separately, not compared as equivalent runtime units.
The $1/1000$ bank needs 11,718 candidate inverse pairs and includes sources
up to $(4^{44}-1)/3$, an 87-bit integer.
Every returned receipt is literally authenticated; no correctness check
uses Python assert, so optimized mode retains all controls.
The companion also compares every odd source below 256, at all six
thresholds, with the independent bounded forward decision from the floor
deadline package: 768 exact interface controls.

The next admissible step is a new independently justified atom floor for
an unresolved supplied source, or a more efficient representation of this
same finite bank. Merely increasing the bank or replacing one local witness
by another does not establish universal coverage.
