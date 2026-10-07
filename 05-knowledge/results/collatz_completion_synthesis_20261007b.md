# Completion-first research with explicit repair obligations

**PROVED:** the guarded two-anchor exponent phases, smaller-child receipts,
fixed-source deletion request interval, and the finite-core/rank mechanism
for digit-factorial dynamics. **FINITE-EXACT:** the explicitly bounded
computations linked below. **CONDITIONAL:** Collatz ROOT completion when a
required smaller child's first-hit certificate is supplied. **OPEN:** a
rule or well-founded family covering every positive integer.

The main advance is actual coverage beyond a declared existing entry bank:
explicit phases in every long first-reset shell admit a paid four-bit
deletion. Their children have reset length two. The completion-first idea
also becomes an exact research method: retain the artificial edges of a
completed model, then compute and discharge their defects at the original
source. The digit-factorial examples give a fully solved test of that method.

## 1. Inheritance and the live board

The starting point was [child closure after ternary reduction](collatz_child_closure_synthesis_20261007.md):
paid exits for guarded child families, genuine supplied ROOT certificates,
and a failure of nonzero-phase canonical reseeding. During the session,
incoming commit `2ba93dab07b1` supplied the two-anchor mechanism in
[THM-4601, two-anchor reduction](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md),
building on [THM-4600, run transparency](../../01-canon/theorems/THM-4600-anchored-debts-are-run-transparent-the-borel-torus-picture-of-pair-chains.md).
We independently rederived the exact affine identities consumed here;
the broader incoming hypotheses are not promoted by this session.

The anchor is **authenticated coverage at an unchanged source**. The
niche is a completely classified digit-factorial model. The wildcard is
Kolakoski self-description as an exact two-clock compiler. The live board
is source, ordered carry, guard, rank, terminal, and marked repair.

The closest proved mechanism is a common-future receipt with a strictly
smaller child. Canonical hostiles are the unpaid `35 -> 23 -> 15` inverse
chain, the factorion fixed point 145, and an infinite word whose every
finite cylinder is populated but whose positive-integer intersection is
empty. The corrected near miss is to replace a supplied source by a more
convenient member of its arithmetic family. The least-used sidecars are
the allowed deletion-depth interval and the exact defects of a modified
completed dynamics.

## 2. New guarded coverage, followed by a change of residual type

Write `M_K=2^K-1`. For each `J>=3`, the
[phase compiler](collatz_completion_anchor_20261007b.md) gives two unique
binary exponent phases:

| Actual post-run head | Shell | Exponent period | Relative size within shell |
|---|---|---|---|
| `(10)` | `v2(K-1)=2J-1` | `2^(2J+8)` | `1/256` |
| `(1,10)` | `v2(K-1)=2J-2` | `2^(2J+9)` | `1/1024` |

Thus there is an explicit phase in every shell `v2(K-1)>=4`, with a
source-preserving common-future receipt from `M_K` to `M_(K-4)<M_K`.
The guards are exact modular equations, including the bit that distinguishes
valuation 10 from a larger valuation. Positivity and first-ROOT legality
hold for every member, not merely the materialized controls.

The comparison bank is **all nonempty words of arbitrary length** over
the three native inverse generators `G1,G5,G17`. On a Mersenne source their
possible first letters require `K=5 mod6` or `K=733 mod1458`; `G1` never
starts. Exactly `485/729` of each new head phase admits neither first
letter. This is not a comparison against all historical Collatz rules,
and it is not the fraction of an entire shell. The new guaranteed parts
of the shells are respectively `485/186624` and `485/746496`.

Two fully explicit progressions outside that bank are

\[
K=18273+11943936s,\qquad K=6129+23887872s,\qquad s\ge0.
\]

The four-bit child changes the residual type: every compiled exponent has
`K=1 mod16`, so `E=K-4` satisfies `v2(E-1)=2`. An arbitrarily long two-run
has therefore become a child with exactly **two** twos. On the two displayed
progressions, the child also admits `G5`, giving

\[
g=\frac{2^{K-1}-13}{9}<M_{K-4}<M_K.
\]

Its actual forward prefix `(1,2)` reaches `M_(K-4)`, so the receipts compose
without routing back to the original source. Grounding the child remains
a separate obligation. An infinite schema of children is not a finite
list of already certified seeds.

## 3. Sufficient depth is an exact request, not an available rule

Suppose a legal inverse word at parent `m` gives child `h<m`. The
[fixed-source request theorem](collatz_completion_paper24_23_20261007b.md)
determines exactly when the numerical deletion target
`z_D=(m+1)/2^D-1` is positive, odd and smaller than `h`:

\[
\operatorname{bitlength}\!\left(\left\lfloor\frac{m+1}{h+1}\right\rfloor\right)
\le D\le v_2(m+1)-1.
\]

The interval is nonempty iff `2*oddpart(m+1)-1<h`. An actual common-future
receipt from `m` to the chosen `z_D` is still required. At the unchanged
parent 35 and child 15 the interval is empty: its only positive odd
deletion target is 17. Increasing the requested depth cannot repair this.
A genuinely supplied ROOT proof for 35 does ground 15, so changing the
certificate type can succeed where the deletion type cannot.

On a fixed native word family, the least native member gives the sharp
uniform numerical budget. In particular, `G5^36` and `G17^64` require
seven bits rather than six. That is a precise request for a new receipt,
not an assertion that a seven-bit deletion has been authenticated.

## 4. Complete a modified model, then retain every defect

The [factorion model](factorion_completion_defects_20261007b.md) proves a
global entry bound for the base-b digit-factorial map `F_b`. With
`M=(b-1)!` and least `D>=2` satisfying `b^(D-1)>D*M`, every source enters
the forward-invariant core `1..(D-1)*M`. Above it, digit count strictly
decreases. Within it, a finite graph rank completes the all-source proof.

Digit histograms provide a cycle-preserving compression because
`F_b=E o H` and the histogram map is `H o E`. This preserves the whole
future after one step; the original source and that evaluated first edge
are retained when emitting a certificate. The same quotient is not
available for Collatz valuation multisets, because changing their order
changes the affine carry and native guard.

The complete portraits, independently checked on the histogram carrier
and the full integer core, are:

* Base 6: exactly four fixed cycles `1,2,25,26` (decimal labels).
* Base 10: fixed cycles `1,2,145,40585`, plus `(169,363601,1454)`,
  `(871,45361)`, and `(872,45362)`.

Now change one edge on each non-root cycle to point to 1. The resulting
map `G` reaches 1 from every positive integer. Let `V` be its exact
remaining time. For the original map define

\[
\delta(n)=\mathbf1_{n\ne1}+V(F_b(n))-V(n).
\]

This residual is zero except at the changed edges, where it equals the
original cycle period. Around any non-root cycle its sum is the period,
for **every** choice of potential. The defects cannot disappear through
reweighting alone. The source 145 is a decisive one-state test of an
invalid transfer from modified completion to original completion.

For Collatz the useful analogue is to start with a proposed completed
proof graph, mark every unsupported edge, and replace those edges with
actual source-preserving receipts or well-founded child dependencies.
The new phase compiler does this for specific infinite families. What is
missing globally is entry into those families and grounding of the
remaining child types.

## 5. Self-description preserves a clock, not universal realizability

The [run-clock package](collatz_self_describing_runs_20261007b.md) corrects
two conventions: the cumulative-one sequence is
[A156077](https://oeis.org/A156077); [A071820](https://oeis.org/A071820)
is the different Kolakoski sequence on `{2,3}`. The first 200,000 terms
of [A000002](https://oeis.org/A000002) have 133,321 runs and 133,320
transitions. Its ten-million-term prefix has 6,666,660 runs and 5,000,046
ones. The limiting one-half frequency remains open.

For A000002, let `E(r)=2*ones(r)-r`, let `S_r` be the end of run r, and
let `D(r)=sum((-1)^(i+1)*k_i)`. Exact self-description gives

\[
S_r=\frac{3r-E(r)}2,\qquad E(S_r)=D(r).
\]

The second expression is an alternating observer, not the original
discrepancy evaluated at the same scale. The run phase and the last
partial-run deficit must remain attached. They yield a lossless compiler
of finite valuation words and their native source cylinders.

There is an especially useful hostile. Every finite `{2,3}` valuation
prefix has infinitely many legal positive odd source realizations. But
`U(n)-1 <= (3/4)*(n-1)` on those steps, while a first arrival at 1 from
above requires a final even valuation at least four. A length-r source
must satisfy `n>=1+2*(4/3)^r`. Thus no positive ordinary integer realizes
the infinite word. Populated cylinders at every finite depth do not prove
an all-depth ordinary-source intersection. Retaining source height blocks
exactly that invalid compactness transfer.

Conversely the `{1,2}` self-describing word has an exact local expansion
certificate from its forbidden factors; this does not establish that a
positive integer realizes the infinite word. Both directions reinforce
the same target: guards and clocks must stay tied to the supplied integer.

## 6. What the four papers contributed

The supplied headlines were accepted as premises for inspiration. The
new elementary Collatz proofs do not depend on transferring those headlines.
The detailed page/diagram readings are in the linked component notes.

| Source | Useful operation | Preserved target and required sidecar | Cheapest hostile |
|---|---|---|---|
| Paper 24, semialgebraic covers | Correct a local model while fixing its base parameters | The original integer; retain its actual threshold and guard | Replacing parent 35 by another progression member |
| Paper 23, full-level cuspidal decomposition | Convert a detector into a finite witness | Actual certificate word; retain ordered data and norm/scale | An assumed support statement supplies no route |
| Paper 22, ordered matrix removal | Recheck obstructions after edits | Chronological legality, including absent patterns | Local edit can create a new forbidden configuration |
| Main 4, binary-information contraction | Track means and conditional information through a channel | Exact guard; retain the original observation model | Full-support noise cannot issue a nontrivial zero-error residue certificate |

The shared design lesson is operational: compression is useful only with
a proof of the predicate it preserves. Factorion histograms preserve the
future; Collatz ordered carries preserve an affine word; noisy statistics
preserve neither exact source membership nor a ROOT receipt by themselves.

## 7. Precise next targets

1. Compile further two-anchor heads into exact exponent guards, using the
   original source parameter throughout. Compare the union of guards,
   with overlaps removed, against the remaining shells.
2. Ground the resulting reset-length-two children or give each a smaller
   authenticated obligation. This is a concrete change-of-type target,
   stronger than asking for another copy of the original long-run rule.
3. For any proposed completion graph, compute its defect support and
   demand a decreasing rank on the actual residual obligations. A finite
   collection of supplied seed certificates is useful only after every
   generated dependency is proved to reach that collection.

No source-average, finite census, conditional child proof, or 2-adic
intersection is promoted here to universal positive-integer coverage.
