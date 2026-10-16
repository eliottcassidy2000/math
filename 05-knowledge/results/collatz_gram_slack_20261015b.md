# Joint-terminal cancellation, observed-edge closure, and paid graph glue

2026-10-15. **PROVED:** exact signed pair-fibre cancellation; a variance
comparison; monotone finite observed-graph closure; a weaker logarithmic
Gram target still sufficient for the inherited nonzero-residue certificate.
**FINITE-EXACT:** eight source-window/controller cases, including failures,
and48 extension blocks compared both with and without retained graph glue.
**OPEN:** both uniform Gram targets, universal ROOT entry, and any uniform
rate of weighted cancellation. No unknown terminal label is sampled.

[Harness](../../04-computation/experiments/collatz_gram_slack_20261015b.py),
[exact output](collatz_gram_slack_20261015b.json),
[quadratic-cycle/Pell companion](quadratic_cycle_phase_20261015b.md).

## 1. Inheritance and the live concept board

The closest proved mechanisms are T2 of
[terminal mesh](collatz_terminal_mesh_20261015a.md), M1--M2 of
[clock mesh](collatz_clock_mesh_20261015a.md), and the checked edge union in
[the terminal-basis construction](collatz_terminal_basis_20261007.md).
We specialize the existing observation-union mechanism to the actual
coarse-current Gram matrix; graph closure itself is inherited.

The hostile is a freely moving lone terminal, which earns zero
cancellation forever. The corrected near miss is treating a stationary
covariance as a pointwise native relation. The least-used sidecar is the
whole authenticated edge union from unsuccessful prefixes, with native
step parity. Retaining just the latest endpoints can discard useful
observations even when the latest Gram matrix is exact.

Concept board: marked coarse means; signed pair collisions; observed
interior edges; linear versus quadratic phase; Q--Pell scale; Bellman slack.
Anchor: source-specific weighted cancellation. Niche: rational six-cycles
and switching generators. Wildcard: enlarge the representation language
while keeping the original observable recoverable, as in the attached
unit-fraction manuscript. Physical source positions and integer interval
lengths never change when native controller depths change.

Incoming at `6e3b04d2f8` initially reserved THM-4658 and HYP-9282/9283.
The later pre-audit checkpoint `378a9b6390` supplies the exact two-octave
identity, a separate conditional logarithmic tail certificate, and a
heuristic against constant budgets. The [incoming-work audit](collatz_incoming_audit_20261015b.md)
repairs one missing density assumption, constructs a mean-zero branch
correlation hostile, and gives an annular Gram interface to that identity.
Our G2 implication does not depend on the new interval run. The incoming
tournament correction MISTAKE-617 reinforces the audit rule: retain every
failing row and state exactly which population was checked.

## 2. The exact Bellman slack is signed overlap, not collision count

Keep the fixed physical partition I, its lengths L_I, and the complete
terminal rows v_I, including ROOT with label1. As in T2,

    G=sum_I v_I v_I^T/(4L_I),    C(P)=z^T G z.

A common authenticated terminal controller maps old terminal e to a new
terminal h(e), with accumulated sign s_e=+/-1. Its signed identification
matrix A has one nonzero per old row, and G'=A^T G A.
For each **ordered** new pair (u,v), let P_uv and N_uv be respectively
the sums of positive and negative absolute contributions
s_e s_f G_ef with (h(e),h(f))=(u,v). Then

    ||G||_1-||G'||_1 = 2 sum_(u,v) min(P_uv,N_uv).          (1)

This is a general real-matrix identity, not a statistical approximation.
If an authenticated non-ROOT cycle makes a terminal's ROOT-supported
label zero, deletion adds exactly the absolute mass of every old ordered
entry touching a deleted column. Such deletion is additional to (1).

**Proof.** Every old entry goes to precisely one pair fibre, with unchanged
absolute magnitude. In that fibre its old mass is P+N and its new mass
is |P-N|. Their difference is2min(P,N). Deletion and transport partition
the old entries. ∎

Thus equality holds precisely when each surviving pair fibre is
one-signed (and no positive mass was deleted). Merely counting collisions,
merged sources, or ROOT hits cannot lower-bound this loss.

Define the conditional sign-variance invoice

    V_pair=sum_(P+N>0) 4PN/(P+N).

Fibre by fibre, with ell=2min(P,N),

    V_pair/2 <= sum ell <= V_pair.                        (2)

Indeed 4PN/(P+N)=(P+N)(1-m^2), m=(P-N)/(P+N), while
ell=(P+N)(1-|m|). Their ratio is1+|m|. The constants are sharp:
balanced fibres give the factor2, and very unbalanced fibres approach1.
The computation checks these inequalities by integer multiplication;
a dyadic enclosure stores the variance without enormous denominators.

For a sequence of blocks, (1) telescopes exactly. A usable arithmetic
lemma would bound their *weighted* opposite-sign overlap in terms of the
marked source rows. Independent random signs are not such a lemma.
Entropy also cannot be substituted at unit cost: for a minority fraction
q tending to0, binary entropy divided by2q is unbounded. This explains
why the attached Bellman paper's exact slack decomposition is a useful
method while its entropy constants do not transfer automatically.

There is also a native **complementarity hostile** to greedy payment.
Take the test row t(3)-t(8) and the two possible observations
3->5 and5->8. With neither edge, or either edge alone, its unnormalized
Gram norm is4. With both edges, t(3)=t(8), so the row and its Gram matrix
vanish. The cancellation value of an edge set is therefore not generally
submodular: two zero-gain observations can have positive joint gain.
More generally, t(2^m)-(-1)^m t(1) has the same all-or-nothing behavior
when the only available edges are a subset of its m-edge halving chain.
These are exact native *test rows*, not the actual coarse-current rows;
they refute a generic greedy-contraction argument while motivating a
controller that can pay for a whole certified block.

## 3. Close the observed graph before forgetting its interior

Fix a finite authenticated graph of native shortcut edges, stop at ROOT,
and keep all the original source weights. Each finite functional component
ends at a missing outgoing edge or a directed cycle. Compile each source
to that sink, retaining (-1)^d and its exact integer native depth d.
A component reaching ROOT is grounded. A witnessed cycle excluding ROOT
has ROOT-supported label zero, for either cycle parity. A missing edge
remains an unconstrained terminal. No mutually conditional cycle is
accepted as a ROOT proof.

**Theorem G1 (PROVED).** This compilation preserves every coarse source
mean exactly. Adding authenticated edges and recompiling cannot increase
the joint entrywise-l1 Gram bound. Among forward-prefix controllers using
only those edges and the witnessed non-ROOT-cycle zero rule, the completely
compiled terminal partition gives the smallest such norm, up to irrelevant
terminal renaming and sign gauge.

**Proof.** Replay each finite path, multiplying its edge signs. All sources
that reach the same sink have the same remaining variable; known cycles
have the stated support rule. A partially compiled endpoint has a further
path to its component's compiled representative, or to a certified-zero
cycle. Therefore complete compilation is a signed identification/deletion
of any partial compilation, and (1) applies. Adding edges is another such
identification/deletion of the old boundary variables. ∎

The minimality is within this representation; solving the terminal cube
or adding other arithmetic constraints may give a smaller bound still.
The theorem does not assert that every missing edge eventually reaches
ROOT. It applies to every finite source graph whether Collatz is true or
false.

For the experiment, initially observe D native steps from **each** source
1<=n<2X, stopping early at1. Freeze their union, then close it. Next observe
q steps from each nonzero terminal column, q in{1,2,4,8,16,32}. Compare:

* transport the terminal rows synchronously by those q steps;
* use the **same new observations**, reinsert them into the retained old
  graph, and close it before transporting the rows.

The second procedure can join trajectories at different times, or reenter
an old ROOT-connected path. It does not invent edges or assign labels to
unresolved endpoints. It uses a longer certified word assembled from
already checked pieces. The physical integer clock in G is unchanged.

## 4. Finite results, with failures and costs retained

All tests use b=2^floor(k/2), X=2^k, and the exact fine-scale variance
allowance V_*. The last column tests the full energy inequality E_X<=2X^2
by the rational square test from T1, not by decimal rounding.

| k | Initial D | Coarse ||G||_1/X^2 after graph closure | Initial energy test |
|---:|---:|---:|---|
| 10 | 8 | 0.477201343 | FAIL |
| 10 | 16 | 0.113808930 | PASS |
| 12 | 16 | 0.191725469 | PASS |
| 14 | 16 | 1.455857909 | FAIL |
| 14 | 24 | 0.618031288 | FAIL |
| 14 | 32 | 0.472391431 | FAIL |
| 16 | 32 | 0.190550115 | PASS |
| 18 | 36 | 0.694866316 | FAIL |

A failure is failure of this sufficient certificate, not failure of the
actual energy inequality. Independent finite ROOT grounding verifies the
actual energies in all these windows. For example at k18 actual S_X/X^2
is0.068484289, much smaller than the coarse upper bound.

Three illustrative repairs using the retained graph are:

| k,D | New q per terminal | Distinct new edges | Coarse bound/X^2 with graph glue | Energy test |
|---|---:|---:|---:|---|
| 10,8 | 1 | 35 | 0.314892918 | PASS |
| 14,32 | 1 | 50 | 0.425711329 | PASS |
| 18,36 | 8 | 3411 | 0.339344099 | PASS |

For k14,D32, synchronous transport keeps0.472391431 through q16 and
still has0.463470 at q32, failing the energy test. One-step graph glue
already passes. For k18,D36, synchronous q32 remains0.679211 and fails;
graph glue at q8 passes. The effect is reuse of unequal-time path joins,
not a stronger assumption about how fast an individual orbit terminates.
The latter certificate still leaves1,908 original sources unresolved;
the k14 one-step repair leaves265. Their nonzero terminal weights remain
in the matrix throughout the successful acceptance test.

Every initial case has zero one-step synchronous loss. Several have
zero loss through16 or32 steps. Therefore a strict uniform one-step
contraction is not supported even by these finite source-weighted data.
Graph glue changes which terminal relations are available and is a
productive route around that obstruction. A universal lower bound on the
newly available weighted overlap remains missing.

The output records the exact Gram norm, small full matrices, hashes of
larger matrices, graph/source/row hashes, both sets of block losses, costs,
terminal counts, unresolved source counts, and native word lengths. The
harness regenerates every matrix; larger matrices are not treated as
uninspected external input.

## 5. A weaker target still pays the amplitude certificate

The old all-scale target ||G_k||_1 <= (4/9)X^2 is sufficient, but need
not be the only useful target. The incoming reserved logarithmic-budget
idea prompts the following unconditional *implication*.

**Theorem G2 (PROVED implication; premise OPEN).** Suppose that for every
integer k>=34, finite authenticated controllers at X=2^k produce

    ||G_k||_1 <= (k^2/100) X^2.                            (3)

Then the nonzero-residue conclusion (8) of
[the amplitude-boundary note, Theorem B3](collatz_amplitude_boundary_20261013c.md)
holds, with its existing interval certificate unchanged. In particular
(3) is a sufficient target permitting a growing energy budget, rather
than a constant budget at every scale.

**Proof, energy step.** The fine-scale theorem gives
S_X/X^2 <= k^2/100+1/27+1/(3X) <= k^2/100+1/16.
For k>=34,

    sqrt(E_X)/X <= 2sqrt(k^2/100+1/16)+X^(-1/2)
                 <= k/5+1 <= k/4.

Thus E_X <= k^2 X^2/16. Here sqrt(a^2+b^2)<=a+b was used with
nonnegative a,b; the deliberately generous constants avoid a numerical
optimization.

**Proof, Mellin step.** Put eta=sigma-1/2 and
c_sigma=sqrt((1-2^(-2sigma-1))/(2sigma+1)). Cauchy--Schwarz on each block
then gives

    integral_X^(2X) |H(x)| x^(-sigma-1) dx
       <= (k/4) c_sigma X^(-eta).                         (4)

For 99/100<=sigma<=1, c_sigma^2<=175/596<(11/20)^2,
whereas integral_1^2 t^(-sigma-1/2)dt >=2-sqrt(2)>11/20.
Consequently (4) is bounded above by

    (1/4) integral_X^(2X) sqrt(x) log_2(x) x^(-sigma-1) dx.

Summing proves exactly the permitted tail estimate (11) in B3 near its
certified zero, whose real part is0.9982092207.... The same finite
boundary term and the independent R-tail are unchanged. For every
sigma>1/2 the series sum k 2^(-k eta) converges, so analytic continuation
also follows on that entire half-plane, uniformly on compact subsets.
The previously certified normalized residue interval is therefore still
0.1193 < 2|r|2^(29(Re(s0)-1)) < 0.5908. ∎

More generally, any certified sequence E_(2^k)<=a_k 2^(2k) yields the
explicit tail invoice

    c_sigma sum_(k>=K) sqrt(a_k) 2^(-k eta).               (5)

This can accommodate bad scales when the total invoice is paid. Merely
writing O(log X), with no constants or threshold, does not certify
noncancellation. Neither (3) nor a sufficiently small invoice (5) has
been proved for all scales. Neither would prove universal ROOT entry.

## 6. How the attached papers and prior controls are used

The attached Paata Ivanisvili manuscript, *A Bellman function for disjoint
partitions: the sharp logarithmic gain for three families*, Appendix A.3,
keeps an exact slack decomposition before estimating its terms. We adopt
that research move in (1)--(2); we do not import its partition-count
inequality as a Collatz theorem or assert an audit of its interval code.

The attached Bloom--van Doorn manuscript, *Optimal bounds for the length
and number of unit fraction representations*, section2.1--2.2, enlarges
integer divisor sums to fractions with controlled smooth numerators and
an extra denominator. The target representation remains recoverable.
Our analogous relaxation is equal-depth terminal transport -> a retained
finite graph with unequal-time joins. The preserved predicate is the exact
source mean; the sidecar is the actual native edge receipt and its parity.
No unit-fraction summand is itself a valid Collatz move. The manuscript's
asymptotic theorems and Lean provenance are not used as dependencies here.

The cycle/Pell companion supplies an especially concrete hostile: two
individually order-six generators compose to the Pell matrix after a
change of basis. Individual invariant forms do not provide a common
changing-generator budget. Q's sampling hypotheses remain open; the
actual signed-identification Gram cone requires no such hypothesis.

The existing [block/glue clock formula](clock_duality_glue_20261015a.md)
keeps the physical seam invoice. Our graph-glue operation is different:
it changes the certificate for fixed source rows, rather than concatenating
new source intervals. Conflating these two meanings of glue would change
the energy being certified.

## 7. Audit, reproduction, and next obligations

    python3 -B 04-computation/experiments/collatz_gram_slack_20261015b.py
    python3 -B 04-computation/experiments/quadratic_cycle_phase_20261015b.py

The declared initial universe is the eight (k,D) pairs in the table; the
extension universe is all six q values for every case. Native edges are
authenticated; each compiled path satisfies its exact signed edge/depth
recurrence. Small cases receive literal source-path replay. ROOT labels
for the independent audit are computed by a separate first-descent method
only after graph and row hashes are frozen. Actual source block means and
energies are checked against the certificates. Gram multiplication uses
int64 only after proving a conservative no-overflow bound; direct Python
integer dot products audit all small matrices and every large diagonal.
Every transported matrix is independently reconstructed from the rows.
No floating-point value participates in acceptance. Normal and optimized
Python runs give identical JSON outputs. All125 partial toy functional
graphs on three non-ROOT vertices and all300 one-edge extensions verify
the compiler and pair-transport agreement, including both cycle parities.
These exhaustive toy graphs are controls, not extra native Collatz cases.

The remaining mathematical obligations are now more concrete:

1. Bound the weighted opposite-sign pair overlap produced when the actual
   terminal boundary reenters retained interior paths; mere source counts
   and equal-time coalescence rates do not suffice.
2. Find a source-addressed block rule preserving both the full edge context
   and physical weights, with either constant budget4/9 or logarithmic
   budget k^2/100. Its test must include the exhibited long zero-loss
   synchronous plateaus and the changing-generator Pell hostile.
3. Estimate the explicit Mellin invoice (5) directly if separate uniform
   bounds discard too much information at rare scales.

The successful move this session is to spend already authenticated
interior information before asking for a new global inequality. Its scope
is a strictly stronger finite certificate and an exact weaker sufficient
all-scale target; the missing arithmetic overlap bound is still OPEN.
