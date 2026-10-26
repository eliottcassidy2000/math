# Native machines, exact clocks, and a configuration carrier

2026-10-26 UTC. **PROVED below:** elementary native block identities and
configuration encoding. **FINITE-EXACT:** `check.py` controls.
**CITED, not locally rebuilt:** the external Goldbach Lean proof.
**OPEN / no transfer supplied:** implications between classical Collatz,
Goldbach, RH, P versus NP, and black-core connectivity.

Closest inherited mechanism is the exact positive-duration acceleration and
quotient criterion in
`../../session25_spectrum_wiring_20261021/session25-paper/PAPER_CONNECTIONS.md`
sections 1–2. Its raw-machine correctness and domain-closure requirements
already address the essential issue; we do not rename them as new results.
The canonical hostile is a control graph that forgets write symbols.
The corrected near miss here is an external exponential-bell transcript.
The least-used sidecar is the actual ordered tape, including head, origin,
and the native clock. The concept board is finite program / unbounded tape /
macro clock / halting quantifier / marked tournament / uniform reduction.
Anchor is exact machine identification; niche is clock certification;
wildcard is the distinction between a fixed conjecture and a complete family.

## Named examples and source status

The practitioner-maintained [nonhalting-machine atlas](https://wiki.bbchallenge.org/wiki/Non-halting_Turing_machine)
identifies these binary tables:

| Name | Exact table | What the name certifies |
|---|---|---|
| Spaghetti sandwich | `1RB1LC_1LA1RD_1RA0LC_1LB0RD` | An observed central irregular region between regular growing patterns |
| Typical irregular bell | `1RB0RC_1LC1RA_1RA1LD_0LC0LA` | A dynamical description, not a general classification theorem |

Both tables are total: every state/symbol pair has a transition. Their
nonhalting is therefore immediate on every tape. It is their dynamical
behavior, not the existence of an eventual halt, that remains interesting.
No exact macro clock for either is certified in this packet. A finite visual
prefix cannot certify membership in a permanent dynamical class.

The [Bouncer page](https://wiki.bbchallenge.org/wiki/Bouncer) supplies
`1RB0LB_1LA0RA`; the [Bell page](https://wiki.bbchallenge.org/wiki/Bell)
supplies the two three-state examples below. Those pages are first-party
community analyses, not formal proof artifacts. Their dates and the separate
formal bouncer-verifier source are recorded in `SOURCES.md`.

## Native clock convention

A configuration is `(q,h,U)`, with state q, integer head position h, and
finite set U of cells containing 1. Blank is `(A,0,empty)`.
`wRq` writes w, moves one cell right, and enters q; L is analogous.
A token such as `C1` in a transcript denotes the state and symbol read
*before one native transition*. A macro step means a full named transcript,
not one native transition. A missing rule attempts the final halting step;
the test uses `None` thereafter. No named total machine has a missing rule.

### Bouncer: a quadratic cumulative clock

Let B_n have state A, head 0, and ones at `-4n,-4n+2,...,-2`, with n>=0.
Then `1RB0LB_1LA0RA` executes

    A0 (B0 A1)^(2n+1) B0 (A0 B1)^(2n+2)

in exactly `8n+8` native steps and reaches B_(n+1), translated two cells
right. This follows directly by the two-cell shift rules for `B0 A1` and
`A0 B1`; the end bounces supply their two outer tokens. From blank, the
k-th checkpoint occurs at `4k(k+1)`, with 2k ones and head position 2k.
This is an inherited example, independently checked here, not a new theorem
about all machines called bouncers.

### Cubic bell: exact return to the same head position

For `1RB0LC_1LA1RB_0RA1LC`, let D_n have state A, head 0, and ones at
`1,...,n`. For fixed n and k=n,n-1,...,1, the intermediate invariant is

    q=A, h=0, U={1,...,k} union {k+2,...,n+1}.

The transcript `A0 B1^k B0 A1 C1^k C0 A1 C0` replaces U's element k by
k+1, retaining all other cells and returning to `(A,0)`. This is verified
symbolically by the right pass over k ones, the reversal at its following
zero, and the left pass over the same k positions. Its length is `2k+6`.
At k=0, `A0 B0 A1 C0` fills cell 1. Consequently

    D_n -> D_(n+1),    duration n^2+7n+4.

Summing gives the blank-tape checkpoint time
`t_k=(k^3+9k^2+2k)/3`. The macro counter is k; the native clock is cubic.

### Exponential bell: correction with an all-order tape invariant

Use precisely `1RB0LC_1RC1RA_1LA0RB`. Let C_n, n>=1, have state A,
head 0, and U=`{-2n+2,...,-1}`. Empty U at n=1 is the blank tape.
The correct identity is

    C_n -> translated C_(2n+2),
    translation 2n+3,    duration 4n^2+10n+5.                 (E)

The literal transcript is

    A0 B0
    product_(k=0..n-1) [(C1 B0)^(2k) (C0 A1)^(2k+1)]
    (C1 B0)^(2n) (C0 A1)^(2n) C0 (A0 B1)^(2n+1).

Here is an all-order proof, separate from the finite checks. After `A0 B0`,
and before loop k, the head is h=2-2k in state C, and the ones are exactly

    {-2n+2,...,h-1} union {h+2j : 0<=j<2k}.                 (I)

For k=0 the second set is empty. `C1 B0` moves right twice, changing `10`
to `01`; its 2k repetitions move to the zero at 2+2k. `C0 A1` moves left
twice, changing the encountered left/right pair `10` to `01`. Its 2k+1
repetitions are legal because k<n leaves the final left one in the solid
interval. They produce (I) for k+1. At k=n the solid interval is empty.
The next 2n right-moving pairs and 2n left-moving pairs leave a zero at
h=2-2n and ones at h+2,h+4,...,2n+2. `C0` adds the one at h and moves
to A at h-1. The last 2n+1 `A0 B1` pairs fill the intervening zeros.
The result is the solid block `{-2n+1,...,2n+2}` and A at 2n+3, proving
the target in (E). Counting tokens proves its stated duration.

Thus from blank,

    n_k=3*2^k-2,
    |U_k|=6*(2^k-1),
    h_k=6*(2^k-1)-k,
    t_k=12*4^k-18*2^k+k+6.

**External-source correction, not silently inherited.** The Bell page as
read on 2026-10-26 gives the recurrence `a_(k+1)=2a_k+3` and uses
`(C0 B1)^(2n)` in the middle of its displayed transcript. For this exact
table and the given a_0=1 both are false. The first mismatch in the
transcript is its fifth token: the native machine reads C1, not C0.
After the correct first 19 steps it has six ones, hence C_4, not C_5.
Changing just the transcript pair does not repair the claimed recurrence.
The correction is scoped to that table, start, and parameterization.

## Goldbach25: a quantified correctness theorem, not a conjecture proof

[Yijun Leng's repository](https://github.com/lengyijun/goldbach_tm) gives
a 25-state, two-symbol machine. The two implications in
[`Content.lean`](https://raw.githubusercontent.com/lengyijun/goldbach_tm/main/GoldbachTm/Tm25/Content.lean)
establish, for the designated blank tape,

    (exists native step i: nth_cfg i = none)
      iff (exists n>=0: 2n+4 is not a sum of two primes).

Goldbach therefore corresponds to **nonhalting**. The formalization was
inspected, not compiled in this session. It cannot be called a locally
verified Lean result. The native table has one missing transition, state
24 reading zero, with states numbered 0..24; the state number is not a
macro counter. Its source proves the first canonical checkpoint at step
45: state 21 reading zero, with four ones immediately to the left.

Our independent native interpreter finds successive even-number checkpoints

    even N:    4,   6,   8,   10,   12,   14,   16,   18,   20,    22
    steps:    45, 132, 301,  596, 1101, 1988, 3263, 5160, 7805, 11522.

These are finite-exact clock measurements, not a universal formula or a
proof of Goldbach. The source table used by the offline checker was read
from `GoldbachTm/Tm25/TuringMachine25.lean`, SHA-256
`9fa72dad7c4b903747ab81debe6cc1a068e9cf548fe1a10ab0b33375e791e9e0`.
Its 49 nonhalting rules were transcribed with the missing 24/0 rule intact.

The current [Cryptids semantics source](https://github.com/google-deepmind/cryptids/blob/main/Cryptids/Semantics/README.md)
also reports a 21-state reduction and proves the same mathematical halting
condition by equivalence to Goldbach25. This is an externally reported
formal result, not an independent build here and not evidence of Goldbach's
truth. It makes “25 is the smallest known” an inappropriate claim.

## Configuration tournament with exact one-step semantics

The root packet separately encodes a **finite program table** in a marked
tournament. Here the different object is an **unbounded sequence of complete
configurations**, using the program as its fixed transition sidecar.

For `(q,h,U)`, take the finite integer interval spanning `U union {h,0}`
and one blank guard at each end. Orient each cell i toward each cell j
exactly when i<j. Color each cell with its bit, a head-state marker if
present, and a marker for coordinate zero. This is a transitive tournament
with a finite color alphabet for each fixed machine.

**Encoding theorem.** Color-preserving isomorphism classes of these
canonical marked tournaments are in bijection with finite-support machine
configurations. The order has a unique linear listing. Its origin marker
recovers integer coordinates; the other colors recover q,h,U. Conversely
these data determine the canonical interval and every arc and color.
The native operation writes the head cell's bit, moves its marker to the
immediate order predecessor/successor, changes its state, and restores one
blank guard at each end of the new active interval. This exactly commutes
with one native transition, including halting. It is a deterministic
graph-rewrite encoding; its size can grow without bound.

The pairwise observable is tape order, not a fabricated competition between
machine states. The colors, origin, canonical trimming, and transition
table are load-bearing. All underlying tournaments are transitive and have
exactly one directed Hamiltonian path. Forgetting colors therefore loses
the computation even if size and H are retained. A bare tournament repair
or isomorphism does not automatically commute with the rewrite. This is
neither a new A000568 quotient nor a black-core connectivity reduction.

The two-state write-symbol hostile used by the root packet is independently
checked here. The same read-labelled next-state graph and movement labels
can give a literal two-step blank-tape loop or a halt on step three. Thus
even a more informative graph than the unlabelled state graph is insufficient.

## Reproduction and limits

From the repository root:

    python 05-knowledge/results/root27_machine_bridges_20261026/machines/check.py
    python -O 05-knowledge/results/root27_machine_bridges_20261026/machines/check.py

The standard-library checker imports no previous engine and uses explicit
exceptions. It checks 194 parameterized native block identities (bouncer
and cubic n=0..64, exponential n=1..64), including each literal read token;
640 exhaustive one-step configuration controls on five tape cells/two
states/five head positions/two hostile tables; 6,000 further native
configuration checks across six named machines; and ten Goldbach checkpoints.
It separately rejects the two external exponential-bell errors.
These finite scopes supplement the proofs above. No sequence from an
unspecified user diagram is asserted to match these sequences.
