---
id: THM-825
title: Quadratic positional moments are a literal mirror-layer codec through size fifteen
status: PROVED + FINITE-EXACT (word-fibre census and staircase geometry)
source: codex-2026-07-15-S13-postjoin
depends_on: [THM-553, THM-809, THM-818]
related: [THM-796, THM-814, HYP-6880]
verification:
  - 04-computation/quadratic_positional_layer_codec_codex_S13_postjoin.py
  - 05-knowledge/results/quadratic_positional_layer_codec_codex_S13_postjoin.out
  - 05-knowledge/results/quadratic_positional_layer_codec_codex_S13_postjoin.json
---

# THM-825 — quadratic positional moments close every B2 layer through `n=15`

For a word `w` on the ordered positions `0,...,s-1` with state set `C`, put

```text
M_j(c;w)=sum_{i:w_i=c} i^j,       c in C, j=0,1,2.       (1)
```

Here `M_0` is the state count already retained by the B2/S2 carrier.  The
joint per-state carrier `(M_0,M_1,M_2)` is injective on every binary or
four-state word of length at most six.  It first fails at length seven.

Consequently the B2 mirror-layer counts with their first two positional
moments reconstruct the entire literal apex-zero tiling at every staircase
size `n<=15`.  In particular this is an unconditional exact repair carrier at
every size relevant to LRC(14), independently of the outcome of THM-818's
`n=9` relation join.  This is a static reconstruction theorem, not a claim
that the carrier is continuation-minimal or preserves the LRC predicate.

The theorem was initially headed for THM-824.  A live pull showed that
THM-824 had just been reserved by the fixed-pair symmetric-radius project, so
the namespace was moved to THM-825 before publication.

## 1. The six-position subset lemma

Let `A,B` be subsets of `{0,...,s-1}`, where `s<=6`, and suppose

```text
|A|=|B|,       sum A=sum B,       sum_(a in A) a^2=sum_(b in B) b^2. (2)
```

Set `P=A\B` and `N=B\A`.  Then `P,N` are disjoint, have the same size `r`,
and satisfy the same sum and square sum.  Since `2r<=6`, only four cases are
possible.

- `r=0` gives `A=B`.
- If `r=1`, equality of sums makes the two singletons equal, contradicting
  disjointness.
- If `r=2`, the common sum and square sum determine the common product by

  ```text
  xy=((x+y)^2-(x^2+y^2))/2.
  ```

  Thus `P` and `N` are the root multisets of the same monic quadratic, again
  contradicting disjointness.
- If `r=3`, disjointness forces `s=6` and `P union N={0,...,5}`.  Equal sums
  would split the total sum `15` into two equal integers, which is impossible.

Therefore (2) implies `A=B`.

For a word `w`, apply the lemma to the position set

```text
A_c(w)={i:w_i=c}
```

for each state `c`.  Equality of (1) gives equality of every `A_c`, hence of
the literal word.  This proves the binary and four-state assertions at once;
no exhaustion is needed for the proof.

## 2. The boundary is sharp

On seven positions the disjoint subsets

```text
P={0,4,5},        N={1,2,6}                              (3)
```

have

```text
(cardinality, sum, square sum)=(3,9,41).                 (4)
```

Colour `P` versus `N` with one state and their complements with a second
state.  The complementary subsets also have common data `(4,12,50)`, so the
full per-state keys agree.  This is simultaneously a binary-word and a
four-state-word collision.

The exact word census confirms the sharp boundary.  At degree two it gives

```text
states  length  words   cells   collision cells/excess   max fibre
   2       6       64      64             0/0                 1
   2       7      128     126             2/2                 2
   4       6     4096    4096             0/0                 1
   4       7    16384   16360            24/24                2. (5)
```

The linear threshold is also recovered exactly.  `(M_0,M_1)` is injective
through length three and first fails at length four by

```text
{0,3} versus {1,2}.                                     (6)
```

This is the positional ambiguity already anticipated in THM-809/818.

## 3. Translation to staircase B2 layers

At a nonfixed THM-553 crossing clock `tau<n`, one representative from every
reflection pair is a tile `(a,b)` satisfying

```text
a+b-1=tau,       a>=b+2.
```

There are exactly

```text
s_tau=floor((tau-1)/2)                                  (7)
```

such positions.  Each is labelled by one of the four ordered bit-pair states
`00,01,10,11`.  The largest nonfixed layer at size `n` therefore has

```text
s_nonfixed(n)=floor((n-2)/2).                            (8)
```

The reflection-fixed clock `tau=n` contains
`floor((n-1)/2)` binary positions.  One is the canonical apex `(n,1)`, whose
bit is zero on the chosen complement-line endpoint, leaving

```text
s_fixed_free(n)=floor((n-1)/2)-1.                       (9)
```

Both (8) and (9) are at most six exactly through `n=15`.  At `n=16`, the
new nonfixed `tau=15` layer has seven positions and admits (3).  Thus

```text
S2 + per-state M1 + per-state M2
```

is literal-exact through `n=15`, and `n=16` is its first layer-level failure.
The verifier constructs the staircase tiles directly for `5<=n<=17` and
checks (7)--(9), rather than assuming the formulas.

At `n=9`, all nonfixed layers have at most three positions and the fixed
layer has only three free positions.  Thus `S2+M1` is already literal-exact,
as THM-818 states.  At `n=10`, (6) can occur, but the quadratic moment repairs
it.  The quadratic carrier then remains guaranteed through the full `n=14`
target and one size beyond.

## 4. Consequences for the relation-join frontier

THM-818's join can have two mathematically different outcomes.

1. If raw `S2` removes every joined pair, the lower codec is already
   injective at `n=9`, and the moments serve only as a proof-independent
   checksum.
2. If raw survivors remain, the first positional moment must separate every
   pair of distinct tilings.  The joined collision genealogy then measures
   exactly how far the economical count quotient is from the literal carrier;
   it is not an obstruction to finite exactness.

For subsequent sizes the same decision tree uses `(M_1,M_2)`.  This turns
“add more positional information” into a bounded exact statement through the
LRC(14) range.  It does not make the full `2^((n-1 choose 2))` tiling universe
cheap: a literal-equivalent static sidecar may still have too many states,
and THM-796's non-lumpability problem remains.

## 5. Tournament Analysis and preservation boundary

The verifier takes `M0`, `M0+M1`, `M0+M1+M2`, and the literal length-seven
word as information-carrier vertices.  Its pairwise observable is the number
of unordered four-state word pairs separated relative to `M0`; retention and
retention per integer field are the two switches, with one edge flip between
their carrier tournaments.  This is planning telemetry, not part of the
subset proof.

The assumption challenged here is that an aggregate curvature or class label
should repair aggregate layer counts.  THM-814 already shows the opposite:
its sixteen black-orbit failures have identical curvature and differ only in
position.  The correct repair is a positional address, and (1) quantifies the
minimum universally sufficient moment depth over the entire relevant size
range.

The carrier preserves every literal bit in each mirror layer through `n=15`,
and hence the apex-oriented tiling.  It deliberately does not preserve a
tournament isomorphism, merged-node ordering, compatibility between distinct
face relations, owner-labelled LRC state, metric gaps, wall chronology,
continued-fraction phase, or future lift/deletion behaviour.  Vertices here
are layer positions and states—not runners or unmarked tournament nodes. ∎
