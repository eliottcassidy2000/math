# Independent review of reader acceptance and actual native repairs

**PASS: proof review, exact independent replay, and hostile-certificate
tests.** Normal and optimized runs produce identical review receipts. No
mathematical correction was required in
[`reader/PROOF.md`](../reader/PROOF.md) or
[`native/reader_repairs/PROOF.md`](../native/reader_repairs/PROOF.md).
All previously pinned files remain untouched.

## Quantifiers and actual objects

The finite acceptance claim is exactly32 targets: both affine families
`63*2^t-35` and `63*2^t-29` at t=1,...,16 in the inherited219-seed reader
bank. Twenty targets, t<=10, have materialized named graphs and native
source-to-target programs. The twelve tail cases have accepted actual-seed
reader descriptions, but no native programs or physical source-name
alignments are claimed. Neither range establishes all-t acceptance,
essential infinite switching, or recursive closure of both next digits.

Each source for the twenty programs is independently reconstructed from
the original five-vertex B2/B4 sibling suffix, its specified port set,
prefix names10,...,t+9, and outer root99. The saved source matrix is required
to equal this graph, not merely have the same H. Each target matrix is
independently reconstructed from its actual seed, retained mark, chronological
reader word, prefix clock, and final root. The independently selected target
order equals the minimum declared within the219-key bank. This is not a
global tournament-order or native-cost minimum.

The name-preserving claim is also exact: the verifier treats physical names
as graph vertices throughout. Original names are never renamed or deleted.
Some original prefix names change structural roles and incidences. There is
no claim that the original homogeneous prefix survives, that every original
arc remains, or that equal-count graphs can be exchanged without a program.
Thus these positive finite certificates do not contradict the fixed
homogeneous-clone exclusion.

## Independent finite reader inverse

The review imports no producer code. First it enumerates literal paths on
every subset of all219 small actual seeds and independently reconstructs
all eight stored seed fields. All seeds are checked strong.

For each fixed target and seed, the prefix clock obeys
`h+2^s F<=N`. At each such clock, write

```
D=N-h-2^s F, K=u-q+2^s p, d=q-e.
```

The inherited identity is `D=K*r+q*A(r)+d*(r mod2)`. Since
`0<=A(r)<=r`, every positive address lies in the complete finite interval

```
max(1, ceil((D-d)/(K+q))) <= r <= floor(D/K).
```

The reviewer enumerates this interval and evaluates A from the chronological
binary reader word, independently of the producer's recursive section
inverse. Address0 is checked separately. Every nonempty word ends0 because
the highest address bit is1; this is the correct strong-ending reader
language for the strong marked seeds. The endpoint reader recurrence and
the actual target construction agree on every saved hit.

Across all32 targets, the independent interval audit examines1,758,518
candidate integers and recovers precisely the176 saved seed/word/clock
certificates. It verifies every reported acceptance count and minimum
bank order. The finite common-(seed,clock) intersections on the four
declared intervals ending at t10 are empty, as claimed. These observations
do not determine an eventual selection rule.

## Native interpretation and the bridge proof

The review represents each named graph by a dictionary of adjacency sets.
A native move first validates that its guard lists each current name exactly
once and follows actual directed arcs. It then reverses every arc and
restores the consecutive guard arcs. This differs from the producer's
pairwise nonpath-edge interpreter.

Each birth checks a fresh name and an actual incoming cut containing only
current vertices. Each two-move arc gadget is compared with a literal
reversal of its specified original pair. Helpers are deleted by physical
name; the exact deletion list must be the reverse of that gadget's fresh
birth list. No original or permanent name can be deleted. Every final graph
must equal both the saved endpoint and the declared target as a named
orientation, independently of matrix index order.

All20 programs pass with the exact reported costs:

```
470 arc gadgets; 940 native complements;
13 permanent births;
491 temporary births and491 temporary deletions;
21 two-helper bridge gadgets; maximum peak order19.
```

The all-order bridge formula is valid for every actual arc a->b, including
the two-vertex boundary where the remaining path P is empty. The inserted
bridge z makes `(a,b,z,P)` directed and ensures b is not the last vertex.
After adding h with sole old outneighbor a, the first guard is
`(a,b,z,P,h)`. Its native complement directs a->h and h->b; the other
consecutive edges of `(a,h,b,z,P)` were retained. The symmetric difference
of the two guard edge sets contains exactly the old pair{a,b} among old
vertices. Deleting h,z therefore effects precisely that reversal. This
proof supplies the fallback independently of bounded guard-search success.

Neither the proof nor the data require intermediate graphs to preserve H,
strongness, phase, root stars, or common frames. Only the original sources
and final targets are checked strong. These certificates do not prove
fixed-order native or black-core connectivity.

Four corrupted copies of a real certificate are deliberately rejected:

| Corruption | Rejection |
|---|---|
|Reverse a saved guard word|guard is no longer directed|
|Reuse an original name for a helper|birth name is not fresh|
|Delete an original vertex instead of a helper|deletion list is invalid|
|Reverse a pair in the declared target|final named graph differs|

## Fresh endpoint, joint, and phase checks

`reader_profile.cpp` uses forward path extension on the actual graph and
its converse. The producer's full profiler uses last-vertex deletion
recurrences, and its native packet uses a separate feedback-gap engine.
The review recomputes all40 actual endpoint profiles independently, including
H, P2, every E1/S1/E2/S2 coordinate, and every ordered J pair. The largest
endpoint graph has17 vertices; all two-slot counts are bounded by18!, so
signed64-bit arithmetic is sufficient.

Every full profile equals the reader's stored profile. The six native
port charts and six ordered native J queries are checked by physical name
against those independently recomputed full profiles. The golden,
F49-plus-charge, and F343 coordinate formulas are then decoded afresh at
every endpoint vertex. This validates observations and the stated inherited
linear coordinate transformations; it does not construct a physical
Frobenius action or identify these graphs with the unrelated81/27 controls.

The all-vertex entry claims also pass independently for all20 targets:

```
S_v(2)-H is never1 or2;
E_v(2)-H is never1 or2;
E_v(2)-H=0 holds exactly at physical outer root99.
```

Consequently the stated immediate whole-graph C2 obstruction is sound:
`H_new=3H+2*(S_a(2)-H)` cannot yield either positive digit at any re-marking,
and the end-defect test gives the same obstruction after taking the
converse. This is a statement about the twenty selected actual targets and
that particular next operation. Other accepted witnesses, preparations,
larger macros, and t-dependent grammars remain open.

The additional logical corollary in [`CLONE_DECIDER.md`](../CLONE_DECIDER.md)
also passes review. For a specified pure clock `a*d^t+c`, the surjection
basis identity gives exactly `A1=ad+c` and `A_r=a*binom(d,r)` for2<=r<=d,
with higher coefficients zero. Positivity and exponential independence
justify the eventual-to-all-t implication; the strong-order bound makes
the remaining marked-base search finite. The necessary restrictions
`(d-1)! | a` and `|R|>=2d-1`, and the separate singleton exception, are
correct. This decides the stated fixed homogeneous-inflation class only.

## Reproduction

```
python 05-knowledge/results/c4_port_returns_20261103/audit/check_reader_review.py
python -O 05-knowledge/results/c4_port_returns_20261103/audit/check_reader_review.py
```

The code compiles its profiler in a temporary directory. Outputs
`reader_review.json` and `reader_review_optimized.json` are byte-identical.
The native producer's normal/optimized receipts are also checked for byte
equality. The audit pins the actual seeds, target and program certificates,
proofs, sidecars, and independent verifier by SHA256.
