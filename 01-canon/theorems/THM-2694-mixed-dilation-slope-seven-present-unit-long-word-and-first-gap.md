---
id: THM-2694
title: "Mixed dilation/slope-seven present-unit long word and first gap"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  One strict
  THM-2680 dilation component refines at its following endpoint to a fixed
  THM-2640 rail/sector/edge packet.  Its canonical slope-seven orbit contains
  95 private vertices on the same initial open cylinder, joined by 94 lawful
  slope steps with grammar (1^6,3,1^4)^8 1^6.  Including the initial
  pre-dilation vertex gives a heterogeneous 96-vertex, 95-edge word retaining
  the literal rail, present factor, delayed Boolean word, predecessor carry,
  future half-digit, private root, and primitive-unit flag.  Every one of the
  twelve canonical next representatives +7*delta/13^6 misses the present
  factor on the whole cylinder; the first root/unit-admissible continuation
  is n=113 and also fails there.  The eleven-state residue grammar cycles
  eight times before this lifted failure, and n=0,91 give a coarse C91 return
  with identical retained metadata but physical drift 49/13^5.  This is a
  finite fixed-configuration word and pumping no-go, not a configuration-
  switch theorem, endpoint current, infinite cycle, row exclusion, or
  LRC(14) conclusion.
source: root/mixed-grammar-long-word-scout-2026-07-28
depends_on:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
related:
  - THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
  - THM-2689-affine-clock-support-typing-tradeoff-and-odometer-phase-locality
  - THM-2691-lawful-odometer-alternating-rail-horizons-and-depth-six-present-collapse
script: 04-computation/lrc14_mixed_slope_long_word_probe.py
output: 05-knowledge/results/lrc14_mixed_slope_long_word_probe.out
script_sha256: 3525992454f488cd4c640e21fe8c0677aaa6e623b5375c331f19e8ab30fc1e0e
output_sha256: 2e426f74f1dc92c0b3282bf8175a6b23e1730a8bca24b8d7c7271a6711dfd1ee
secondary_script: 04-computation/lrc14_mixed_slope_long_word_independent_referee.py
secondary_output: 05-knowledge/results/lrc14_mixed_slope_long_word_independent_referee.out
secondary_script_sha256: a79bcd4a1532780dcd72b7528d2006de8de20f144f80bf327c406a8c2922f758
secondary_output_sha256: 1d261bc9646c1e0bc789ca97b08a1b1665f0550af931d3ced2367fef2df19e3e
hash_basis: LF-normalized bytes
---

# THM-2694 -- one mixed D/slope word survives 95 edges, then present support fails

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2680 proves that dilation can join two physically typed positive atoms,
but its three-event chronology is nilpotent.  THM-2640 provides the formal
slope-seven carry/root covariance, but warns that the present packet is not
translation-invariant.  The point of this theorem is to compose the two
operations before discarding either sidecar.  A surprisingly long finite
word survives, and its first failure identifies the missing coordinate
exactly.

## 1. Frozen physical D edge

Let

```text
R=13^6=4826809,
x=649039434905733/1304692766858936,
z=D(x)={13x}=46873542509301/100360982066072.              (1)
```

Use the following exact THM-2680 atoms:

```text
current:   (source,rail,j,h,epsilon,kappa)=(1,2,5,2,0,0),
following: (source,rail,j,h,epsilon,kappa)=(1,0,2,6,1,1). (2)
```

Their delayed labels glue because `h_current=j_following=2`.  The component
of their literal D-fibre product containing `x` is the strict open interval

```text
I=(960117507257/1930018885886,
   324519717452867/652346383429468),                      (3)
```

of length

```text
|I|=1/652346383429468.                                   (4)
```

The midpoint of `(3)` is exactly `x`.  The independent referee reconstructs
both endpoints by intersecting the current base atom, following pulled-back
base atom, current delayed prefix, and following pulled-back delayed prefix.
The following delayed factor is binding on both sides; this is not a sampled
or rounded interval.

## 2. Fixed following configuration and slope covariance

Refine the following endpoint to the THM-2640 configuration

```text
(rail,sector,edge,kappa,h,shallow,carry,root)
                    =(0,0,0,1,6,1,2,6),                 (5)
```

where rail `0` has the literal key `(1,0,12)`.  For an integer address `n`,
put

```text
z_n=z+7n/R.                                               (6)
```

No wrap occurs on the finite word below.  More generally all support tests
in the complete orbit are made modulo one.  Translation by

```text
tau_delta=7 delta/R,             delta in {1,...,12},     (7)
```

is the canonical displayed THM-2640/2657 slope-seven representative.  It
preserves the delayed coordinate because `R*tau_delta=7 delta` is integral,
and it changes the two private digits by

```text
carry -> carry+7 delta,          root -> root+delta        (mod 13). (8)
```

For the fixed coefficient row in `(5)`, exact finite-field evaluation gives
the private root/unit residue bank

```text
A={0,1,2,3,4,5,6,9,10,11,12}.                            (9)
```

Thus residues `7,8` are the two local algebraic holes.  Equation `(9)` alone
does not address the moving rail or present supports; those are tested
physically below.

## 3. Complete orbit census and the common 95-vertex branch

Scan all `n in Z/RZ`, retaining `(9)` together with the literal rail and
present factors at `z_n`.  There are exactly

```text
3346                                                        (10)
```

retained addresses.  Sort them cyclically and join successive points only
when their positive forward gap is at most twelve.  The gap histogram is

```text
1^3016, 3^295, 45^3, 46^12, 47^3, 48^3,
3564^4, 3565^1, 3566^1, 3567^1,
685282^1, 685283^5, 685285^1.                            (11)
```

There are `35` components.  The largest has `229` midpoint addresses and
runs from `1379246` to `1379514`.  This complete-orbit census is a midpoint
support classification only; it is not promoted to a 229-vertex common
cylinder.

The circular component containing address zero has size `104`.  Its positive
forward branch from the frozen D endpoint is exactly

```text
P={n: 0<=n<=110 and n mod 13 in A}.                       (12)
```

Hence `P` has `95` private vertices.  Its `94` slope edges have spelling

```text
(1^6,3,1^4)^8 1^6,                                      (13)
```

namely `86` steps of size one and `8` steps of size three.  Prepending the
single D edge from `x` to `z_0` gives

```text
96 vertices,                 95 edges.                   (14)
```

The distinction in `(14)` is load-bearing: there are 95 private slope-orbit
vertices, not 95 total vertices.

## 4. One common open cylinder, not 95 unrelated points

For every `x' in I`, transport the following state by

```text
z_n(x')={13x'}+7n/R.                                     (15)
```

Exact interval containment proves that for every `n in P`, the same open
interval `I` lies in all of the following pulled-back factors:

```text
literal rail,
current-present factor,
delayed Boolean prefix,
predecessor carry cell,
future half-digit,
private deep half-tooth,
primitive coefficient-unit row.                         (16)
```

The common symmetric radius about `x` is

```text
1/1304692766858936,                                      (17)
```

so `(16)` has precisely the positive length `(4)`.  The primary companion
computes every exact factor margin and then rechecks containment as an open
set.  The independent referee instead reconstructs `(3)` from four affine
inequalities and tests both endpoints against every factor at all 95 private
vertices.

## 5. The first gap is the present coordinate

The terminal private address in `(12)` is `n=110`, with residue `6`.  Apply
each of the twelve canonical next representatives `(7)`.  They land at

```text
n=111,112,...,122.                                       (18)
```

For every address in `(18)`, the pulled-forward open cylinder still meets the
literal rail, but is disjoint from the present factor.  Addresses `111` and
`112` additionally fail the root/unit gate.  Address

```text
n=113,                 residue 9                          (19)
```

is the first root/unit-admissible continuation, and its entire cylinder is
still present-free.  Therefore `(13)` is maximal for this inherited D
component, this fixed configuration, and the twelve displayed lifts.

This statement does **not** range over every THM-2657 lift `k/R` satisfying
`k=7 delta mod 13`; those lifts may differ by multiples of thirteen and move
the low-speed factors differently.  It also does not rule out switching to a
different rail, sector, edge, present clock, or D component.

## 6. Cyclic quotient, finite lift

Modulo thirteen, `(13)` projects to the genuine directed cycle

```text
0->1->2->3->4->5->6->9->10->11->12->0,                  (20)
```

with edge steps `(1,1,1,1,1,1,3,1,1,1,1)`.  The lifted path completes this
cycle eight times and then reaches the prefix `0->...->6` at odometer height

```text
floor(n/13)=8.                                           (21)
```

The next projected edge `6->9` is precisely the failed lift `110->113` from
`(19)`.  Thus a cycle in the residue transition graph cannot be pumped after
forgetting the integer height.  This is the same controlled-forgetting
phenomenon that THM-2646 exposes for modular braids versus central height,
but here it is an exact LRC carrier witness rather than a knot invariant.

There is also a sharp graph-theoretic form.  Give the selected chronology
only its 94 displayed consecutive edges.  Its lifted graph is the path
`P_95`, hence a partial cube.  Quotienting its vertices by `n mod 13` and
forgetting multiplicities gives exactly the cycle `(20)`, namely `C_11`.
That quotient is odd and therefore not bipartite; since every partial cube is
bipartite, it is not a partial cube.  Thus forgetting odometer height can
destroy partial-cube structure even while preserving a perfectly coherent
cyclic word.  This is complementary to THM-2606's partial-cube/origin lesson.
The claim concerns the selected chronology graph, not the larger graph with
every admissible slope chord installed.

Gracefulness does not detect this loss.  The standard alternating labelling
`0,94,1,93,...,47` makes `P_95` graceful.  The cyclic labelling

```text
(0,4,6,5,8,3,9,2,10,1,11)                              (21a)
```

makes `C_11` graceful because its cyclic edge differences are exactly
`(4,2,1,3,5,6,7,8,9,10,11)`.  Both certificates are checked internally by
both companions.  Hence gracefulness survives this height-forgetting
quotient, while bipartiteness and partial-cube structure do not; gracefulness
alone cannot certify chronological liftability.

There is an even sharper coarse return inside the positive branch.  Addresses
`0` and `91` agree modulo both seven and thirteen.  Direct reconstruction
shows that their fixed configuration, carry, root, unit flag, delayed word,
present flag, rail weight, and inherited atom weight are identical.  Yet the
physical points differ by

```text
z_91-z_0=7*91/R=49/13^5.                                 (22)
```

An integer-frequency endpoint character has the same phase at these two
points exactly when `13^5` divides its frequency.  Hence the coarse `C91`
return is not an endpoint/current return.  The destroyed coordinate is the
exact odometer address (already visible through `floor(n/13)`), and an
endpoint-frequency or equivalent height sidecar is necessary.

## 7. Relation to the nearby boundaries

This theorem does not conflict with the existing no-go results:

* THM-2672 and THM-2687 concern simultaneous slope faces across missing
  labels and arbitrary configuration switching.  The word here repeats only
  the eleven labels in `(9)` and never supplies a thirteen-label face.
* THM-2689 restores intrinsic affine support but loses a global seven-clock
  action.  No such action is asserted here.
* THM-2691 follows the alternating `+/-(13^5+1)` rail germ and finds its
  dynamically typed present failure at horizon six.  The present theorem
  instead starts with one D edge and then uses monotone small slope-seven
  representatives; it is exactly one of the mixed faces left open there.

No semantic endpoint current, target action, terminal word, configuration-
switch continuation, infinite cycle, all-row statement, row exclusion, or
proof of LRC(14) follows.

## 8. Exact reproduction

Run

```bash
python3 04-computation/lrc14_mixed_slope_long_word_probe.py
python3 -O 04-computation/lrc14_mixed_slope_long_word_probe.py
python3 04-computation/lrc14_mixed_slope_long_word_independent_referee.py
python3 -O 04-computation/lrc14_mixed_slope_long_word_independent_referee.py
```

Normal and optimized executions of both companions match their stored
outputs byte-for-byte.  The primary script performs the full `13^6` orbit
scan and common-cylinder verification.  The secondary script does not import
the primary: it reconstructs the D interval from local affine bounds and
checks the 95-vertex word and twelve terminal failures directly.  A separate
hostile audit independently enumerated the complete orbit and reproduced
`(10)`--`(11)`, all 35 component sizes, the largest component, and maximum
gap `685285`.
