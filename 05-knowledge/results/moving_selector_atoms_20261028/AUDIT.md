# Audit and reproduction

2026-10-28. Written proofs plus exact computations; no formal build or
independent-agent audit. Local C1--C6 refer to THEORY.md. Sources are pinned
by normalized text hashes, and the inherited computations are credited.

## Proof gates

C1 declares a read-only vertex selector, the rooted observation, cursor
transport under insertion and reset, and the retained index. It proves
closure and necessity, without claiming the quotient is discrete at all
sizes. The nine-pair is the previous fixed-boundary obstruction.

C2 changes the input to ordered path covers inside children. All singleton
runs are legal; individual cursor responses use possibly different
automorphisms. The assembled Hamilton paths themselves are legal on the
same source, and the three children and resulting partners are strong in
the primary example. The targets differ by an ordinary invariant H. The
five-child minimality concerns two distinct fixed outer ports and arbitrary
complete run orders; thirteen is not asserted to be the least total size.
One bit is complete only on the displayed two-class fiber.

C3 retains a simultaneous source-target alignment. Agreement recovers
directed paths; an automorphism fixing each such path setwise is identity.
The action on complete component orders is free. Oddness makes all these
component permutations even, but the choice of parity zero requires a
reference. The base criterion quantifies over every complete ordering,
including the transposition hostile. Faithful classification is distinct
from the unresolved behavioral minimum under supplied-run quotient-walk
contexts. Incoming persistent C3 tags settle the richer selected-substitution
menu on Hamilton presentations; CONNECTIONS.md audits the all-size proof.

C4's base lower bound uses independently rotating bottom triangles; its
upper bound uses the full wreath-group theorem and one vertex per triangle.
F64/F256 are orientation-permutation guards here, not new field census runs.
Paley27 is constructed over an explicit irreducible cubic. The displayed
Frobenius subgroup is not assumed to be its full automorphism group.
The wreath27 reversal count and the distinct Paley27 example are not identified.

C5 minimizes equitable state partitions preserving level and final
acceptance. It does not minimize arbitrary arithmetic programs or require
the intermediate predicate v->pole to be observed. Guarded variants add
exactly the stated current-vertex observations. Cell initialization is the
total number of labelled initial states; nonfree cells forbid naive division.

C6 allows every induced oriented spanning tree. Its generators are always
available in that language; this is precisely the assumption unavailable
in the legal-Hamilton metagraph. The rank proof includes n=2 separately.
Unlabelled connectivity uses an odd relabeling, not a claimed labelled
translation between the two parity sheets.

## Finite universes and independent paths

* probes.py: padded controls at7,9,11,13,15 with every cursor position on
  both complement sheets; all39 presentations of wreath9; the exact13-
  vertex two-path examples. Hamilton counts use independently written
  forward and backward subset recurrences. The transitive-separator example
  also uses a full canonicalizer. No large random sample is called exhaustive.
* frames.py: all593 directed spanning linear forests on the complete strong
  atlas through5, with every component permutation and parity sheet; full
  minimum-base searches on the five named atoms; all six middle orders of
  the five-child fiber; actual81/27 groups, attachment orbits and transporters;
  a concrete nonidentity even holonomy; all field elements and Frobenius
  arcs of Paley27. The large wreath base theorem is proved, not inferred
  by exhaustive groups at27 or81.
* folded_minimum.py: all531 normalized lifts of core size0..3 and wreath9.
  Iterated partition refinement and an independently implemented bottom-up
  signature calculation yield identical partitions. Every intermediate
  cell mass equals the sum of a separate labelled-state recurrence; totals agree
  with the inherited labelled DP. The complete16-cell quotient and both
  selected refinements are exported, including states and transition weights.
* tree_actions.py: every18248 labelled tree at orders2..7, and all1050
  tournament/tree pairs through4 for operation conjugacy. Gaussian
  elimination checks the proved ranks; this is not a Hamilton-legality census.
* Incoming check_strong_tags.py: separately replayed118 contexts and944
  active checks in normal and optimized Python. This is the incoming
  implementation, not an independent algorithm. Its universe is exhaustive
  through order4, with four chosen presentations and two ordered pairs each
  at order5. The all-size strong-module argument is audited in CONNECTIONS.md.

The atlas, field primitives, canonicalizer and DP routines have disclosed
shared ancestry. Replaying an inherited routine is not described as an
independent implementation. All checks raise exceptions under python -O.

## Commands

From the repository root, standard-library Python3 suffices:

```text
python -B 05-knowledge/results/moving_selector_atoms_20261028/verify.py
python -B -O 05-knowledge/results/moving_selector_atoms_20261028/verify.py
python -B 05-knowledge/results/joint_port_germs_20261028/check_strong_tags.py
python -B -O 05-knowledge/results/joint_port_germs_20261028/check_strong_tags.py
```

The verifier checks source/artifact seals, reruns every producer, compares
all saved receipts, and rechecks the seals. Individual producers overwrite
their own receipt but never reseal it. seal.py is an explicit maintainer
step, never invoked by verify.py. Text hashes normalize CRLF to LF; gzip
hashes retain original bytes. Optimization reports are outside the sealed
manifest because their optimization flag differs.

The previous16/24/64 counts have different scope from the three base sizes
or the component-order fiber counts. No aggregate coincidence is used as
a proof. No A000568 generation job, large self-converse census, or full
general-composition minimality theorem is claimed.
