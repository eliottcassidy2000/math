# The runner world has a tiling model: steps are the path, balanced relations are the tiles

*mac-mini-2026-07-07-S43 (HYP-4887, THM-639). Owner: "a tree on 8 events has 7 edges and
that is just the Hamiltonian path connecting each element in an 8 player tournament. the
tiling model leverages viewing a tournament from one of its Hamiltonian paths to reveal
the symmetry in the isomorphism class structure forced by the nature of intersubjective
binary relation itself, apply similar analysis to lonely runners." This note records what
that analysis yields — three small theorems, one refuted route, and a structural
correspondence I did not expect.*

## The decoding

A spanning tree on the event set is a FRAME: it splits the C(n,2) pairwise data into
n−1 tree coordinates (the frame itself, carrying no invariant content) and the rest,
which the tree turns into CYCLES (each non-tree pair closes a unique loop). In the
tournament world the frame is Rédei's Hamiltonian path, the tree data is the score
hierarchy, and the cycle data is the tiles — the binary tiling whose S_n-symmetry
IS the iso-class structure. The owner's instruction: find the same split for runners.

It exists, and it is exact (THM-639):

- **The frame** = the sorted order of the 13 speeds; its 12 tree coordinates are the
  steps `s_i = e_{i+1} − e_i`. Every pairwise "relative runner" `e_j − e_i` is an
  interval sum of steps — the tree reconstructs all pairs, exactly as the path does.
- **The tiles** = the balanced lattice `L_0(E) = {m : Σm_ie_i = 0, Σm_i = 0}`: the
  CLOSED LOOPS of the difference structure. Its two defining constraints are precisely
  the invariances forced by the relation itself — translation (Σm = 0: only relative
  positions are real) and homogeneity (Σm_ie_i = 0: only relative rates are real).
- **The law** (THM-639-C): the entire gap-process law of the family is Haar measure on
  the annihilator loop of the relation lattice; it factors through `L_0(E)`, and `L_0`
  is a COMPLETE invariant — lattice class = rational-affine class = primitive step
  sequence up to reversal. *The isomorphism classes of runner families are the step
  sequences up to reversal*, as the tournament classes are the tilings up to the fibration.

## The correspondences that turned out to be theorems

1. **AP = transitive tournament** (THM-639-B). The wall count — the number of cyclic-
   order changes ("tile flips") the configuration undergoes as x sweeps a period —
   equals `Σ_r s_r·r(k−r)`, uniquely minimized by the all-steps-equal family. The AP is
   the minimal-complexity frame of its world exactly as the transitive tournament is of
   its. Its three-distance witness structure, its Farey-cell tractability, and its role
   as the master extremal object are all faces of this minimality.

2. **Step reversal = transpose/complement symmetry** (THM-639-A). Gap laws are invariant
   under `E ↦ −E` = reversing the step sequence. Unique extremals must be palindromes —
   and the census delivers: the AP, death-star's family, monad's record, the S41
   PZ-minimizer are ALL palindromic step sequences. The tiling model's transpose-self
   (SC) classes have runner siblings, and they are exactly the record families.

3. **Spontaneous symmetry breaking = the SC/NS split.** The E[U]-minimizer is NOT a
   palindrome — so (by A) it is not unique: it comes as a mirror pair, and the
   palindrome-constrained floor is strictly higher (0.0988 vs 0.0938). Different gap
   functionals resolve the reversal symmetry differently: μ's extremal is
   self-dual (an "SC class"), E[U]'s breaks into a conjugate pair (an "NS pair").
   This is the runner-world shadow of the metagraph's spine/sea structure — the
   division of extremals into self-complementary and paired classes — appearing not in
   a finite graph but in a variational landscape.

## What the frame does for the open program

- It UNIFIES the week's scattered invariance arguments (dilation, translation, boxeph's
  speed↔co-offset reflection, "pairs never contribute", kps-S59's difference-mode
  uniformity) into one statement with one proof (Haar on the loop).
- It gives cheap PRUNING and CHECKS: palindromicity halves symmetric adversarial
  searches; any claimed unique extremal that is not a palindrome is wrong.
- It NAMES the enemy precisely: all deviation from the iid gap law is a function of the
  balanced lattice graded by vector norm — "how many short cycles at each scale". The
  c-adic cascade records are the families that tune their short-cycle spectrum.
- It does NOT (yet) close the sparse lane: the single-parameter girth bound is refuted
  (THM-639-D — 13 integers are pigeonhole-forced to girth 4, and the rank-11 packing
  count is astronomically generous). The honest repair is the successive-minima profile
  of `L_0` — the runner analog of asking not "is the tournament transitive?" but "what
  is its full flip-distance spectrum?".

## The reflection proper

Both worlds run the same program: an intersubjective binary relation, a canonical
spanning frame guaranteed by the relation's own combinatorics (Rédei's path / the total
order of ℤ), and the discovery that everything invariant lives in the loops the frame
closes. The tournament project spent hundreds of sessions learning that the metagraph's
structure is the cycle space seen through S_n; the runner project keeps re-learning that
its floors are functions of the balanced lattice seen through Haar. The two are one
lesson: **the tree is the observer; the cycles are the observed.** What made the
tournament side tractable was refusing to study tournaments one at a time and studying
the tiling space's symmetry instead; the corresponding move here — studying the space of
step sequences with its reversal involution and its cycle-spectrum stratification,
rather than one family at a time — is where the sparse lane's floor should eventually
come from.

## Pointers
THM-639 (proofs); `lrc14_hampath_frame_macmini_S43.py/.out` (census, constrained
descent, T(W) table); S41 HYP-4837 (the balanced lattice's first appearance + the
truncation dead ends); klein-S155 (the successor-digraph = the tournament ON the orbit —
the same dictionary read at fixed x instead of averaged).
