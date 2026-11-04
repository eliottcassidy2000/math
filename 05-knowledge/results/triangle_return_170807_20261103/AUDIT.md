# Scope, proof audit and reproduction

2026-11-03. Same-agent audit with independent exact algorithms. No separate
reviewer or Lean formalization is claimed. A computation establishes only
its declared finite universe; all-size extensions use the supplied proofs.

## Proof obligations

* **Native interface, FINITE-EXACT:** source names and full rows are pinned
  to return18979_palindrome/returns.json. Both helper programs are replayed
  as actual path complements. Only the helper is deleted; the changed
  old-edge set is exactly the declared singleton. Root50/spare6 frames
  agree between each pair of paths. Completed repairs and target are strong.
  H is independently56935,68533,67599,170807.
* **R3--R5 in RETURN, PROVED with exact coefficients:** at most two visits
  to either module follows from its actual outgoing incidence and the
  single exit31, for arbitrary nonempty module tournaments. Four TT1/TT2
  values are computed by both forward and backward subset DPs. The full
  sixteen TT1,TT2,TT3,C3 input pairs agree. Positivity C2>=(n-1)H handles
  every larger input, so R5 is not inferred from samples.
* **TR1--TR4 in TRIANGLE, PROVED:** uses the inherited maximal-run bijection,
  free action on ordered covers, rooted indegree defect, and stronger
  positive-cut bound w_uv>=H(T-{u,v}). Endpoints must be fixed individually;
  arbitrary weights need not be invariant. The all-mark result excludes
  m=1. The birth obstruction requires a proper triangle module and a
  strong original source. Monotonicity at the first triangle-splitting
  birth extends the obstruction to every finite induced extension, without
  requiring strongness at subsequent prefixes.102 labelled strong quotient/module controls at quotient
  orders3,4 pass. A whole C3 has+2 births, testing the proper-module guard.
* **PC1--PC4 in PERFECT, PROVED:** the inverse retains positive deficiency,
  prime q, coprimality, integral prime-power ratio and Euler shape when
  applicable. The q interval is half-open; an integral A/D forces the
  correct endpoint.2000 finite cores are independently checked against
  direct prime/exponent loops. The C0-family proof bounds all q>=17 by
  p<=53 and treats q=13 separately; it has no unbounded-exponent sampling
  assumption. The Descartes candidate fails primality and has the exact
  true-sigma excess28856318400.

## Complete menus and positive controls

All searches start with the same named16-vertex H56935 source. A birth mask
records **old vertices pointing into the newborn**. Complementing this mask
gives the out-neighbor cut used in the cut-theorem documents.

| Receipt | Finite universe and conclusion |
|---|---|
| one_edge_search.json |0/1 old-edge flips, every mark, then the fixed reset; also U2/U3/U4 and both offsets; no hits |
| two_edge_search.json |All7140 pairs of distinct old edges, every mark, fixed reset; no hits |
| two_births_all.json |All65536 first births,46 retained below the target, every second birth; no hits |
| three_births.json |All192 retained second-birth prefixes, every third birth; no strong or weak target |
| edge_birth.json |All120 old-edge flips, every birth; no target |
| edge_two_births.json |All120 old-edge flips,6583 retained first-birth prefixes, every second birth; no strong target |
| two_edges_birth.json |Discovery stopped at first strong hit, pair6771; unvisited pairs are not classified |
| positive_birth_control.json |Same unrestricted two-birth code at target170809 finds actual strong hits |

Positive-cut monotonicity justifies pruning only a sequence of births. All
192 retained third-stage input counts are independently replayed by the
backward DP. The local search does not identify iso-classes or quotient
old labels. It cannot exclude a different56935 source or a later deletion.
Root re-marking does not change H; no numerical phase match is used as a
legal transition.

## Arithmetic and phase safety

Discovery/inspection use n<=18; all ordered two-path covers are bounded by
(n+1)n!<=19!, and every19-vertex Hamilton count is at most19!, within signed
64-bit range. The independent backward counter allows n<=20 and H<=20!,
also within range. The endpoint checker has explicit overflow checks.
Python handles the number-theory and field arithmetic with exact integers.

For every completed interface state, audit.json stores the full relation,
ordered pair responses on all three physical ports50,6,51, and the current
golden triple. It recomputes all five Q4 presentations and both field charts.
The F343 vector codec is not the older interval spectral clock. Earlier
ancestral controls are pinned in return6325_completion/transport.json;
their source is not relabeled into the new target. The81/27 literal pair
and its actual cut response are replayed independently.

The module and arithmetic proofs concern actual operations. Neither a
matched count nor an observer quotient is substituted for a source.

## Run commands

```sh
python -B -X utf8 05-knowledge/results/triangle_return_170807_20261103/reproduce.py
python -B -X utf8 05-knowledge/results/triangle_return_170807_20261103/reproduce.py --full
```

The first runs independent controls normally and under-O, comparing exact
receipts. The second also recompiles and replays every discovery menu.
Source hashes normalize CRLF to LF. reproduction.json records whether the
full mode was run; a default rerun legitimately changes that field.

Inherited results, incoming work and literature provenance are separated
in RECOVERY. Global spectrum completeness, a reusable ternary source
grammar, and odd-perfect existence remain OPEN.
