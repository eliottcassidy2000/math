# LRC(14) z216 costly gcd-8 probe — lrc-probe-agent

**Status: FINITE-EXACT SCRATCH RESULT + INDEPENDENTLY AUDITED; NOT CANON.**

The two gcd-eight wall rows omitted only by the intrinsic work cutoff in
THM-3264 both close under the existing THM-3139 exact screen plus complete
one-high terminal bank.  A wrapper-free independent audit confirms the
closure.  If promoted, this removes
the entire gcd-eight stratum at `z1=216`, changing the current composed ledger
from `373186` to `373184` and the occupied layer from `382` to `380` wall rows.
It does not touch the gcd `24/36/72` families, arbitrary `k<=1`, the rung, a
physical-cover classification, or LRC(14).

## Inheritance pass

- Closest proved mechanism: THM-3264,
  `01-canon/theorems/THM-3264-projected-k3-z216-low-cost-gcd8-seventeen-row-terminal-descent.md`,
  which closes the other 17 gcd-eight rows with the THM-3139 upper-relaxation
  screen and high-slot terminal bank.
- Canonical hostile: zero-high relaxations are deliberately retained and then
  excluded only by the wall gate; duplicate-permitting two-high assignments
  must be killed by a strictly positive exact gap.
- Corrected near miss: MISTAKE-331/333.  Exact LP validity does not make the
  solver-selected Farkas basis canonical.  The probe hashes only `row[:19]`,
  direct/inherited branch counts, and witness-free terminal semantics.
- Least-used relevant sidecar: the pointed divisor quotient.  Every surviving
  mask is at `D/first_d=8`, the top of the four-point divisor chain; it is used
  as a typed residual coordinate, not promoted to a `C2` action.

Meta-patterns used: “Turn certificate failure into an address, then change
sidecars,” “Use redundant paths as detectors,” and “Exactify before
interpreting a decimal.”

## Session portfolio and concept board

| lane / object | predicate | inherited invariant / extremal | operation | lost coordinate / obstruction | cheapest next test |
|---|---|---|---|---|---|
| Anchor: costly gcd-8 rows 238/370 | projected row empty? | wall status, `L=560560`, `first_d=70070` | exact THM-3139 screen then terminal lift | screen alone leaves 115/19 masks | run full terminal bank |
| Niche: two residual mask banks | is one a quotient/subbank of the other? | common five-label body tail | set intersection | terminal semantics depend on the leading label | compute literal containment, do not assume symmetry |
| Sidecar: top `q8` divisor address | do residuals descend below the top? | `D/first_d=8` for every survivor | divisor stratification | no action, phase, owner, or common carrier | histogram every residual divisor |
| Remaining wall families | can the same mechanism close gcd 24/36/72? | THM-3281 empties three whole ruler families | rank families by `L*r`, screen complete families | 380 rows remain after candidate composition | batch the next full low-invoice ruler family |
| Wildcard: THM-3285 horn / THM-3277 backbone | can phase data supply wall semantics? | 63-label local `R-M-R`; 12 target-class geodesics | compare loss ledgers | both forget endpoint origin/current or physical action | no probe: first require a typed map into projected wall rows |

The phase wildcard did not overtake the anchor: THM-3285 retains literal
ancestry only on three address cylinders and THM-3277 is an internal critical
quotient of vertex-simple paths.  Neither currently maps to the two projected
wall rows while preserving the physical predicate.  No syntax-only bridge was
manufactured.

## Exact result

The universe is the complete maintained `z1=216` projected atlas of 480 rows.
The target rows are exactly

```text
238: E=(1,8,10,11,13,14), L=560560, r=34, high_floor=55207, wall=True,
370: E=(2,8,10,11,13,14), L=560560, r=34, high_floor=55207, wall=True.
```

They are the two rows excluded from THM-3264 solely because each has intrinsic
invoice `L*r=19,059,040 > 2,000,000`.

The first draft of this scratch note mislabeled the atlas field `high_floor`
as `components`.  The independent audit caught the tuple-position error.
The corrected component count is `r=34`; the invoice, screen fields, and
every downstream terminal calculation were already indexed correctly.

The exact screen gives

| row | states | crude | exact status | residual | residual quotient |
|---:|---:|---:|---:|---:|---|
| 238 | 1035 | 237 | 683 | 115 | `D/first_d=8` for all 115 |
| 370 | 354 | 213 | 122 | 19 | `D/first_d=8` for all 19 |

All 805 status exclusions are inherited full-table Farkas contradictions
checked over exact rationals.  The residual banks are genuinely asymmetric:
they are not equal, but the complete 19-mask bank of row 370 is exactly the
intersection with (hence a subset of) row 238's 115-mask bank.

For row 238, the duplicate-permitting two-high gap is

```text
66973458766897 / 29716657283313000 > 0.
```

The terminal enlargement retains 109 zero-high hostiles and enumerates 137
one-high cases over nine body-local low-label sets.  It closes 131 by coarse
cardinality and the remaining six by exact complete-cell cardinality, with no
max-gap fallback, no failure, and minimum support slack one.

For row 370, the corresponding gap is

```text
3668224417997 / 1297249202749500 > 0.
```

The terminal enlargement retains 16 zero-high hostiles and enumerates 19
one-high cases over two body-local low-label sets.  It closes 16 coarsely and
the remaining three exactly, again with no max-gap fallback, no failure, and
minimum support slack one.

### Mechanism and logical direction

The THM-3139 screen is a necessary upper relaxation of each projected row.
For each residual body, strict positivity of the duplicate-permitting
two-high gap rules out two high suffix slots.  The wall condition requires at
least one, so every actual survivor has exactly one high slot.  Replacing that
high band by its ray supremum enlarges the feasible set.  Every case in that
complete enlarged one-high bank violates the necessary projected support
cardinality.  Therefore both original projected rows are empty.  No converse
from quotient feasibility is used.

### Exact digests

```text
row 238 mask:      f7db818098c70097349459e472e3da4a70daf21627d5f8158b2dc4963ff7210a
row 238 screen:    ab213499b0b774678caa4b3658539f192e2b423e2b4832d28eb384e7b3884de7
row 238 terminal:  83dcd9a3c35c3b6a1da315881b3ba0cd469239b9f8007d19146fc628ead2c9bf
row 370 mask:      4c383892c25798129ad19799bee95996d78c3bb7506501ab9ea0f9f34662d05d
row 370 screen:    f724cac04852034b8d813e0b1c9259579494511e836ba4f5ab2cd7e8462de971
row 370 terminal:  428ea1129fcbd8c1e607a26b9322aee7f2a82553d6dee1603a5329887b0ce3b3
```

The probe pins the LF-normalized source/output/semantic triples of THM-3139
and THM-3264.  Normal and optimized two-process runs printed identical bytes.
A prior one-process run independently reproduced row 238's screen counts and
canonical data, but this is execution redundancy, not a second implementation.

## Dependency manifest and inherited-code trust

The audited baseline checkout was `c2cd7143e661d3b270ea0684cc9107ddd0d92691`
with Python `3.13.14` (64-bit Windows).  Direct status/provenance documents are

```text
01-canon/theorems/THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent.md
01-canon/theorems/THM-3264-projected-k3-z216-low-cost-gcd8-seventeen-row-terminal-descent.md
```

The scratch wrapper directly pins these source/output/semantic triples:

```text
THM-3139 source:
  04-computation/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py
  92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36
THM-3139 output:
  05-knowledge/results/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out
  4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5
THM-3139 semantic:
  1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a

THM-3264 source:
  04-computation/lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.py
  d9841a86850c6609e62d0522b50ea38722cd5850f7173374a3b143ef577e0e3e
THM-3264 output:
  05-knowledge/results/lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.out
  2eebe5f7acb2d1c02fa126ee03166cd274ed209cc52d4b8a729a8fbc5f0a9782
THM-3264 semantic:
  64604458a454f7da3468b07ec5697a6dd62ac4413625a92dbd8e8ffff15e1a7e
```

For an independent audit, freeze the full dynamic import chain as well:

```text
04-computation/lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.py
  8de1e3d03b5070a84b040ac13a173a598646107f85e7ba0defc2ca070808f162
04-computation/lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.py
  1e23ec19fa147c55fb6d38a965eedae0132f5e069b9f820bfd5c300dce4d8f89
04-computation/lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.py
  42323171481deba2371eed9947b2079976cb367dac340cf58b8f1f0c0afb5082
04-computation/lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.py
  1f74f2b2368c04f514f2c388b54c70a9ee66c9387fbc437093884b807b3eb23c
04-computation/lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.py
  f6f64ab8d8ea9b04a1a03e26fc6026efc864e44518e9cb40df4fe8471a4a7991
04-computation/lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.py
  2a051babe109f56056fe61476870f8e2e13cfc99b2f9bb7ac122b8780c8fa168
05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out
  cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda
```

The current scratch wrapper does **not** independently rederive or individually
hash-gate that transitive chain.  It trusts the promoted THM-3139 implementation
of `screen_worker`, the inherited THM-3078 implementation of `terminal_probe`,
and their theorem-level audits.  In particular:

- the proof that the screen and high-slot terminal objects are genuine upper
  relaxations is inherited, not reconstructed here;
- floating HiGHS may select the candidate Farkas basis, while inherited code
  performs the load-bearing rational verification; this wrapper checks counts
  and canonical witness-free fields but is not a second rational verifier;
- the wall flag and row order are trusted from the pinned THM-2941 atlas;
- normal/optimized and serial/parallel agreement test execution stability, not
  implementation independence;
- the 19-mask containment observation is diagnostic only and is not used to
  close either row.

Those were the primary wrapper's trust boundaries.  They are now discharged
by `.scratch/lrc14_z216_gcd8_independent_audit_gcd8_independent_audit.py`,
which does not import this wrapper.  It reparses the complete atlas, checks all
`805` Farkas contradictions with a separately spelled exact verifier (and a
rejected negative-alpha mutation), rebuilds the signed rays with
`524599+70839` primitive checks, and reconstructs every terminal cell on the
full `Z/LZ` grid.  It reproduces both screen packets and both terminal hashes,
with minimum slack one and no failure.  It also proves disjointness from the
17 low-cost gcd-eight, 33 order, and 47 natural-family closures.  Its exact
transcript is
`.scratch/lrc14_z216_gcd8_independent_audit_gcd8_independent_audit.out`.

## Connection contract

```text
source:
  the two high-invoice gcd-eight wall rows in the complete z216 projected atlas;

target:
  complete enlarged terminal banks indexed by residual divisor masks and a
  unique high-slot position;

map:
  exact ray/status screen, then wall-induced exactly-one-high reduction and
  ray-supremum terminal enlargement;

preserved:
  row body, ruler, wall status, complete residual masks, divisor address,
  support-cardinality necessity, and exact open-cell inequalities;

destroyed:
  physical cover realization, owner/phase/current, endpoint origin, and all
  information outside the projected necessary atlas;

sidecar:
  positive duplicate-permitting two-high gap plus complete low-label cell bank;

decisive test:
  every enlarged one-high cell has support slack at least one.
```

## Reproduction

```powershell
python .scratch/lrc14_z216_costly_gcd8_probe_lrc_probe_agent.py --processes 2
python -O .scratch/lrc14_z216_costly_gcd8_probe_lrc_probe_agent.py --processes 2
```

Stored transcript:
`.scratch/lrc14_z216_costly_gcd8_probe_lrc_probe_agent.out`.

LF/current-byte artifact hashes at close-out:

```text
script  1eacc80c17baa70f9b11d99bade90b3cd470e9c11aee0b15e8426e7e4ceb0234
output  936b5871ac8cc616b1d52a1bd1771c2cb796e8baea593b8021533c6544fa6956
```

## Recommended integration

The requested independent reconstruction is complete and passes.  Promotion
as a scoped two-row continuation of THM-3264 is justified.  The audited
disjoint composition updates the ledger to `373184` and the `z216` wall count
to `380` without changing the cap.

The next economical wall probe is a *complete intrinsic ruler family* among
the remaining gcd `24/36/72` rows, ranked by total `sum Lr`; avoid another
externally chosen row cutoff unless it is stated only as a scratch work batch.
