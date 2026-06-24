# HYP-2996: LRC14 Residual Section Characterization And Packet Grid Verification

**Status:** PROOF-INTERFACE / executable bounded audit; not a proof.

## Claim

The HYP-2963 labelled-packet bank admits an explicit residual-section map into
the HYP-2989 Haar-product grid:

- `Q-WITNESS` packets are exact `0`-cochain exits and map to
  `orthogonal_zero`;
- AP/Goddyn-Wong boundary atoms are declared boundary cohomology and map to
  `same_tile_indicator`;
- C27/unit-petal packets are owner-strip coboundaries and map to
  `horizontal_owner_strip`;
- positive open-Haar packets map to `vertical_owner_strip`;
- covering boundary-moment packets descend by `nested_refinement`;
- K33 packets are named THM-572 state-lift residuals and map to `cross_handoff`;
- F7 is now a testable section: hard, zero-open, non-AP/GW, and not discharged
  by owner-strip, cross-handoff, or nested-refinement grid exits.

## Evidence

Script:
`04-computation/lrc14_residual_section_packet_grid_codex_s168.py`

Stored output:
`05-knowledge/results/lrc14_residual_section_packet_grid_codex_s168.out`

The default run imports the HYP-2963 packet bank
(`single_limit=180`, `two_swap_limit=36`, `alias_depth=4`,
`lcm_tail_max=5`) and the HYP-2989 depth-3 Haar-product grid.  The Haar grid
baseline is unchanged:

- `same_tile_indicator=225`
- `vertical_owner_strip=1020`
- `horizontal_owner_strip=1020`
- `cross_handoff=2312`
- `nested_refinement=2312`
- `orthogonal_zero=43736`

Packet-grid verification:

- audited packets: `21913`
- hard `q>=14` packets: `7237`
- hard non-AP/GW packets: `7235`
- hard non-AP/GW packets with owner/cross/nested grid exit: `7235`
- zero-open hard non-AP/GW packets: `0`
- candidate F7 harmonic sections: `0`
- validation errors: `0`

Section counts:

- `direct_q_witness_section`: `14676`
- `ap_gw_boundary_section`: `2`
- `unit_petal_c27_strip_section`: `4`
- `open_haar_witness_section`: `6040`
- `covering_boundary_moment_section`: `1188`
- `k33_state_lift_section`: `3`

Packet-grid counts:

- `same_tile_indicator`: `2`
- `vertical_owner_strip`: `6040`
- `horizontal_owner_strip`: `4`
- `cross_handoff`: `3`
- `nested_refinement`: `1188`
- `orthogonal_zero`: `14676`

## Tournament Analysis

Vertices are residual sections and packet-grid exits, not runners.  The
pairwise observable is retained LRC predicate, exact scale, boundary data, Haar
grid class, state-lift visibility, moment-route visibility, named residual
status, proof exit, and anti-scalarization guard.  The switch orients toward
the stronger weighted retention vector, with the declared Hamiltonian path over
section names.

Fingerprint:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`
- `directed_3cycles=0`
- `SCC_sizes=[1,1,1,1,1,1,1]`
- `Hamiltonian_path_count=1`

## Interpretation

Within the bounded HYP-2963 bank, the residual stack has no anonymous remainder:
all hard non-AP/GW packets either expose a positive owner-strip witness, descend
to a covering boundary-moment chart, or route as named K33/THM-572 cross-handoff
debt.  This does not prove the global LRC14 reduction, but it turns "F7" from a
vague failure bucket into an executable missing-section predicate.

## Anchors

HYP-2995, HYP-2994, HYP-2993, HYP-2992, HYP-2991, HYP-2989, HYP-2988,
HYP-2987, HYP-2963, HYP-2956, THM-572, OPEN-Q-108.
