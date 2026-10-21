# Balanced deletion: group fibres and two cut ports

2026-10-21. **PROVED, family-specific; independently reviewed.**

This packet studies literal centre-flanking deletion in the two-return family
T_k. It does not classify the unrestricted tournament metagraph and does not
identify the deletion fibres with quadratic class groups.

## Results and scope

- [Main proof](BALANCED_FIBRE_GROUP_AND_PORT_DEFECT.md): the 4/2/1 deletion
  fibres are Klein-four orbits with two distinct size-two stabilizer lines.
  Reflection exchanges those lines and extends the action to a faithful
  dihedral group of order eight.
- The actual path-reversal target does not descend through an undecorated
  deletion fibre. Two explicitly retained cut-port pairs give the exact
  labelled correction. A same-fibre hostile has nonisomorphic deleted
  targets, even after converse-merging.
- [All-depth proof](ALL_DEPTH_TWO_SEAM_ADDENDUM.md): arbitrarily many
  consecutive balanced insertions still require at most those same two
  cut ports. Exact defect counts are proved for every k >= 2 and d >= 1.

These are path-set actions and reconstruction laws, not claims that every
path-set permutation is a legal metagraph move. The labelled-neighbour
obstruction must not be strengthened to an unlabelled non-neighbour claim.
The independent audit retains a countercontrol to precisely that overstatement.

## Verification

The [independent review](review/AUDIT.md) checks 5,709 literal DFS paths,
7,457 all-depth path/deletion cases, and 174,720 general tournament/deletion
cases. Producer and independent checks replay identically normally and with
Python optimization. These are finite controls of the proofs, not their
all-order justification.

From this directory:

    python3 -B check_balanced_group_ports.py
    python3 -B -O check_balanced_group_ports.py
    python3 -B check_all_depth_ports.py
    python3 -B -O check_all_depth_ports.py
    cd review
    python3 -B independent_literal_review.py
    python3 -B -O independent_literal_review.py

The immutable producer proofs, programs, and receipts are retained. The
review folder includes its own frozen sources, normal/optimized receipts,
[AUDIT_PINS.json](review/AUDIT_PINS.json), and SHA256SUMS. Exploratory probes
not covered by the final audit are deliberately excluded.

## Inheritance

- [Two-return decoder](../strong_blue_games_20261020f/double_return_counts.md).
- [Earlier pointed-product and 4/2/1 deletion laws](../zero_heavy_ports_20261020g/microcosm.md).
- [Earlier interior thinning and endpoint-growth kernel](../balanced_fibre_endpoint_extensions_20261020/BALANCED_DELETION_EXTENSIONS.md).

SOURCE_PINS.json records the source checkout and historical evidence hashes.
No reserved theorem or hypothesis stub is used as a proved dependency.
