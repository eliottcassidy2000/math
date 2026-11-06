# Exact integer index for four paired O/Z orbits

**Proved and independently audited.** Every word in each of four physical O/Z carrier orbits has an exact one-integer arithmetic parameterization. The resulting complete finite inverse improves the upper density bound for their union to less than 0.75 percent of odd bases.

This result does not prove the staggered profile criterion, prime-track coverage, or any eventual Hamilton-path threshold. The upper bound applies only to the four specified orbits.

- [Proof and scope](PROOF.md)
- [Independent review](audit/REVIEW.md)
- [Frozen physical and arithmetic dependencies](dependencies/PROVENANCE.json)
- [Original reviewed proof](audit/reviewed-inputs/PROOF.md)
- [Status-only promotion diff](STATUS_PROMOTION.diff)

## Reproduction

From this directory:

    python3 verify.py > /tmp/paired-index-producer.json
    cmp /tmp/paired-index-producer.json NORMAL.json

From the audit directory:

    python3 check.py > /tmp/paired-index-independent.json
    cmp /tmp/paired-index-independent.json NORMAL.json
    g++ -O3 -std=c++17 check_physical.cpp -o /tmp/paired-index-physical
    /tmp/paired-index-physical > /tmp/paired-index-physical.json
    cmp /tmp/paired-index-physical.json PHYSICAL_NORMAL.json

The preserved audit manifest covers all 18 reviewed files, with no local binaries. The original producer manifest is preserved within the frozen reviewed inputs. The package-level manifest additionally covers the promoted proof and release documentation.
