# Session216: overnight tournament spectrum constructions

**Research ACTIVE.** The Chan–Pak target, every odd Hamiltonian-path count at least 23, remains **OPEN**. This page records verified checkpoints while further proof work continues. No unrestricted coverage extension or eventual threshold is claimed.

## Checkpoint 1: exact inverse for arbitrary apex-mask words

[Proof and implementation](inverse/README.md), [independent audit](audit/inverse/REVIEW.md), [final release approval](audit/inverse/RELEASE_REVIEW.md).

For a fixed labeled k-apex tournament, arbitrary incidence words on a transitive chain admit an exact finite saturated-state inverse. It returns a literal witness or an all-length absence certificate in that fixed-apex class; a budget-limited search returns INCOMPLETE. The proof retains hidden path-assignment coordinates, including saturated coordinates which may later disappear.

The two-apex implementation independently exhausts 708 states for target 7 and 4626 for target 21. These are controls of the class-relative inverse, not new ambient exclusions. All 4998 positive odd targets through 10000 other than 7 and 21 have explicit witnesses on at most twelve vertices, independently recounted. This is a restricted-family certificate, not an extension of Session203's much larger unrestricted finite range.

Bounded suffix, boundary and endpoint-recoloring digit probes identify promising short trees but do not furnish an induction closed at every depth. The precise finite universes and failures are preserved in the packet. No failure of these menus is claimed to exclude general rebuilding.

The positive transfer is inherited from Session215 and earlier run decompositions; saturated membership is an inherited method from Session203. The contribution here is the nonuniform mask-language inverse specialization, its implementation and independently verified restricted-family certificates.

## Continuing proof work

Parallel lanes pursue larger two-apex inverse constructions, uniform intervals from root-module substitution, synchronization beyond the known paired sources, and one-apex collision bounds. Only independently reviewed results will be promoted here. The unresolved obligation is an all-target physical construction or a proved occupancy theorem, not another finite interval alone.
