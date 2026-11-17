# Let the target choose a core: decimal scales and their exact phase families

This session continues the [two-flank inverse](../chanpak_archive_target_windows_20261117/README.md). It produces actual alternative sources for several prescribed misses, and exact structural families related to powers of ten. Neither the catalog nor the decimal constructions imply unrestricted Chan–Pak.

## Take the core from the target

The fixed even-flank formula retained ordinary marked covers but set adjacent signed responses to zero. The complete [signed two-flank read](FLANKS.md) permits either parity on each flank. Once the target fixes a parity choice (u mod2,w mod2), matching two repairs needs only the same ordinary marked reads L,R and corrected intercepts two apart:

    intercept = H(D+apex) − [u odd] R^-(D) − [w odd] L^-(D).

Here L^-,R^- are the signed counts by parity of the marked path length. If the full signed profile matches, the +2 choice survives both flanks at every height. The new catalog does not freeze D: it tries exact candidates for its physical core masks, using the smaller flank to bound the binary inverse. When asked for an entire unit window, it additionally requires a second mask with the same complete signed profile. A signed phase resolution is always followed by checking the exact full equality.

- The old two-flank core missed 65537. A directed triangle center and the address pair (65,83) realize it exactly.
- The five-vertex center (18,12,17,5,10), masks 5 and 24, with flanks (0,8064), realizes **458759 and 458761 on the same frame**. Here 458759=7·65537, another old miss. A direct full subset DP independently verifies the two graphs.
- The catalog also constructs 10007 and 10009 with different addresses on the triangle center, and chooses exact cores for 100003 and 1000003. These are simultaneous families of separate witnesses, not all unit switches. The recovered billion-prime pair still uses one pair of flanks.

[Catalog replay](catalog.py) selects directly from those prescribed targets and [literal replay](replay_flanks.py) checks the core profiles and individual test graphs. To run the independent large checks: `g++ -O3 -std=c++20 full_dp.cpp -o full_dp && python replay_flanks.py | ./full_dp` (the C++ reader ignores nonnumeric status lines).

## What the powers-of-ten work adds

[The decimal packet](DECIMAL.md) has three exact results:

1. Several literal cores have the same P₂=10000 but different H. The target can choose which cores to multiply into its decimal offset while holding the leading decimal scale fixed. For example the prime **100098587=10^8+311·317** selects those two cores exactly.
2. Repeated (11,100) blocks with uniform apex contacts have normalized counts converging to a product of two Cantor sets. The same grammar leaves a provable absolute gap after each power of 100, so self-similarity is not coverage.
3. A different, binary-apex seam realizes whole families based on 10^k for every k. Its candidate prime offset has normalized values distributed in nested dyadic bands. This extends the old decimal-window lesson: keep exponent phase, not a heuristic about the digits becoming independent.

Run `python replay_decimal.py` for the literal and arithmetic checks. The operative distinction is between freedom to choose **which family** realizes a target and proof that some family must work. The former is now implemented; the latter is the remaining proof obligation.
