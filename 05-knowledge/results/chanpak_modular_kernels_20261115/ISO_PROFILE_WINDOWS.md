# A fixed response, with freely chosen low bits in the reflected count

This complementary construction addresses the “change the subtractive coordinate without moving the baseline” requirement exactly, but only **locally at 2**. Four physical strong modules from the order-nine [profile census](HYPERBOLIC_MOVE.md) suffice. Each pair shares size and the complete two-visit positioning coordinate (Π₂), while its (H) differs:

| Size | (H) | (Π₂) | lexicographic tournament code |
|---:|---:|---:|---:|
| 8 | 217 | (2^{12}) | 27280129 |
| 8 | 223 | (2^{12}) | 10760449 |
| 9 | 695 | (2^{14}) | 10342203650 |
| 9 | 707 | (2^{14}) | 31344104449 |

A code bit is one when the smaller vertex beats the larger, with pairs ordered lexicographically and the first bit least significant. The small modules are checked independently by direct subset DP. Put one order-eight choice and (m) order-nine choices in an ordinal join, and optionally pad with (t) transitive vertices. Its (Q) is exactly (2^{12+14m+t}) regardless of every module choice. Apply the literal anti-seed from THM-5130(D4) as in the reflected cone of the fixed-kernel note. Then for (e∈{0,1}), (0≤j≤ m), we have the physical counts
$$K=22+14m+t,\qquad N_{e,j}=2^K-217^{1-e}223^e\,695^{m-j}707^j.$$
All the tournaments at fixed ((m,t)) have the same order. The split modules themselves need not be consecutive on paths through the outside anti-seed; the two-visit profile formula, not that assertion, computes the count. The placement baselines are literally identical.

## Exact local decoder

Modulo (2^r), (r ≥ 2), the ratio (u=707/695) has (v_2(u-1)=2). Successive squaring shows that (u) generates the full subgroup (1 mod 4) of order (2^{r-2}). The other ratio (v=223/217 ≡ 3 mod 4) supplies its missing coset. Given a target odd residue for (N), divide its negative by the common odd factor (217 · 695^m). The class mod 4 selects (e); a discrete logarithm to the base (u) selects a unique (j mod 2^{r-2}). As soon as (m≥2^{r-2}-1), take the least nonnegative (j) and realize it *without altering (K)* or the exterior readout.

In particular for every (K≥22+14(2^{r-2}-1)), choose (m=⌊(K−22)/14⌋) and pad with the remainder. Every odd residue mod (2^r) then has a witness in the same upper band, at relative distance at most
$$\frac{223}{2^{22}}\left(\frac{707}{2^{14}}\right)^m2^{-t}.$$
below (2^K). For (r=1) all outputs are already odd. Thus the highest few bits and a growing set of lowest bits can be prescribed at once, with a single common factored response. This complements the broader fractional windows of [JOINT_WINDOWS](JOINT_WINDOWS.md), where the baseline itself moves.

## Why this is not an interval of absolute offsets

At fixed (m,t) the deflections are (217 · 695^m) times the union of the two geometric meshes (u^j) and (v u^j) on the positive real line. Here (1<u<v<u^2), and in fact the three ratios (u,v/u,u^2/v) are all greater than (201/200), by integer cross-multiplication. Any two distinct adjacent deflections therefore differ by more than (217 · 695^m/200). Arbitrarily many controlled low bits need not put even consecutive odd counts into the same reflected window. This reduces the open constructive step to a sharper one: find exact-response replacement blocks whose **real** deflection differences shrink rather than grow under composition, or change the baseline while retaining these low-bit controls.

Run `python verify_iso_profile.py` for the literal four module profiles, the exact real separation inequalities, and the finite residue controls. The valuation argument above supplies the all-order extension.
