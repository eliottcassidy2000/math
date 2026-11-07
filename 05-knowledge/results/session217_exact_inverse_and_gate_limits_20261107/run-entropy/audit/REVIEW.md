# Independent review: Session217 run-entropy packet

Verdict: PASS within the explicitly scoped claims in ../session217-occupancy/PROOF.md as read on2026-11-07. No all-target occupancy theorem is proved or implied.

The digit recurrence and summed run formula check algebraically. Leading zeros preserve v and delta. The prefix before a one-run ends in0, giving0<=u<1/2 and the safe3r/4 bound. The alternating-word formula gives a genuine unbounded logarithmic defect. Expansion of alpha gamma/Q gives all four terms of B with the stated signs and denominators; empty padding is handled correctly.

The low-run count is valid for all legal template lengths, not just a fixed cap. Positivity bounds k=a+b<=floor(log2 X). The inserted zero separator preserves exactly the sum of the two run counts and is injective for a fixed split. The number of n-bit strings with j one-runs is binom(n+1,2j); summing splits produces exactly the stated40(L+1)(L+2) upper factor. Dropping the separator-zero condition overcounts, safely. Entropy with rho<1/4 gives exponent H2(2rho)<1; polynomial logarithmic factors enter the o(1) exponent. At rho=0 the fixed-run bound handles the endpoint. The quantifier is a zero-density subset of targets with at least one low-run recipe, not an assertion that every individual integer needs many runs or that high-run padding suffices.

Independent checks use check_entropy.py:8191 binary words through length12; exact one-run count histograms;100 alternating lengths;76,880 rational core-product instances. All PASS. Receipt ENTROPY_CHECK.json. These are independent algebra/combinatorics checks, not new whole-tournament HP counts. The inherited graph construction is independently exercised by the separate combined inverse packet's60 small subset-DP controls.

Novelty is correctly limited: the packet cites the pre-existing one-apex reversal-defect mechanism and presents the new scope as the fixed-three-core run-budget obstruction, without a general literature-priority claim. No issue found requiring theorem repair.

The added finite-bank extension in section4 also passes: for every nonempty tournament core q>=2h>=2; nonempty-cut A>=h>=1. Hence the same padding-length bound holds. Replacing80 by the fixed number B of core-cut states gives the exact prefactor B(L+1)(L+2)/2. Strong-core legality has the same sufficient nonempty-cut/source-prefix condition. The entropy exponent and zero-density conclusion are unchanged, with bank-dependent constants only.
