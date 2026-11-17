# Chan–Pak: what has been proved, and where the next proof can happen

## The verdict

The work has gone from asking which odd counts occur to engineering **physical arithmetic on tournaments** and identifying exactly where that arithmetic stops composing. The positive existence statement at all heights remains open. No percentage of a conjecture is meaningful here: there are at least three different obligations, and progress on one does not discharge another.

1. **Exclusions and finite witnesses.** Rédei gives oddness. The in-repo [7/21 proof](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) uses insertion monotonicity, the strong-component floor, and the few small strong orders it leaves. This reproduces the obstruction in the original [2016 discussion](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament), not a new priority claim. The [corrected independent replay](../../../session203/coverage-recertification/README.md) certifies every odd through **4,294,967,295** except 7 and 21. It reran the full range without the defective F0 pass after repairing the endpoint issue. It proves neither an eventual threshold nor density one. The open statement is [Chan and Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).
2. **Reduction and exact operations.** Ordinal join multiplies Hamiltonian counts. The only infinite obligations needed for full coverage are odd primes other than 7 and 7p with p an odd prime other than 3 or 7. The singleton count and the finite exceptions for factoring, 49, 63 and 343, are witnessed. For a proof of the reduction: with two or more factors 7, use powers 7² and 7³; with one factor 7 and an odd prime other than 3,7, use 7p and factor the remainder; with only 3's accompanying a single 7, use 63 whenever there are at least two 3's. Do not infer that realizing p automatically realizes 7p. The [literal lifts](../forty_two_7p_plane_20261116/README.md) give actual sevenfold examples, and [the modular controller](../forty_two_7p_plane_20261116/ALL_MODULI.md) attains all eligible residues, but not every prescribed integer.
3. **State, not just value.** A tournament's Hamiltonian count does not tell us which operations it can survive. The [curried rooted reader](../curried_births_tooth_diagram_20261117/README.md) recovers the exact state (H,g,E,O) from four simple reads for the sink/near-sink language. Equal current counts can have different futures. The [rooted restart](../ready_restarts_rounding_20261116/README.md) already had literal 13/15 and 41/43 windows with matching continuations; the [later one-port work](../one_port_neutral_continuations_20261117/README.md) compresses the first to five vertices, classifies the smallest deletion profiles, and removes the three-run penalty for an arbitrary binary tape. That is local closure, not a universal way to choose a legal core.

## The story the failed approaches tell

The repeating mistake was to solve the wrong inverse problem. Three equations must hold **on the same physical witness**:

    land on the prescribed odd integer
    retain the required continuation data
    have enough independent choices to place the next requested integer

- **Local repair is real.** The [triangle switch](../three_visit_triangle_switch_20261117/README.md) gives a literal edge reversal for an exact +2, even under unbounded shared continuations. Changing the tile fixed an actual prescribed midpoint at 1090. The [signed ordinal construction](../signed_ordinal_unit_centers_20261117/README.md) pays a physical unit by −161+3·54=1, and reordering the tiles moves the midpoint without changing that unit.
- **The price matters.** An arbitrary binary internal separator in that fixed two-tile construction has 2^m addresses but forces at least 3^m actual paths, giving O(X^log_3(2)) possible centers through X. Coordinating two addresses can stabilize a signed gap, but the synchronized branch [pays the growth twice](../two_address_gap_selector_20261117/SELECTION.md). The [four-contact kernel](../neutral_contact_flag_interactions_20261117/README.md) further proves that neutrality for the gap does not make center digits additive. Its finite-set flag model makes the extra path runs visible rather than hiding them in a scalar.
- **Getting rid of the price is progress, not the end.** The one-root alphabet has at most two runs, and unlike the coordinated seam each bit is independently usable. Yet even here there is a difference between controlling the sign, moving a center, and covering all centers. The old [literal ternary digit](../fresh_unit_prime_gluing_20261116/ternary_deposit/README.md) realizes 3H, 3H+2 and 3H+4 with synchronized continuations on actual families. At the next round it can retain the old digit while failing to deposit an independent fresh one. A change of physical encoding can reopen the next step, as the [fresh-continuation handoff](../fresh_unit_prime_gluing_20261116/target_permuted/README.md) did at 8053, but isolated handoffs are not an all-order dispatch rule.
- **Macroscopic abundance isn't support.** The [sharp-cut Fourier estimates](../chanpak_true_cut_fourier_20261115/README.md) control presentations in physical domains. They are not asymptotics for distinct attained counts. Likewise all-moduli steering and relative bracketing do not place a target. An exact sublevel cost and a collision estimate are different tests. This is the common thread connecting the residue, Collatz and finite-field analogies to the actual problem.

## A new diagnostic: what does insisting on a permanent repair cost?

This is an exact projection test of **selected one-root continuation families**, not a new tournament coverage record. Take complete unrooted representatives at orders seven and eight from the [archived census](https://github.com/eliottcassidy2000/math/tree/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities). Root at every vertex, use literal subset DP, and retain normalized positive profiles with O−E=g and E>0. Deduplicate scalar output laws, not root labels. A law is (h,A,g) with A=2E, read by

    H(z) = h + A z + g(f(z)−1),
    f(0)=1, f(2u+1)=f(u), f(2u)=f(u)+u+(u mod 2) for u>0.

The binary address is read from low bit to high as the chronological near-sink choices. All solutions through X are exhausted within z≤(X−h)/A. Call a law matched if the checked order bank has h−2 or h+2 at the same (A,g). A global source can pad smaller hosts without changing any of these data: no rooted-start path can include it, and adding it to the complement preserves that complement's count. Thus the order banks include the projected smaller profiles. A single law counts once even if it participates in two pairs.

On the shell (2^19,2^20), the exact results are:

| Bank | Distinct laws | Odd hits / 262144 | Prime hits / 38635 | 7p hits / 6444 | Mean address presentations per odd |
|---|---:|---:|---:|---:|---:|
| Two smallest explicit paired banks | 4 | 83658 | 10627 | 1507 | 0.503 |
| Matched, order ≤7 | 39 | 174329 | 25221 | 4033 | 1.797 |
| Matched, order ≤8 | 257 | 244086 | 35801 | 5935 | 4.084 |
| Unrestricted normalized laws, order ≤8 | 2808 | 262144 | 38635 | 6444 | 35.657 |

The matched order-eight bank misses **18,058 odds**, including **2,834 primes and 509 values 7p**, in that shell, despite averaging more than four presentations. The unrestricted bank also hits every odd from 23 through 2^20 in this calculation. This last statement is strictly finite and is not a claim of unpaired coverage beyond the stated range, nor a statement that the unrestricted bank uniformly supplies repair-ready outputs. Similar fractions for the matched banks occur in the smaller measured shells. This diagnoses a cost of imposing a permanent repair at the *entrance*. It does not show that all possible matched banks are sparse. Compare [the full ledger](ledger.out).

## What I would test next

**1. Separate landing from rearming.** Use the broader exact source bank to choose an absolute count, then test whether that count admits a second physical encoding with useful continuation. Do not require permanent neutrality at the entrance. The concrete starting point is the score-preserving, same-H same-E changing-g recontact developed in the fresh-continuation work. On each of the missed prime and 7p targets above, retain different actual root/mark states, not just the first hit, and ask whether a tested exact rewrite can land on the required defect-zero profile. Then test whether that physical state accepts *all* required digits. Reaching defect zero alone does not certify readiness. A useful falsifier is a target where the full declared rewrite orbit has no eligible output state, not merely a failed random search.

**2. Treat choice capacity as a theorem to check before searching.** For a fixed core in the [one-live-port two-tape compiler](../one_port_neutral_continuations_20261117/BALANCED_TAPES.md), a chosen signed gap gives one affine equation C_R P_R−C_L P_L=d per low-bit pair. With g0=gcd(C_L,C_R), primitive steps are a=C_R/g0 and b=C_L/g0; fixing parities permits only steps (2a,2b). The absolute count has the positive cross contribution R P_L P_R, with R the optional ordered two-cover count of the interior deletion core. Counting cuts of full paths, R≥(m+1)H(core) for m interior vertices. Thus each nonempty line has at most 1+√(X/(4Rab)) candidate addresses through X, regardless of carry details. Deduplicate identical absolute-read formulas before summing capacities. If an exhaustive growing family has o(X/log X) capacity on a prime-scale shell it cannot fill either prime channel. If capacity is large and actual support remains low, examine collisions and the distribution of uncovered fibers, not another modular image. For presentation counts R_x the elementary bound |support|≥(ΣR_x)^2/ΣR_x² identifies exactly the second moment that must improve.

**3. Attack the physical cofactor inverse where arithmetic actually descends.** In the existing sevenfold route a target quotient p has form 30v+11u or 30v+13u, with u=H(B), v the optional ordered two-cover count of a physical block. A congruence fixes u mod 30, and v≥2u gives u≤p/71 or p/73. Computing v after selecting u is easy; proving an actual B exists on that line is the unresolved part. Build a bank of physical (u,v) envelopes at prescribed u, using the fact that ordinal joins multiply *both* readings. Test those lines on independent growing shells, and try to prove overlap of exact support or a source-changing induction. Mere residue surjectivity of (u,v), or bounding v above and below, would leave holes on the line.

The principle that survives these tests is **exact integer selection with renewable state**. The proof should not be judged by how many large examples it can manufacture or how elegant the arithmetic identity looks. It needs a mechanism that takes a prescribed target, selects a physical continuation state for it, and closes at the next height on every allowed branch. I would prioritize the landing/rearming separation because the diagnostic isolates a measurable cost there, and the archive already contains legal same-count re-encodings to start from.

## Reproduction of the diagnostic

Use the two files from the archived census (see link above), with the order-eight CSV converted to outgoing rows by unpacking each decimal code into eight little-endian bytes:

    python -c 'import csv; [print(*(int(r["code"])>>(8*i)&255 for i in range(8))) for r in csv.DictReader(open("n8_classes.csv"))]' > n8_rows.txt
    g++ -O3 -std=c++17 ledger.cpp -o /tmp/chanpak_ledger
    /tmp/chanpak_ledger n7_representatives.txt n8_rows.txt

The executable validates the 456 and 6880 input lengths and checks every row pair is a tournament. Completeness of those class lists is inherited from the cited archive, not re-proved by the diagnostic.
