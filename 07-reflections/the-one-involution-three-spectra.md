# One involution, three spectra: the metagraph, the cap, and the witness

*mac-mini-2026-06-29-S11. The owner asked to consider the metagraph eigenvalue multiplicities as per-level Burnside orbit-counts, and to connect prior tournament work to the LRC. klein-2026-06-29-S1 proved the metagraph half (THM-584); this connects it to the LRC and the witness, and notes why the multiplicity sequence is hard.*

## The hypothesis was right, and klein proved it

The owner's hypothesis — metagraph eigenvalue multiplicities are per-level Burnside orbit-counts — is now a theorem (klein THM-584). Encode a labeled tournament as a vertex of the arc-hypercube `Q_d`, `d = C(n,2)`, one bit per pair. The complement `T ↦ T^op` flips every bit: it is the **antipodal map** of `Q_d`, acting as `(-1)^k` on hypercube level `k`. The iso-class metagraph is the `S_n`-quotient, its eigenvalues are `d − 2k`, and the multiplicity of each is the dimension of `S_n`-invariants at level `k` — a Burnside orbit-count, per level. The `R`-even / `R`-odd (`ε = ±1`) split is the even/odd-level split, with `dim R-even = (A000568 + SC)/2 = V_merged` and `dim R-odd = (A000568 − SC)/2 = #NS pairs`, and the Perron mode `d` always `R`-even.

This is the *metagraph-side* realization of the `R`-eigenspace organizing principle that the LRC cap obeys (HYP-3538). So the project now has the same spectral skeleton on two sides — and, with the witness half-system (THM-583), three.

## One `R`, three spectra

Every one of the project's "two-index" phenomena is the `ε = ±1` eigenspace split of a single involution `R`, realized as a *spectrum* in three places:

| realization | the space | `R` is | mult / dim | `R`-even (the bulk) | `R`-odd (the obstruction) |
|---|---|---|---|---|---|
| **metagraph** (klein THM-584) | arc-hypercube `Q_d` / `S_n` | the **antipodal** map | per-LEVEL Burnside orbit-counts | even levels, `(A000568+SC)/2`, Perron `d` | odd levels, `(A000568−SC)/2`, `#NS` pairs |
| **LRC cap** (HYP-3538) | the 6 inner sectors | the **sector reflection** `(1 5)(2 4)` | `#R`-orbits | `dim M_even = 4 = #orbits`, Perron | `dim M_odd = 2`, the signed obstruction |
| **witness** (THM-583) | the `(p−1)/2` pairs | the **reversal** `φ∘reverse` | the half-system | the half-data | stored in `φ`, regenerates the second half |

In all three the `R`-even part is the SOS/Brouwer bulk (the Perron mode is `R`-even everywhere), and the `R`-odd part is the Borsuk–Ulam/sign obstruction (THM-582's two indices). The merged metagraph `G_n/Z₂`, the LRC half-domain `(0,1/2)`, and the witness half-system are the *same move* — the projection onto the `R`-even invariants — and the warning is always the same: do not silently discard the `R`-odd coordinate (it is the obstruction, and on the witness it is stored in `φ`). I verified the LRC entry directly: the sector reflection has four orbits `{1,5},{2,4},{3},{6}`, and `dim M_even = 4` is exactly that Burnside count.

So "the metagraph eigenvalue multiplicities are per-level Burnside orbit-counts" is not an isolated fact about tournaments — it is the spectral form of the principle that has been organizing the LRC proof, and the LRC cap's obstruction `M_odd` is the LRC instance of the metagraph's odd-level (signed-Burnside) spectrum. Bounding the cap *is* bounding that odd-level spectrum.

## Why the multiplicity sequence is hard (and the right target)

The natural guess — `mult(d − 2k) = #(S_n\text{-orbits of }k\text{-subsets of the }d\text{ pairs}) = #(\text{graphs with }k\text{ edges}) = A008406`-row — is wrong, and the reason is illuminating. Those counts sum to `A000088(n)` (`11, 34, 156` for `n=4,5,6`), but the metagraph multiplicities sum to `A000568(n)` (`4, 12, 56`). The discrepancy is the **bit-flip**: a vertex swap `(i\,j)` does not merely permute the pair `{i,j}`, it *reverses* it. So the group acting on `Q_d` is not `S_d` on coordinates but a **signed** permutation group — the vertex-induced subgroup of the hyperoctahedral `B_{C(n,2)}`. The metagraph multiplicity sequence (klein's open HYP-3540) is therefore the per-level evaluation of the **cycle index of that signed action**, not a plain graph-by-edges row — which is exactly why it does not match a standard OEIS sequence. The right target for HYP-3540 is the Pólya/Burnside generating function over the signed `S_n`-action, level by level. The level data (`n=5`: even `1,1,4,3,1`, odd `1,1`; `n=6`: even `1,1,5,10,13,4,1`, odd `1,5,8,6,2`) is the signature of that signed cycle index.

The same `±√n = ` Gauss-sum that appears in the Paley DRT spectrum (THM-586), in the imaginary `Q(√−7)` of the LRC apex, and in the `(-1)^k` antipodal eigenvalues here — the involution `R` and its sign are one object wearing the project's many hats.

See [[the-pm-one-eigenspace-of-reversal-is-the-whole-split]] (HYP-3538, THM-583), [[the-two-indices-redei-is-odd-lonely-is-even-half-tiling-is-the-quotient]] (THM-582). klein: THM-584 (metagraph antipodal spectrum), HYP-3540 (the open level-multiplicity sequence). This: HYP-3543.
