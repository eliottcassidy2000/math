---
id: HYP-3543
title: The tournament metagraph spectrum (klein THM-584: eigenvalues d-2k, multiplicities = per-level Burnside orbit-counts of the SIGNED S_n action on the arc-hypercube) and the LRC cap co-emptiness matrix M (HYP-3538: R-even/R-odd split) are ONE structure -- the project's involution R governs the spectral multiplicities of both as Burnside orbit-counts; merged metagraph G_n/Z_2 and the LRC half-domain (0,1/2) are both the R-even projection; the level-multiplicity sequence is a SIGNED-Burnside count (not the plain graph-by-edges A000088), which is why HYP-3540's closed form is non-standard
status: SYNTHESIS (connects klein THM-584 / HYP-3540 to mac-mini HYP-3538); the LRC co-emptiness Burnside dims VERIFIED, the signed-Burnside observation explains HYP-3540's difficulty.
source: mac-mini-2026-06-29-S11
related:
  - THM-584   # klein: complement = antipodal map; metagraph = S_n-quotient of Q_d; eigenvalues d-2k
  - HYP-3540  # klein: metagraph level-multiplicity = Burnside count (closed form OPEN)
  - HYP-3538  # mac-mini: the LRC cap M splits R-even (+) R-odd; obstruction is R-odd
  - THM-583   # the witness half-system (the third realization)
  - HYP-3242  # cap = Euler char of the cover nerve
  - OPEN-Q-108
---

# HYP-3543 -- the metagraph and the LRC cap are one Burnside spectrum

## The user's hypothesis, proved (klein) and connected (here)
klein's THM-584 proves exactly "the metagraph eigenvalue multiplicities are per-level Burnside
orbit-counts": the labeled arc-flip graph is the hypercube `Q_d` (`d=C(n,2)`, one bit per pair);
the complement `R` (T->T^op) is the ANTIPODAL map (all bits flip, acting `(-1)^k` on level `k`); the
iso-class metagraph adjacency is the `S_n`-quotient of `Q_d`, with eigenvalues `d-2k` and multiplicity
`mult(d-2k) = dim(S_n-invariants at level k) = a Burnside orbit-count`.  The `R`-even/`R`-odd
(`eps=+-1`) split = the EVEN/ODD-level split, `dim R-even = (A000568+SC)/2 = V_merged`,
`dim R-odd = (A000568-SC)/2`.

## The LRC cap is the SAME structure (the connection)
HYP-3538's LRC cap co-emptiness matrix `M` (6x6 on the inner sectors) commutes with the project's
involution `R` -- here the SECTOR reflection `(1 5)(2 4)` fixing `3,6` (the `t->-t` complement on
sectors).  `R` has exactly **4 orbits** on the 6 sectors (`{1,5},{2,4},{3},{6}`), and by Burnside
> `dim(M_even) = 4 = #(R-orbits on sectors)`,  `dim(M_odd) = 2 = #(2-element R-orbits)`.
So the LRC cap's `R`-eigenspace dimensions ARE Burnside orbit-counts of `R`, exactly as klein's
metagraph `R`-even dimension is `#(even-level orbits)`.  Three realizations of one `R`-spectral
structure:
- **metagraph** (klein THM-584): `R` = antipodal of `Q_d`; mult = per-LEVEL Burnside counts; merged
  `G_n/Z_2` = `R`-even projection.
- **LRC cap** (HYP-3538): `R` = sector reflection; `dim M_even` = `#R`-orbits; the half-domain
  `(0,1/2)` = `R`-even projection; the obstruction is `R`-odd (`M_odd`).
- **witness** (THM-583): `R` = reversal; `f` = the half-system count, `R`-odd stored in `phi`.
In every case `R`-even = the SOS/Brouwer bulk (the Perron mode is `R`-even), `R`-odd = the
sign/Borsuk-Ulam obstruction.  This is the SPECTRAL form of the two-index split (THM-582).

## Why HYP-3540's closed form is non-standard (a contribution)
The naive guess `mult(d-2k) = #(S_n-orbits of k-subsets of the d pairs) = #(graphs with k edges) =
A000088 row` is WRONG: those sum to `A000088(n)` (11,34,156 for n=4,5,6), but the metagraph
multiplicities sum to `A000568(n)` (4,12,56).  The reason: the `S_n` action on TOURNAMENTS carries
BIT-FLIPS -- a vertex swap `(i j)` reverses the pair `{i,j}` -- so it is a SIGNED permutation action
(a subgroup of the hyperoctahedral `B_d`), not the plain coordinate-permutation `S_d`.  Hence
`mult(d-2k)` is a SIGNED-Burnside / hyperoctahedral orbit-count, which is why it is not a standard
graph-by-edges OEIS row.  Klein's data (level-mult, summing to A000568):
- n=4: `k=0,2,3,4 -> 1,1,1,1`
- n=5: even `k=0,2,4,6,8 -> 1,1,4,3,1`; odd `k=3,5 -> 1,1`
- n=6: even -> `1,1,5,10,13,4,1`; odd -> `1,5,8,6,2`
Target for HYP-3540: the cycle-index of the SIGNED `S_n` action on `Q_d` (Polya over `B_{C(n,2)}`
restricted to the vertex-induced signed permutations), evaluated per level.

## Why it matters for LRC
This is the spectral backbone of the whole tournament<->LRC bridge: the same `R` that makes the
metagraph spectrum a Burnside count makes the LRC cap's obstruction the `R`-odd eigenspace.  The
metagraph `R`-even bulk (provable, Perron, `(A000568+SC)/2`-dim) is the tournament-side avatar of the
LRC cap's S75e cosine-SOS bulk; the metagraph `R`-odd part (`(A000568-SC)/2` = #NS pairs) is the
avatar of the LRC obstruction `M_odd`.  Proving the LRC cap bound = bounding `M_odd`, and `M_odd` is
the LRC instance of the metagraph's odd-level (signed-Burnside) spectrum.
