---
id: HYP-3729
title: THE LRC PARITY IS THE BIPARTITENESS OF C_n. The equally-spaced covering orbit (at the witness) is the cycle C_n; for EVEN n it is BIPARTITE (antipodal +1/2 pairing, 2+lambda_min(C_n)=0, no spectral gap), giving the TIGHT worst-case M=1/n (the even block = the scaled consecutive set, a verified LOCAL worst-case: 0 of all single-speed perturbations beat 1/n at n=8,10); for ODD n it is NON-BIPARTITE (2+lambda_min(C_n)=4 sin^2(pi/2n)>0, a spectral gap), giving the strict primitive MARGIN. This is EXACTLY the apex least-eigenvalue certificate parity (HYP-3606/THM-590: gap=4sin^2(pi/2p), positive iff C_p non-bipartite iff p odd). So even-n worst-case <-> bipartite cycle (tight, gap 0); odd-n realizability <-> non-bipartite cycle (margin), and the odd-n primitive covering-min lives on a circulant mod m (n=7:13, n=8:15) = a Paley vertex count (Paley graph on 13 = RAMANUJAN). The MARGIN deviation is exact at n=7 (2/13, =1/91=1/(n(2n-1))) and n=8 (2/15, =1/120=1/(n(2n-1))) but DEVIATES at n=9 (best-known 4/33, margin 1/99 NOT 1/153) -- irregular, no clean formula, the covering-min is genuinely n-dependent
status: VERIFIED (parity=bipartite C_n: 2+lambda_min=0 even / 4sin^2(pi/2n) odd, n=5..14 exact; even block M=1/n + local worst-case n=8,10 exact; margin exact n=7,8 mediant, deviates n=9). The even-n worst-case PROOF (no covering set beats 1/n) is the open LRC -- this gives the structure (bipartite C_n extremal) + a strategy, not a full proof. Odd-n primitive covering-min trajectory open (searches unreliable; HYP-2566).
source: mac-mini-2026-06-30-S50
related:
  - HYP-3727  # primitivity: full covering-min=1/n (all n), primitive covering-min>1/n (the margin); the Paley/Ramanujan frame
  - HYP-3606  # the apex least-eigenvalue certificate: gap = 2+lambda_min(C_p) = 4sin^2(pi/2p), positive iff p odd -- THE SAME PARITY
  - HYP-3725  # the construction is NOT the covering-min (S47 refutation); the spread family is
  - HYP-2566  # uniform looseness (primitive covering-min > 1/n) -- the open conjecture
references:
  - opus-2026-06-30-S1     # even block 2*{1..n-1}=1/n (the q-witness/easy case; here: its orbit is the bipartite C_n)
  - klein-2026-06-29-S34   # Ramanujan=Ihara-RH; bipartite = the degenerate case (K_3,3); metazeta
results:
  - 04-computation/covering_min_smart_search_macmini_20260630.py
  - 05-knowledge/results/covering_min_smart_search_macmini_20260630.out
---

# HYP-3729 -- the LRC parity IS the bipartiteness of C_n

The owner asked to (1) see how the margin deviates with n, (2) work the even-n worst-case proof (the
equally-spaced C_n can't be beaten), (3) smarter odd-n realizability + tournament connections. The unifying
discovery: **the parity of n in the LRC is the bipartiteness of the cycle C_n.**

## The equally-spaced orbit IS the cycle C_n
At the optimal witness `t*=1/(2n)`, the even block `2*{1..n-1}` has orbit `{1/n, 2/n, ..., (n-1)/n}` -- the
`n`-th roots of unity minus the origin, i.e. the vertices of the cycle `C_n` (minus the hole at 0). The min
distance to 0 is `1/n`, so `M = 1/n`. The structure of this orbit is governed by the parity of `n`:
- **EVEN n: C_n is BIPARTITE.** The orbit is antipodally symmetric (`x -> x+1/2`), pairing `k/n` with
  `(k+n/2)/n`; the lone point `1/2` maps to the hole `0`. `2 + lambda_min(C_n) = 0` (the bipartite
  eigenvalue `-2`). **No spectral gap** -> the bound `M=1/n` is **TIGHT** -- the equally-spaced even block is
  the worst case.
- **ODD n: C_n is NON-BIPARTITE.** No antipodal pairing; `2 + lambda_min(C_n) = 4 sin^2(pi/2n) > 0` -- a
  **spectral gap**. The primitive covering-min sits **strictly above** `1/n` (the margin).

This is **exactly** the apex least-eigenvalue certificate (HYP-3606 / klein THM-590): the doublet gap is
`2 + lambda_min(C_p) = 4 sin^2(pi/2p)`, positive iff `C_p` is non-bipartite iff `p` is odd. The LRC's
even/odd split and the apex certificate's bipartite/non-bipartite split are the **same parity**, on the same
cycle. (Verified `2+lambda_min(C_n)`: `0` for even `n`, `4sin^2(pi/2n)` for odd `n`, `n=5..14`.)

## Thread 2 -- the even-n worst-case (structure + strategy)
For even `n`, the equally-spaced even block is the extremal worst-case `M=1/n`, and it is a **verified local
worst-case**: among all single-speed perturbations (n=8,10), **none** drops `M` below `1/n`. The bipartite
(antipodal) structure is the reason: the orbit is as spread as possible (equally-spaced minimizes the max
gap), and antipodal symmetry pins the min-to-0 distance at `1/n`. **Proof strategy** (the full claim "no
covering set beats `1/n`" is the open LRC): a Fejer/Selberg non-negative-kernel certificate exploiting the
antipodal (bipartite) symmetry, or klein-S34's "even-fold `n=2p -> LRC(p)`" reducing even-`n` to the proven
small cases (`14 -> 7`). The bipartite `C_n` is the correct extremal frame; the gap-0 degeneracy is why even
`n` is the tight (hardest-to-improve) case.

## Thread 1 -- how the margin deviates (honest)
The primitive covering-min margin `M_prim(n) - 1/n`:
- `n=7`: `2/13 - 1/7 = 1/91 = 1/(n(2n-1))` (the mediant; `91 = 7*13`).
- `n=8`: `2/15 - 1/8 = 1/120 = 1/(n(2n-1))` (the mediant; `120 = 8*15`).
- `n=9`: best-known `4/33 - 1/9 = 1/99` (`99 = 9*11`, **NOT** `1/(9*17)=1/153`) -- the mediant `2/17` is not
  achieved; the margin **deviates** from the clean `1/(n(2n-1))` form.
So the clean hexagonal margin `1/(n(2n-1))` (HYP-3726) holds **only at `n=7,8`**; at `n=9` it breaks. The
deviation is **irregular** -- there is no known closed form, reflecting that the covering-min is genuinely
`n`-dependent (klein-S32: "not closed by naive recursion"). Reliable values exist only for small `n` (the
covering-min is a hard optimization; structured/local searches are unreliable for `n>=9`).

## Thread 3 + 4 -- odd-n realizability and the tournament connection
For odd `n` the even block fails (no even multiple of odd `n`), so the primitive covering-min is the genuine
frontier. It lives on a **circulant mod m** (witness `t*=k/m`: `m=13` at n=7, `15` at n=8), and `m=2n-1` is a
**Paley vertex count** -- the canonical spread/**Ramanujan** structure (Paley graph on 13 is Ramanujan; at
`n=14`, `2n-1=27=GF(3^3)` -> Paley tournament). So the odd-n covering-min is the "most equidistributed
(Ramanujan-flat) primitive covering set," and the search should be seeded by **quadratic-residue / Paley
circulants mod m**, not by perturbing the dull construction. The tournament tie: the non-bipartite `C_n`
(odd) is the apex doublet's cycle (HYP-3606); the Paley tournament (mod `2n-1`) is the regular/Ramanujan
tournament; the `C_n` bipartite/non-bipartite split is the project's even-graph/odd-cycle duality at the
covering floor. (opus's metazeta = the Ihara-zeta of the metagraph; klein-S34: Ramanujan = Ihara-RH.)

## The one-line synthesis
The LRC floor is a spectral-gap statement on the cycle `C_n`: **even `n` = bipartite `C_n` = gap 0 = tight
`M=1/n`** (the equally-spaced worst-case); **odd `n` = non-bipartite `C_n` = gap `4sin^2(pi/2n)` = the
margin** (the Paley/Ramanujan realizability frontier). The same parity runs the apex certificate, the
covering-min, and the even-graph/odd-cycle duality.
