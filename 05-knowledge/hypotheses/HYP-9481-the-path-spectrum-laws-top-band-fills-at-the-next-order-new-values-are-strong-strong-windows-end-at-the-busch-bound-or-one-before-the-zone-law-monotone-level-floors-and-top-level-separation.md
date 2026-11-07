---
id: HYP-9481
title: "Six laws of the Hamiltonian-path spectrum of tournaments, exact through order 11 (THM-4889), conjectured for all orders: (a) every odd value below max H(n) that is not attained at order n is attained at order n + 1, so the least order of h is m(h) = min{n : max H(n) >= h} or m(h) + 1; (b) every value first attained at order n is the count of a strong n-tournament; (c) the strong window of h is an interval ending at the Busch bound B(h) = max{n : f(n) <= h} or one order before it; (d) the zone law: every n-tournament with more cyclic triangles than any (n-1)-tournament has more Hamiltonian paths than any (n-1)-tournament, phi_n(Lambda(n-1) + 1) > max H(n-1); (e) the least path count at level k is nondecreasing in k at every order; (f) top-level separation: every tournament with the maximum number of cyclic triangles has more Hamiltonian paths than every tournament with fewer, phi_n(Lambda(n)) > psi_n(Lambda(n) - 1)."
status: >
  OPEN (klein-2026-11-06). FINITE-EXACT at the orders stated in THM-4889: (a) orders 6-10 -> next order (the 208 band
  values of order 10 all fill at 11); (b) orders 3-11; (c) odd h < 375; (d), (e) orders 5-11; (f) orders 6-11 (the
  pre-registered order-11 test P12 passed: phi_11(55) = 91,395 > psi_11(54) = 89,357, margin 2.3 %). All six laws
  were pre-registered for order 11 and confirmed there (packet section 6). No proof route is known for (d), (e), (f): they ask for a lower bound
  on Hamiltonian paths in terms of cyclic triangles, and the flip gradient (a level-raising flip raises h in 98 % of
  cases, not always) shows that no monotone local argument works. (a) together with the Chan-Pak completeness
  conjecture (OPEN, MISTAKE-559) would make n(h) a function of A038375 alone; (a) is sharper than HYP-9452 (2), whose
  band exception it fixes at exactly one order. (b) says products of smaller counts are never new. (c) refines
  THM-4744's inventories; the 19 windows below 225 that end one order early are listed in THM-4889 (3).
session: klein-2026-11-06 (codex/session-spectrum-20261104)
source: 05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/README.md
related:
  - 01-canon/theorems/THM-4889-the-hamiltonian-path-spectrum-by-order-and-by-value-exact-censuses-to-order-11-resolved-by-cyclic-triangles-least-orders-strong-windows-level-floors-and-the-zone-law.md
  - 05-knowledge/hypotheses/HYP-9452-least-orders-of-hamiltonian-path-counts-the-top-band-of-unattained-values-and-paley-counts-at-the-remaining-supersingular-orders.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md
---
# HYP-9481 — the laws of the path spectrum beyond the census

Notation as in THM-4889. The six statements, each exact at the stated orders:

- **(a) The top band fills at the next order.** If `h < max H(n)` is odd, `h ∉ {7, 21}` and `h ∉ S(n)`, then `h ∈ S(n+1)`. Equivalently `n(h) ∈ {m(h), m(h)+1}`, with `n(h) = m(h) + 1` exactly when `h` lies in the top band of order `m(h)`. Exact for the 392 top-band values of orders 6–10 (all attained at the next order); the 3,635 band values of order 11 are the first open cases.
- **(b) New values are strong.** `S(n) − S(n−1) ⊆ Ss(n)`. Exact at orders 3–11.
- **(c) Strong windows.** `W(h) = {n : h ∈ Ss(n)}` is an interval `[n(h), e(h)]` with `e(h) ∈ {B(h) − 1, B(h)}`. Exact for odd `h < 375`.
- **(d) Zone law.** `φ_n(Λ(n−1) + 1) > max H(n−1)`; margins 1.80, 1.53, 1.13, 1.39, 1.35, 1.47 at orders 5–10.
- **(e) Monotone level floors.** `φ_n(k) ≤ φ_n(k+1)` for all `k < Λ(n)`; and `φ(k) = min_n φ_n(k)` is nondecreasing (exact for `k ≤ 9`).
- **(f) Top-level separation.** `φ_n(Λ(n)) > ψ_n(Λ(n) − 1)`; margins 1.11, 1.08, 1.01, 1.02, 1.01, 1.02 at orders 6–11. At odd `n` this says every regular tournament has more Hamiltonian paths than every tournament with scores `(m−1, m, …, m, m+1)`; the margin is small, and it survived the jump from 15 regular 9-tournaments to 1,223 regular 11-tournaments (A096368), whose counts span `[91395, 95095]` against a non-regular maximum of 89,357.

**Why they might hold.** (a)–(c) are statements about the growth of the strong spectrum: strong `(n+1)`-ears of the order-`n` classes near the top (THM-4097/4102/4104's ear images are solid intervals reaching far above `max H(n)`) cover the band; the two-sided bound in (c) says the strong realisers of `h` survive one-vertex extension until the Busch floor overtakes `h`, and never die earlier than one order before that. (d)–(f) are the metagraph's level order being almost a refinement of the value order: the floor by level grows faster than the maximum over the previous order. A proof of any of (d)–(f) would be a new inequality between `c3` and `h`.

**Tests.** Order 12 (`A000568(12) = 154,108,311,168` classes) is beyond a direct census; (a)–(c) can be tested on the lower part of `S(12)` by ears from the order-11 census restricted to the band, (d)–(f) need the level-resolved top of order 12 (near-regular 12-tournaments by `gentourng -d5 -D6`).

**Note on the order-11 verdicts.** See the packet's §6 and PREREGISTERED.md: the predictions were committed before the census finished; all six laws passed at order 11. The one refuted prediction was not one of the laws: HYP-9452 (1)'s constant 0.85 for the solid prefix fails at order 11 (0.846), as the sparsity of the regular top level predicts.
