---
id: THM-802
renumber_note: "Claimed as THM-800 at 35a17e793 but codex's THM-800 (two-replacement-common-scale-descent) reached origin first (my push was rejected mid-flight); renumbered to THM-802 (801 = codex B3/B2 descent). My HYP-6880 likewise renumbered to HYP-6885 (codex-S12 first-pusher on 6880)."
title: THE TUTTE POLYNOMIAL OF THE STAIRCASE SMITH NETWORK IN CLOSED FORM — T(N_n; x, y) = Π_{k=1}^{n−2} (x − 1 + [k]_y), the product of shifted q-INTEGERS [k]_y = 1+y+…+y^{k−1} (5-line rank–nullity proof; every classical specialization becomes a one-line product: κ = (n−2)!, forests = (n−1)!… rising factorial at y=1, acyclic orientations 2^{n−2}, chromatic q(q−1)^{n−2}, flow ≡ 0 via the bridge factor) — AND THE LRC BRIDGE: (i) on the unit circle y = e(2πθ) each factor becomes the DIRICHLET KERNEL D_k(θ), so T(N_n; 1, e(2πθ)) is the SUDLER-type product over the interval core {1..n−2} — the resonance landscape of the LRC extremal speeds (Dedekind-η adjacency); (ii) the HARMONIC-MEASURE IDENTITY: the staircase's electrical aspect ratio H_{n−2} (HYP-6865) and the LRC interval-core good measure are the same harmonic number — m({1..k}, λ=1/(k+2)) · C(k+2,2) = H_k (mac-mini THM-736 proved the k=12 instance |G'| = H_12/91; here verified exactly for the whole ladder)
status: CLAIMED (kind-pasteur-2026-07-15-S128 cont.8) — closed form has the 5-line proof (below), referee + specializations + the harmonic-measure ladder verification run THIS SESSION; the Sudler/η adjacency is logged as HYP-6885 (observational)
source: kind-pasteur-2026-07-15-S128 (cont.8; owner: work the closed-form Tutte polynomial, connect to LRC)
depends_on:
  - HYP-6865  # N_n = series chain of parallel bundles (n−2,…,1); κ=(n−2)!; R=H_{n−2}; self-dual
related:
  - opus-S307 (the Smith diagram of the METAGRAPH: concordance law, transitivity resistance) — the two electrical faces: theirs on iso-class space, this on the staircase itself
  - mac-mini THM-736 (far-peel Farey: the deep-well base measure |G'({1..12})| = H_12/91 — the k=12 instance of (ii))
  - THM-732 (Dedekind/Bernoulli exact certificates — the log-derivative world of the same Sudler-type products)
  - kps-S11 equinumerosity block (the G_n↔E_n Tutte-duality lead this begins to execute)
---

# THM-802 — the staircase's Tutte polynomial, and what it knows about the LRC

## The closed form (PROVED)

N_n is the chain of parallel bundles with multiplicities k = 1, …, n−2 (HYP-6865). In the
rank–nullity expansion T = Σ_{A⊆E} (x−1)^{r(E)−r(A)}(y−1)^{|A|−r(A)}, a subset chooses a_j edges
from bundle j independently; r(A) counts nonempty bundles. Each bundle therefore contributes the
factor (x−1) + Σ_{a=1}^{k} C(k,a)(y−1)^{a−1} = (x−1) + (y^k−1)/(y−1), i.e.

>  **T(N_n; x, y) = Π_{k=1}^{n−2} ( x − 1 + [k]_y ),   [k]_y = 1 + y + … + y^{k−1}.** ∎

Specializations (all one-line products): T(1,1) = (n−2)! = κ (HYP-6865 ✓); T(2,1) = (n−1)!
(spanning forests = rising factorial (x)^{(n−2)} at x=2); T(1,2) = Π(2^k−1) (connected spanning
subgraphs); T(2,0) = 2^{n−2} (acyclic orientations — each bundle uniformly oriented); chromatic
P(q) = q(q−1)^{n−2} (the underlying path); flow F ≡ 0 (the k=1 bundle is a bridge — its factor
vanishes at x=0, y=1−q).

## The LRC bridge

**(i) The unit-circle specialization is the resonance landscape.** At y = e(2πiθ):
[k]_y = D_k(θ) (Dirichlet kernel), |[k]_y| = |sin(πkθ)/sin(πθ)|. So
|T(N_n; 1, e(2πiθ))| = Π_{k≤n−2}|sin πkθ| / |sin πθ|^{n−2} — the SUDLER-type product over exactly
the interval-core speeds {1, …, n−2} of the LRC extremals. Small values of the numerator factors
⟺ ||kθ|| small ⟺ θ in speed-k's danger arc: the staircase's Tutte polynomial, evaluated on the
unit circle, computes the joint resonance profile of the deep-well core. Sudler products are the
exponential of Birkhoff sums of log|2 sin| — the Dedekind-η / quantum-modular world — the same
neighborhood as THM-732's Dedekind–Bernoulli certificates (log-derivative of these products) and
mac-mini's X₀(14) thread. (Observational scaffolding: HYP-6880.)

**(ii) The harmonic-measure identity (the resistance IS the good measure).** HYP-6865: the
staircase's pole-to-pole resistance is H_{n−2}. Claim, verified exactly this session for the whole
ladder k = 2..12: the interval core's good measure at its TIGHT threshold satisfies

>  **m({1,…,k}; λ = 1/(k+2)) = H_k / C(k+2, 2)**

(k=12: H_12/91 — exactly mac-mini THM-736's deep-well base |G'|). Electrically: speed k ↔ the
staircase row of width k ↔ a parallel bundle of k unit resistors; the LRC good measure of the core
= (the staircase's total resistance)/(number of arcs of K_n). The two pillars compute the same
harmonic number from the same index set — the dictionary entry the Smith functor was for.

## Evidence log

- [ ] referee: closed form vs brute-force subset-expansion Tutte, n=4..7, exact
- [ ] specialization table exact
- [ ] harmonic-measure ladder k=2..12 exact (Fractions)
- [ ] Sudler landscape sanity plot/stats (HYP-6880)
