---
id: THM-4504
title: "Moran pressure, fair-word rise bounds, finite integer windows, and model families of 27"
status: >
  PROVED WITH STATED HYPOTHESES (M, repaired R, R', S for
  1<beta<log_2 3, G, D with the open hitting-set definition, B under
  regular variation and a limiting residue fraction). CITED MODEL
  (Lagarias-Weiss and Kontorovich-Lagarias). FINITE-EXACT attributed to
  the incoming exhaustive run through 2^32; that scan was not rerun in
  the poset session. EMPIRICAL for full-orbit and individual-tree rates.
  Correction 2026-09-26: finite conditional means require positive hit
  probability; a supremum need not be attained; prefix-cover and integer
  height exponents differ by the scale conversion. Actual backward-tree
  densities and actual record asymptotics remain OPEN.
source: collatz-procgen-20260922 session, family27 lane (2026-09-26), answering the owner's question about numbers beyond 27 and how their occurrence rate governs fractal recursion; audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md
related:
  - 01-canon/theorems/THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality.md (g_q(2) = (1+q)/4)
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (3 + 1 = 4)
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md (peak discount)
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md (the owner's inverse tree)
note: 05-knowledge/results/procgen_family27_20260926_long_orbit_families.md
scripts: 04-computation/experiments/procgen_family27_20260926_{run,theory,bfiles,reference}.py and procgen_family27_20260926_scan.c
script_audit: 04-computation/experiments/procgen_family27_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_family27_20260926.out
output_audit: 05-knowledge/results/procgen_family27_20260926_orchestrator_check.out
output_sha256: 19c8ea0f09fd33fe9e693fbbdf547145236ae132b62b7f0d2c50e58439783365
hash_basis: raw bytes
audit: >
  The following records the incoming orchestrator audit before the scoped
  2026-09-26 corrections documented in the body and correction audit.
  The orchestrator read Proposition M, Theorems R, R', S, G, D and
  Proposition B and found them sound. Theorem R is optional stopping for
  the AM-fair martingale; R' is the affine sandwich
  M_j <= T^j n/n < M_j + (3/4)^k; G and D are the Terras bijection plus
  the cycle lemma.
  Independent code (procgen_family27_20260926_orchestrator_check.py,
  written without reading the lane's scripts) confirms:
  * the Moran numerics (min g, lambda*, beta = 0.0239937, 1/beta = 41.677648);
  * Theorem R's identity and bounds exactly for k <= 14;
  * Theorem R' on all integers of blocks k = 10, 14, 18;
  * Theorem G exactly for k = 12, 16, 20;
  * 27's branch density 0.3929 +- 0.002 by sampling [2^24, 2^25), with a
    1/3 split over classes mod 3;
  * 27 maximises ln t(n)/ln n below 10^6.
  The OEIS path record 10709980568908647 (A006884, index 77) has
  t(n) > n^2, verified. That Kontorovich-Lagarias's Table 3 omits it
  rests on the lane's reading of their paper.
  The lane's full pipeline was re-run (519 s; 84 checks). Its output is
  identical up to timing lines. The OEIS b-files were fetched at run time
  with the generic user agent.
---

# THM-4504 — Moran pressure and the scoped families of 27

**Correction lineage, 2026-09-26.** The incoming proof and scan are retained,
with the finite-horizon, topology, scale and model qualifications below.
These repairs do not retract the fair-word martingale argument or its
finite-window integer transfer. See the [full note](../../05-knowledge/results/procgen_family27_20260926_long_orbit_families.md)
and [independent correction audit](../../05-knowledge/results/crossroads_poset_20260926_moran_audit.md).

## 1. Exact identities and their domains

For the shortcut map, write M_j=3^(o_j)/2^j. The averaged branching function
is g(s)=2^(-s)+(1/3)(3/2)^s=phi(s-1), where
phi(t)=((3/2)^t+(1/2)^t)/2. Its roots g(s)=1 are s=1,2;
its minimum is 2^(-(1-h)), h=h(log_3 2). The tangent optimization is
max_(s>0) -ln g(s)/s=0.0239937..., whose reciprocal 41.677648... is a
**model constant**, not a proved upper bound on actual stopping times.
The legality weight 1/3 is a model average, not the conditional law of
branches in a prescribed integer's inverse tree.

For fair independent parity letters, M_j is a martingale. For W>1 let
 tau=min{j:M_j>=W}, eps_k=E[M_k;tau>k]. For every finite k,

    E[M_tau;tau<=k]=1-eps_k.

Only when P(tau<=k)>0 may one divide by its conditional mean; then
W<=E[M_tau|tau<=k]<3W/2. When no hit is possible, eps_k=1 and the
conditional mean is undefined. In the infinite fair-word law,

    2/(3W)<P(W):=P(exists j:M_j>=W)<=1/W.

The hitting set is open. It agrees almost surely, but not pointwise,
with {sup_j M_j>=W}. The infinite-riser set {sup_j M_j=infinity}
has Haar measure zero and Hausdorff dimension h. No ordinary-integer
exclusion follows from this Haar-null statement.

## 2. What transfers to integers

Theorem R' gives the exact dyadic-block sandwich for j<=k using the
carry error (3/4)^k. Hence the **window** densities tend to P(W).
For complete trajectories only the lower natural-density bound

    lower_density{n:t(n)>=Wn}>=P(W)>2/(3W)

is proved here; the matching upper bound and density equality are OPEN.
Theorem S's window rise exponent is 2-beta for 1<beta<=3log_2(3)/4,
and h(beta/log_2 3) afterwards, with beta<log_2 3 throughout.

Theorem G gives, for m=L-1<=k log_3 2,

    #{n in [2^k,2^(k+1)):glide(n)>m}=2^(k-m) W_m.

W_m is the number of length-m cylinders meeting the positive-slope set Bad,
and log_2 W_m/m tends to h. Thus if m~ck, the integer-height counting
exponent is **1-c(1-h)**, while the prefix-cover exponent is h.
For c=log_3 2 these exponents differ. The separate c=1 conclusion needs
the cited THM-4495 argument, outside Theorem G's exact window.

Proposition B is conditional: for a backward-closed set with only finitely
many missing forward-image roots, a regularly varying count of index s>0,
and limiting class-2-mod-3 fraction kappa, one has
2^(-s)+kappa(3/2)^s=1. It does not establish those hypotheses for a tree.

## 3. What remains a model or finite observation

| statement | retained status |
|---|---|
| exponent 1 for backward-tree growth from the averaged function | cited model; actual tree regular variation is not proved |
| fair-word hitting probability approximately 0.83/W | numerical/model asymptotic; actual full-orbit equality is open |
| path-record count approximately 2 ln X | proved in the stated independent Frechet model; empirical for Collatz |
| delay constant 41.677648 | cited stochastic-model constant; no actual stopping-time bound asserted |
| 27-branch fraction approximately 0.3927 | incoming finite block counts, not an established natural density |
| 104 of 148 scanned delay records in 27's branch | incoming finite exact census |
| 27 maximizes ln t(n)/ln n through 2^32 and among the listed records | incoming finite exact universe only |

The inverse-tree ladder identities are exact set identities. Their averaged
weights do not determine densities of individual trees. The numerical and
heuristic claim that small members explain a tree's density is not promoted
to a density theorem. Collatz remains OPEN.
