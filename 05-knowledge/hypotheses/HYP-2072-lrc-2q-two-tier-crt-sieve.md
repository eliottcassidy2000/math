---
id: HYP-2072
status: METHOD + validated prototype — per-component savings measured; end-to-end feasibility for n=14 needs pipeline integration
source: opus-2026-06-02-S561
related:
  - HYP-2063
  - HYP-2056
  - THM-369
---

# HYP-2072: a two-tier CRT sieve for k+1=2q that scales with the residual, not B_k

Turns the S559/S560 structural ease (n=14 per-tuple EASIER than n=10; difficulty is
pure scale) into an algorithm whose work tracks the shrinking residual.

## Method
Base correction `s·v+r·(1,…,2q-1)∈{1,…,2q-2}` (units mod 2q):
- **Tier 1 (mod-q oracle):** mod-2 forced ⇒ decide via `s'w_i+r'c_i≠f_i` mod q
  (`f_i∈{0,q-1}` by parity). Modulus q, not 2q.
- **Tier 2 (residual only):** the characterisation (apex-stuck ∨ ratio-cover) lets
  us GENERATE the hard set directly; peel apex-stuck (apex base-safe at 1/2).

## Validated (`lrc_2q_crt_sieve_prototype_s561.py`)
- tier-1 mod-q oracle == brute mod-2q: **EXACT, 0 mismatches**, q=3,5,7,11,13.
- direct residual generation == brute-found residual: **identical** (q=3 full enum,
  4465=4465).
- tier-1 coverage GROWS with q: 43%,81%,**91%(n=14)**,96%,96% — the structural ease
  is the speedup.
- apex-stuck fraction of residual grows: 30,52,**76(n=14)**,96,98% ⇒ genuinely-hard
  ratio-cover residual ~**2% at n=14**, ~0.2% at n=22.
- per cheap tuple ~2.2× lighter (modulus q vs 2q).

## Net for n=14
~91% settled cheaply; ~9% to lift, ~3/4 apex-peelable ⇒ ~2% truly hard (~45× fewer);
that 2% generated directly. Every factor improves with q — opposite scaling to the
brute c=(k+1) lift.

## Honest scope / open
Fractions over random mod-2q tuples (proxy); oracle exactness + direct generation +
apex peel are structural. Accelerates the base tier (the paper's I(k,p,1)
bottleneck); end-to-end n=14 feasibility needs integration with the Rosenfeld-
Trakulthongchai higher-c lifts. Not a proof of LRC(14). Next: (1) wire mod-q oracle
+ generation into the pipeline, remeasure on real workload; (2) full analytic apex
peel; (3) clear the ratio-cover residual via r/p higher-lift freedom.

**See:** `07-reflections/lrc-2q-two-tier-crt-sieve-exploiting-structural-ease-s561.md`,
`04-computation/lrc_2q_crt_sieve_prototype_s561.py` (+.out); HYP-2063 (apex/2q),
HYP-2056 (even-fold), THM-369; arXiv:2604.23906.
