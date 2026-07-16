---
id: THM-887
title: THE GENERAL-CLUSTER RESONANCE LAW (the multi-comb assembly; THM-886's named next step executed) — for ANY cluster E, any section, the signed endpoint sum decomposes by owners S = Σ_e S_e with: (I) PER-OWNER COMB IDENTITY, exact for every owner e and every integer a: S_e(ea) = ν̂_e(a) = Σ_{c mod 7} N_c^{(e)} e(ac/7) — every owner carries a 7-periodic mode comb on eZ (phase collapse, two lines); (II) PER-OWNER OFF-COMB PROFILE |S_e(m)| ≤ min(M_e, R_e/sin(π‖m/e‖)) (run-geometric; the honest constant — the /2-form fails per-owner and THM-886(II) is corrected in place); (III) MULTI-COMB PROFILE |S(m)| ≤ Σ_e profile_e(m), refereed 0 violations on five cluster types with worst ratio 1.000 (SHARP on a balanced section); (IV) combs meet on lcm-lattices — compound resonance needs simultaneous Diophantine proximity (CRT); (V) THE ASSEMBLED LAW Q_s(w) ≤ C₁M + 2π²Σ_ℓ k̂-bound(ℓ)·[Σ_e profile_e(ℓw)]², covering-verified on every cluster × w battery (bounds 39–89·M uniformly, incl. resonant w). STRUCTURAL READING: the resonance classifier is w's SIMULTANEOUS Diophantine position against the cluster's own comb system {eZ} — Q_s(w) is governed by how lonely w is against E itself at mode scale; the family (THM-886/klein cont.4) = one dominant comb; the balanced regime (death-star S19 in flight) = all combs comparable
status: (I)–(III) PROVED (elementary; machine-refereed: (I) exact on every owner of five cluster types; (II) 0 violations at the proven constant; (III) 0 violations, worst ratios 0.524–1.000); (V) upper law PROVED with the stated constants and covering-verified (deliberately lossy, uniform); the sharp off-resonance constant and the compound-hit (CRT) optimization are the named refinements
source: boxeph-2026-07-16-S26 (owner directive: take the general-cluster resonance law)
depends_on: [THM-886 (family case + k̂ bound), THM-880/881 (frame), klein cont.4 (family mode law), death-star S18/HYP-7017 (independent refutation), death-star S19/HYP-7018 (balanced-regime limit measures, in flight — complementary)]
script: 04-computation/lrc14_general_resonance_law_boxeph_S26.py -> 05-knowledge/results/lrc14_general_resonance_law_boxeph_S26.out
---

# THM-887 — the general-cluster resonance law

**Setup.** Cluster E, section s, swing set R_s, P = 7·lcm(E). Every R_s endpoint is
owned (THM-881 P1): S(n) = Σ_{e∈E} S_e(n). For owner e: endpoints at j/(7e); residue
class c = j mod 7; N_c^{(e)} = signed count in class c; M_e = #endpoints;
R_e = total run count of the per-(class, sign) active index sets in Z_e.

**(I) Per-owner comb identity.** For every e ∈ E and every integer a:
> S_e(ea) = ν̂_e(a) := Σ_{c mod 7} N_c^{(e)} · e(ac/7).
*Proof.* At n = ea the phase of the endpoint j/(7e) is e(ea·j/(7e)) = e(aj/7), which
depends only on j mod 7; collect classes. ∎
Each owner carries a mode comb on eZ whose teeth are 7-periodic in a (at most four
distinct strengths: |ν̂_e| at a ≡ 0, ±1, ±2, ±3 mod 7); the a ≡ 0 tooth is Σ_c N_c^{(e)}.

**(II) Per-owner off-comb profile.** For m ∉ eZ:
> |S_e(m)| ≤ min(M_e, R_e / sin(π‖m/e‖)).
*Proof.* Group by class: S_e(m) = Σ_c e(mc/(7e))·Σ_{k∈A_{c,±}} ±e(mk/e); each run of
consecutive indices in Z_e is a geometric block with |1−ζ^L|/|1−ζ| ≤ 2/|1−ζ| =
1/sin(π‖m/e‖), ζ = e(m/e); R_e runs total; the M_e cap is trivial. ∎
(The /2-strengthening FAILS per-owner — 278/1488/8624/4447 violations across the
battery — while the honest constant gives 0 violations; THM-886(II) corrected in place.)

**(III) Multi-comb profile.** |S(m)| ≤ Σ_e profile_e(m) (triangle inequality), where
profile_e = |ν̂_e| on the comb, (II) off it. Refereed over full Z_P on: the family
{1..6,60}; two-owner {1..5,36,60}; two-large-sharing-a-factor {1..5,56,84}; the
BALANCED compact-core-like {12,15,20,21,28,30,35} (all seven sections); near-AP
{8,9,10,12,14,15,18}: **0 violations everywhere, worst ratios 0.524–1.000** — sharp
at a balanced section: the profile is the right envelope, not an artifact.

**(IV) Comb transversality.** eZ ∩ e'Z = lcm(e,e')Z; within any window of N
frequencies the double-comb points number N/lcm(e,e') + O(1). Compound resonance of w
against two owners requires SIMULTANEOUS Diophantine proximity ‖ℓw/e‖, ‖ℓw/e'‖ small —
CRT makes this rare exactly to the extent that gcd(e, e') is small. (Near-AP clusters:
pairwise gcds small ⟹ compound hits negligible; gcd-rich clusters concentrate them.)

**(V) The assembled law.** With |k̂_P(n)| ≤ 1/(2P²sin²(πn/P)) (THM-886(III)):
> Q_s(w) ≤ C₁·M + 2π² Σ_{ℓ≥1} k̂-bound(ℓ) · [Σ_e profile_e(ℓw mod P)]²,
C₁ = 12 suffices on the battery (low-band absorber). Covering-verified on every
cluster × w battery (random / per-owner resonant / gcd-sharing / compound): bounds
land at 39–89 × M uniformly — the law is O(M)-scale off-resonance and grows exactly
through the profile at resonant w.

## The structural reading

The classifier of (V) is the vector of Diophantine positions
(‖ℓw/e‖)_{e ∈ E, ℓ small}: **w is resonant against the cluster exactly as a runner
is lonely-or-caught against a torsion system — Q_s(w), the second moment of the
good-set geometry at dilation w, is governed by how lonely w is against E itself at
mode scale.** The LRC structure recurses into its own second moment. Special cases:
one dominant comb = THM-886/klein-cont.4's family law; all combs comparable = the
balanced regime whose limiting tooth strengths death-star S19 is computing
(their A* = 3B* limit measures are the ν̂-amplitudes' balanced limits); every
regime in between is a comb-strength profile.

## Consequences

1. The per-instance decidability (THM-881 P2) + this profile give CHEAP a-priori
   Q_s(w) certificates: no scan needed, just the comb bookkeeping (M_e, N^{(e)}, R_e).
2. The peel-sequence limit w → ∞ stays safe for every cluster (profiles are bounded
   by Σ M_e = M; the k̂-weights kill the tail) — the density route's two-scale error
   closes with EXPLICIT per-cluster constants: the constant-propagation seam (kps
   task (b)) can consume (V) directly.
3. The named refinements: (a) the sharp off-resonance constant (the battery says ~1–6
   vs the law's 39–89 — a factor ~10–40 is on the table via the exact k̂ instead of
   its bound and the CRT compound-hit count); (b) the balanced-regime tooth-strength
   asymptotics (death-star S19's lane — complementary, not duplicated).

## Evidence log
- [x] (I) exact, every owner, five cluster types; (II) 0 violations at the proven
      constant (and the /2-form's failure censused); (III) 0 violations, sharp at 1.000
- [x] (V) covering on all cluster × w batteries
- [ ] sharp off-resonance constant; CRT compound-hit optimization
- [ ] plug (V) into the THM-727/729 row-closure chain (task (b) consumer)
