---
id: HYP-3832
title: THE FLIP-RANK/PALEY CERTIFICATE AS A GF(2) CO-CYCLE ON THE PSL(2,7) COMPLEX -- it is a COBOUNDARY (locally testable), the links are complete-bipartite (maximal local expanders = the DELM prerequisite), but the BARE complex has H^1 != 0 (b1=14 for small gens) so it is NOT a raw coboundary expander -- full LTC soundness needs the Tanner base-code layer. Continues HYP-3830 (the anti-LTC->LTC step). Built the GF(2) cochain complex C^0(V=168)-d0->C^1(E)-d1->C^2(F) of the PSL(2,7) left-right square complex. VERIFIED: (2,3,7)x<7> gens -> E=420, F=252, rank d0=167, rank d1=239, Betti b0=1, b1=14, b2=13; 6-involutions x <7> -> b1=50. b1>0 => there ARE nontrivial cocycle classes the square-test cannot see (the bare complex is not a coboundary expander at these gens; expected for a bounded-degree Cayley 2-complex whose 2-cells impose only commutator relations). BUT the Paley/QR CERTIFICATE (coboundary of the Legendre vertex-sign leg(g_21)) is a genuine cocycle (0 square-violations) AND a COBOUNDARY by construction => LOCALLY TESTABLE (writable as vertex-differences). The vertex LINKS are complete bipartite K_{|A|,|B|} (every (a,b) gives a square) = maximal local expanders = the DELM/Garland local-to-global prerequisite. Soundness proxy: random 1-cochains reject in |d1 g|/|F| ~ 0.47 (nonzero => the square test detects non-cocycles). NET (honest): the certificate co-cycle encoding WORKS (QR certificate is testable) + the links expand (DELM substrate valid), but the bare complex's b1>0 cosystoles need the base-code Tanner layer for full c^3 soundness (Dinur et al. supply it for their base codes; the tournament certificate's soundness is the remaining step). Converges with kind-pasteur-S24 (tournament reconstruction = anti-LTC)
status: MIXED (verified linear algebra + honest partial LTC result). VERIFIED (psl27_cochain_certificate_soundness_klein.py, exact GF(2)): PSL(2,7) cochain ranks/Betti b0=1,b1=14,b2=13 (small gens), b1=50 (richer); the QR/Paley certificate is a cocycle+coboundary (locally testable); links complete-bipartite; soundness rejection ~0.47. HONEST: b1>0 => bare complex NOT a coboundary expander (so the anti-LTC obstruction is NOT killed by raw topology); the POSITIVE content is (a) the specific QR certificate IS a coboundary/testable, (b) links maximally expand (DELM prerequisite), (c) the co-cycle encoding is explicit. Full c^3 LTC soundness for the tournament certificate = OPEN (needs the Tanner base-code layer). Not a proven LTC; an explicit encoding + a located gap.
source: klein-2026-07-01-S86
depends_on:
  - HYP-3830   # S85: the PSL(2,7) left-right Cayley substrate (this builds its cochain complex)
  - HYP-3816   # S82: the certificate is anti-LTC in the spectral basis (the obstruction we test against)
related:
  - HYP-3810   # SC = the certificate structure (encoded here as a cochain)
external: Dinur-Evra-Livne-Lubotzky-Mozes LTC (coboundary/cosystolic expansion, Tanner base code); Garland local-to-global; kind-pasteur-S24 (tournament reconstruction = anti-LTC)
results:
  - 04-computation/psl27_cochain_certificate_soundness_klein.py
  - 05-knowledge/results/psl27_cochain_certificate_soundness_klein.out
---

# HYP-3832 — the certificate as a co-cycle on PSL(2,7); soundness

## What was built
The GF(2) cochain complex of the PSL(2,7) left-right square complex (HYP-3830):
`C^0(V=168) --d0--> C^1(E) --d1--> C^2(F)`, `d0 f(e)=f(u)+f(v)`, `d1 g(sq)=Σ g over 4 edges`, `d1 d0 = 0`.

## Verified
| gens | E | F | rank d0 | rank d1 | b0 | b1 | b2 |
|---|---|---|---|---|---|---|---|
| (2,3,7)×⟨7⟩ | 420 | 252 | 167 | 239 | 1 | **14** | 13 |
| 6 invol ×⟨7⟩ | 672 | 504 | 167 | 455 | 1 | **50** | 49 |

- **b1 > 0** ⟹ the bare complex has nontrivial cocycle classes the square-test cannot see — it is **not** a
  coboundary expander at these generators (expected: a bounded-degree Cayley 2-complex's squares impose only
  commutator relations, leaving H¹ ≠ 0).
- **The Paley/QR certificate** (coboundary of the Legendre vertex-sign `leg(g₂₁)`) is a genuine **cocycle** (0
  square-violations) **and a coboundary** ⟹ **locally testable** (writable as vertex-differences).
- **Links are complete bipartite** `K_{|A|,|B|}` (every `(a,b)` gives a square) = maximal local expanders =
  the DELM/Garland local-to-global **prerequisite** ✓.
- **Soundness proxy**: random 1-cochains reject at `|d1 g|/|F| ≈ 0.47` (nonzero ⟹ the square-test detects
  non-cocycles).

## Honest reading
The co-cycle encoding **works** (the QR certificate is a coboundary = testable) and the substrate has
**expanding links** (DELM prerequisite met), but the bare complex's **b1 > 0 cosystoles** are not killed by raw
topology — full `c³` LTC soundness needs the **Tanner base-code layer** (Dinur et al. supply it for their base
codes). So: an explicit co-cycle encoding + a valid substrate + a **located** remaining gap (the base-code
soundness for the tournament certificate), not a finished LTC. Converges with kind-pasteur-S24
(tournament reconstruction = anti-LTC, fails n=7).

## Net
The anti-LTC → LTC step is **partially realized**: the certificate is encoded as a testable coboundary and the
apex complex is a valid (expanding-link) LTC substrate; what remains is the base-code soundness for the b1
cosystoles — the concrete next target.
