---
source: opus-2026-07-31-S4 (FC(2)/Kontsevich-Zagier, repo JC/Keller/period state, idea-translation map)
status: >
  SYNTHESIS + verified connection + caveats. The Factorial Conjecture FC(2) exponential-period approach the
  owner's friend describes maps DIRECTLY onto existing repo machinery, and one connection is exact: Long's
  July-2026 GM(3) counterexample (E[P^m]=0, E[QP^m]=m!) is byte-for-byte the repo's "no-pole" object
  (deathstar-S61e), since E[Q e^{tP}] = sum m! t^m/m! = 1/(1-t), a pole at t=1. Three external July-2026
  results (Alpoge JC(3) counterexample = repo THM-1300; Wilson GM(2) true = repo THM-2022; Long GM(3) false =
  repo S61e) show the repo is synchronized with the frontier. FC(2) is the live target; FC(n>=3) collapses
  with the JC. CRUCIAL CAVEAT: planar JC(2) needs GM(4), which is FALSE (Long), so the FC/GM moment route does
  NOT close JC(2) -- the repo's JC(2) progress is a separate rigidity cage. Translation: period-irreducibility
  (my S(k) work) <-> period-nonvanishing (FC) <-> no-pole/monodromy, all governed by Kontsevich-Zagier.
tags: [factorial-conjecture, FC2, jacobian-conjecture, keller, mathieu-zhao, gaussian-moments, exponential-periods, kontsevich-zagier, no-pole, furter-rigidity, sk-series, translation, repo-state]
related: [THM-2022, THM-1300, THM-1510, THM-1370, Sk-elementary-locus, central-binomial-edge]
external: ["Alpoge JC(3) counterexample (Fable)", "Wilson arXiv:2607.23887 GM(2) true", "Long arXiv:2607.18186 GM(n>=3) false", "Adamus-Hajto arXiv:2607.09340 Pascal-finite", "Edo-vdE arXiv:1304.3956 FC<->Furter", "Derksen-vdE-Zhao arXiv:1506.05192 GM=>JC"]
---

# The Factorial Conjecture, exponential periods, and the repo state: a translation map

## 1. The July-2026 landscape (external), and the repo is already on it

| Result | External (July 2026, source-flagged) | Repo | Match |
|---|---|---|---|
| **JC(3) is FALSE** | Alpoge, weighted grading `(-1,1,2)`, `det=-2`, 3-to-1 (found w/ "Fable") | THM-1300 (same map, `det=-2`, weights `(1,-1,-2)`, Dixmier A3) | **SAME object** (repo annotates it as Alpoge/Fable) |
| **GM(2) is TRUE** | Wilson arXiv:2607.23887 (face isolation) | THM-2022 (NC2/GMC(2) via de-factorialization, Lean-checked) | **SAME theorem**, two proofs |
| **GM(n>=3) is FALSE** | Long arXiv:2607.18186: `E[P^m]=0, E[QP^m]=m!` (3-var quartic) | deathstar-S61e "no-pole": `E[Z e^{tP}]=1/(1-t)`, pole at `t=1` | **SAME object** (verified below) |

The repo has been shadowing the research frontier. **Verified exactly:** Long's counterexample gives
`E[Q e^{tP}] = sum_m E[QP^m] t^m/m! = sum_m m! t^m/m! = sum_m t^m = 1/(1-t)` -- precisely the repo's no-pole
failure mode. So "the Gaussian/Factorial moment conjecture fails" `<=>` "the exponential generating function
`E[Q e^{tP}]` develops a pole at `t=1`." That pole is the whole phenomenon.

(arXiv:2607.09340, Adamus-Hajto, is adjacent but off-topic: planar automorphisms = diagonal-Jacobian scaling
composed with Pascal-finite factors; its "Exponential Generators Conjecture" is exp-of-LND generators, NOT
exponential PERIODS -- do not conflate.)

## 2. FC(2): statement, exponential-period form, and the friend's Kontsevich-Zagier route

**FC(n):** for `f in C[x_1..x_n]`, if `L(f^m)=0` for all `m>=1` where `L(x^alpha)=alpha!`, then `f=0`. Since
`int_0^inf x^a e^{-x}dx = a!`,

```
   L(g) = int_{[0,inf)^n} g(x) e^{-(x_1+...+x_n)} dx  = E[g],  x_i ~ Exp(1) iid;
   L(f^m) = int_{[0,inf)^n} f^m e^{-sum x} dx  are EXPONENTIAL PERIODS (for f in Qbar[x]).
```

`FC(2) <=> [ int_{[0,inf)^2} e^{t f} e^{-x-y} dx dy == 1 for all t  =>  f == 0 ]` (all moments `m>=1` vanish
`<=>` the exponential integral is constant in `t`). The **friend's route**: `L(f^m)` are periods; the WEAK
Kontsevich-Zagier (exponential) period conjecture says a period vanishes only for a geometric reason
(additivity / change of variables / Stokes). Classify, in dimension 2, all geometric mechanisms giving
simultaneous vanishing for every `m`, and show they force `f=0`. This is exactly the repo's **no-pole
reformulation** (S61e), from the period side: FC `<=>` no pole `<=>` no geometric obstruction.

**Two cautions (both already visible in the repo):** (i) `f` is complex-valued, so the pushforward of the
positive measure `e^{-sum x}dx` is a COMPLEX measure -- ordinary moment-determinacy fails, cancellation is
possible; this is exactly why the repo needed `p`-adic (THM-2022) not real analysis, and why K-Z (not real
moment theory) is the right tool. (ii) Any `K-Z => FC(2)` argument MUST use an `n=2`-special feature, because
`n>=3` genuinely collapses (Long, Alpoge) -- the same `n=2 | n>=3` wall that governs the JC.

## 3. THE CRUCIAL CAVEAT (do not overclaim FC(2) => planar JC)

The chain is `GM(2n) => SIC(n) => JC(n)`. **Planar JC(2) needs GM(4)** -- and GM(4) is FALSE (Long,
`n>=3`). Wilson's GM(2)=true only feeds JC(1) (trivial). So **the moment / factorial route to the planar
Jacobian conjecture is DEAD**, yet JC(2) remains open -- no contradiction. Consequently:

> FC(2) falling (via K-Z) would settle the Factorial Conjecture in the plane and (via Edo-vdE) bear on
> **Furter's Rigidity Conjecture**, but it does NOT close planar JC(2) through the moment chain. The owner's
> "implications for Keller pairs / planar Jacobian" are real but INDIRECT -- through Furter rigidity and the
> shared `n=2` geometry, not the `GM => JC` implication. The repo's actual JC(2) progress is the independent
> **rigidity cage** (THM-2063/2102/2113/1345/1370, THM-2696-2722): every accessible planar Keller stratum is a
> tame triangular automorphism; no counterexample; JC(2) framed as an Abhyankar-Moh leading-form descent.

## 4. Idea translation across the repo's frontiers (the owner's real ask)

The unifying object is the **period** and its vanishing/irreducibility, which Kontsevich-Zagier governs. This
lets ideas move between problems that look unrelated:

- **My S(k) work `<->` FC.** I proved `S(k>=4)` are IRREDUCIBLE hypergeometric-motive periods (they do NOT
  vanish into elementary pieces), and pinned the elementary locus by the quadratic-transformation signature
  `b-a in {0,+-1/2}`. FC asks the dual: whether exponential periods `L(f^m)` all VANISH. Same toolkit --
  monodromy of the period family, PSLQ-certified (non)relations, the "edge value / no-pole" analytic
  signature. My PSLQ-irreducibility certificates for `S(4)` are a template for certifying that a candidate
  FC(2) `f` does NOT have all-vanishing moments (find one `m` with `L(f^m) != 0` via an exponential-period
  computation), and the no-pole/monodromy view is shared.
- **AMM golden `<->` Paley tournament `<->` these periods.** The central-binomial/semicircle object
  (`xF^2+F-1=0`, `F(1)=1/phi`) is a period too; the AMM capacity floor and Paley path-ratio are its edge and
  exponential evaluations. FC's generating function `int e^{tf}e^{-sum x}` is the exponential evaluation of
  the FC period family -- the SAME "edge vs exponential resummation" duality I used for AMM<->Paley.
- **The `n=2 | n>=3` wall is universal.** JC (open n=2, false n>=3), GM (true n=2, false n>=3), and the S(k)
  level bound (`C_{1/4}={1,2,3}`, irreducible at 4) all have a small-parameter regime that dies at a threshold.
  This is the same "capacity holds until it doesn't" shape as the AMM golden floor. The transferable lesson:
  the FC(2)/K-Z classification of dim-2 vanishing mechanisms is a CAPACITY argument (how many geometric
  cancellations can a plane curve `f=const` support), structurally like THM-3009's `(ARCH)` and my S(k)
  motivic-dimension gate.

## 5. State of the repo's open problems, and where a fresh idea lands

Closed (repo, matching frontier): FC(1)=EMP (THM-1510); GM(2)/NC2 (THM-2022 = Wilson); positive-graded JC all
dims (THM-1370). Open + tractable, in the exponential-period habitat (repo agent's ranking, endorsed):
1. **No-pole control of `E[Q e^{tP}]` for FC(2)** (S61e) -- the cleanest K-Z hook; classify when a dim-2
   exponential integral can acquire a pole/Stokes ray. This is the friend's approach in the repo's language.
2. **Effective Strong-Factorial cutoff `(k-1)R`** (HYP-8765) -- bound the vanishing ORDER of `L(H^j)` from the
   Newton/support data of `H`; the "return-2 death-at-6" fakers are exact test data. A period-theoretic
   vanishing-order bound is the missing piece (and is where my S(k) monodromy/PSLQ methods transfer most
   directly).
3. **Confluent-Vandermonde reality wall** (THM-2033) -- a Laguerre-Polya question, natural for
   exponential-period asymptotics.

FC(3): expect false (SFC false, GM(3) false, JC(3) false), but the BARE factorial FC(3) is not directly
refuted anywhere -- the DEZ dictionary is "exponential in `n` <-> Gaussian in `2n`" only for special symmetric
polynomials, so Long's Gaussian counterexample does not mechanically transfer. A worthwhile concrete target:
**transport Long's `E[QP^m]=m!` object across the DEZ dictionary to test whether bare FC(3) is literally
false** -- a finite computation the repo is equipped for.

## Honest scope

Reconstructions flagged: the `K-Z => FC(2)` implication is a research STRATEGY (no published proof); the
July-2026 arXiv items and the Alpoge counterexample are post-cutoff and single/few-source. The exact
verification here is only the algebra `1/(1-t)` linking Long and S61e, and the repo-vs-external identifications
by matching invariants. The value is the map and the translations, not a new theorem on FC.
