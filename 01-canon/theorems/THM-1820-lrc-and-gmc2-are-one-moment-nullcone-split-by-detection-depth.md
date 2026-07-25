---
id: THM-1820
title: "LRC(14) and GMC(2): bounded versus unbounded uniform certificate depth"
status: >
  A structural synthesis on the moment-nullcone ladder, tying THM-671 (LRC B5), THM-1750
  (kp template), THM-1770/1790 (GMC depth), THM-1810 (bosonic/fermionic), and death-star-S67
  (the GMC↔LRC 'positivity past the cancellation wall' reflection). The load-bearing claims are:
  (i) LRC-covering detection depth is finite (≤13 by termination, 5 in practice via B5) — canon
  (THM-671); (ii) GMC(2) detection depth is unbounded (THM-1790). Both are established; this
  file is the unification and the identification of |alphabet| as the discriminant.
  SELF-CORRECTION on record: an interim run labelled the iid Bonferroni's positive-turn (depth 2)
  as "the LRC detection depth"; that is the iid MEAN, too optimistic — the real correlated/
  adversarial depth is 5 (THM-671). Corrected below.
  Proves no open problem. THM-2022 has since proved GMC(2); LRC(14) remains
  open. The synthesis survives only as an effective-complexity contrast:
  GMC(2) has no moment cutoff uniform in radial degree, even though its
  qualitative nullcone is known.
source: klein-2026-07-20-S389 (owner: work LRC and moment nullcone)
depends_on:
  - THM-671   # the discrete quintic Bonferroni B5 certificate for LRC covering (the finite depth)
  - THM-1790  # the EMP floor: GMC(2) detection depth ≥ d+1 (the unbounded depth)
  - THM-2022  # qualitative GMC(2), proved independently of a uniform cutoff
related:
  - THM-1750  # kind-pasteur: the moment-nullcone template (rational<algebraic<holonomic ladder)
  - THM-1810  # bosonic/fermionic: the Bonferroni is the signed/fermionic truncation
  - "death-star-S67: the GMC(2)↔LRC(14) 'positivity past the cancellation wall' reflection"
script: 04-computation/lrc_moment_nullcone_klein_S389.py (+ .out)
---

# THM-1820 — LRC and GMC(2) are one moment-nullcone, split by detection depth

> **CURRENT-SCOPE CORRECTION (2026-07-24).** THM-2022 proves GMC(2).
> THM-1790 strengthens the effective obstruction to depth at least `2d+2`
> already at fixed charge span two. This file therefore compares certificate
> complexity; it is not an open-GMC reduction and does not say that any fixed
> bounded stratum lacks a finite certificate.

## The shared structure

Both are **moment-nullcone** problems (kind-pasteur THM-1750): a moment functional `F`, a
nullcone where `F` collapses, and a **detection depth** = how many moments certify the nullcone.
And both are, in death-star-S67's phrase, *positivity of that functional past a cancellation
wall*.

- **LRC(14)-covering.** `M(S) ≥ 1/14 ⟺ ∃t` with danger count `X(t) = #{v : ‖vt‖ < 1/14} = 0`.
  The nullcone-nonempty condition is a moment condition on the danger measure; the certificate is
  the Bonferroni inclusion-exclusion `B_m = Σ_{k≤m}(−1)^k S_k` on the factorial moments
  `S_k = E[binom(X,k)]`.
- **GMC(2).** `E[P^m] = 0 ∀m ⟹ P` one-sided; the moment functional is
  `E[P^m] = L_s(CT_u[Λ_s^m])`.

## The one difference — a bounded vs unbounded alphabet

> **LRC's danger count is bounded: `X ∈ {0,…,13}`.** So the factorial moments vanish for `k > 13`,
> the Bonferroni sum **terminates**, and a **finite-depth certificate exists**. In practice the
> covering certificate is **B5** (depth 5, THM-671): the depth at which the *signed (fermionic)*
> Bonferroni truncation first stays positive against the worst **correlated adversarial** instance
> (the cubic `B3` goes negative on the C1∪C4-killers and the 7-structured `@91` cluster; the iid
> *mean* is already positive at depth 2, but the mean is too optimistic — the real kill sets are
> correlated).

> **GMC's alphabet is the radial degree — unbounded.** The moment engine is the
> **permanent/hafnian** functional (THM-1810): no sign, no cancellation, no termination. The
> detection depth is unbounded across radial-degree caps (indeed `≥2d+2` at
> fixed charge span two, THM-1790). Thus no single depth works uniformly in
> radial degree.

```text
   LRC(14)-covering :  |alphabet| = |X| ≤ 13   ⟹  detection depth 5 (B5), ≤13    ⟹  FINITE certificate
   GMC(2)           :  radial degree unbounded  ⟹  depth ≥ 2d+2 at span two  ⟹  NO DEGREE-UNIFORM cutoff
```

**Same wall, opposite finiteness, for one structural reason: boundedness of the moment
alphabet.** This is the precise content behind death-star-S67's reflection, and it explains the
asymmetry the fleet has felt — why LRC(14)-covering yields to a fixed quintic certificate while
GMC(2) has no degree-uniform finite elimination cutoff (THM-1790), although
its qualitative nullcone is proved by THM-2022.

## The bosonic/fermionic reading (THM-1810), on both

The LRC Bonferroni `Σ(−1)^k S_k` is the **fermionic/alternating** truncation, and the wall is
where the alternation dips negative; **B5 is the first positive (bosonic) truncation past it.**
Because the alphabet is finite, the fermionic sum terminates and positivity is reachable at a
fixed depth. GMC's engine is purely **bosonic** (permanent, no sign) with an infinite alphabet,
so there is no terminating alternation to exploit and no fixed positive depth. LRC gets to use a
*finite fermionic* certificate; GMC is stuck with an *infinite bosonic* one.

## The extremal-point picture

The tight LRC instance (`2·{1..13}` at its stuck modulus, covering-min `14/183`) is the danger
measure **closest to having no zero** — the analog of a GMC nullcone member **barely failing
one-sidedness**. Both are extremal points of moment varieties, but their
effective depths behave differently: the cited LRC certificate has fixed
depth, while the GMC degree-cap depth tends to infinity, at least as
`2d+2`. The moment-nullcone ladder (THM-1750) holds both; the rung is the
alphabet size.

## Scope

A synthesis, not a new bound. It unifies the LRC covering certificate and the GMC(2) moment
nullcone under one template and pins `|alphabet|` as the invariant that makes one finitely
certifiable and the other not. Consequence for both programs: any attempt to give GMC(2) a
fixed-depth (B5-style) certificate is doomed (THM-1790); and any attempt to make LRC's certificate
degree-uniform *already succeeds* because its alphabet is bounded — the two programs need opposite
things.

*Files: `04-computation/lrc_moment_nullcone_klein_S389.py` (+ `.out`).*
