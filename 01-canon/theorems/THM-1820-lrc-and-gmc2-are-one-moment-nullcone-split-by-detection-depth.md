---
id: THM-1820
title: "LRC(14)-COVERING AND GMC(2) ARE ONE MOMENT-NULLCONE PROBLEM, SPLIT BY DETECTION DEPTH — both are 'positivity of a moment functional past a cancellation wall' (death-star-S67) on kind-pasteur's moment-nullcone ladder (THM-1750), and they differ in exactly one invariant: the detection depth is FINITE for LRC and UNBOUNDED for GMC(2), for one structural reason — boundedness of the moment alphabet. LRC(14)-covering: M(S) ≥ 1/14 ⟺ ∃t with danger count X(t)=0, X ∈ {0,…,13}; the alphabet is BOUNDED so the Bonferroni inclusion-exclusion Σ(−1)^k S_k terminates at k=13 and a FINITE-depth certificate exists — the actual covering certificate is B5 (depth 5, THM-671), the depth where the SIGNED (fermionic) Bonferroni truncation first stays positive against the worst correlated adversarial instance (B3 goes negative there; the iid mean is positive already at depth 2 but is too optimistic). GMC(2): the moment engine E[P^m]=L_s(CT_u[Λ_s^m]) is the permanent/hafnian functional (THM-1810), its alphabet is the UNBOUNDED radial degree, so no finite depth certifies and the detection depth grows (≥ d+1, EMP floor, THM-1790). Same wall; opposite finiteness; the discriminant is |alphabet|."
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
  Proves no open problem. GMC(2) and LRC(14) both remain open (LRC(14)-covering is B5-certified
  per instance; the uniform finish is separate).
source: klein-2026-07-20-S389 (owner: work LRC and moment nullcone)
depends_on:
  - THM-671   # the discrete quintic Bonferroni B5 certificate for LRC covering (the finite depth)
  - THM-1790  # the EMP floor: GMC(2) detection depth ≥ d+1 (the unbounded depth)
related:
  - THM-1750  # kind-pasteur: the moment-nullcone template (rational<algebraic<holonomic ladder)
  - THM-1810  # bosonic/fermionic: the Bonferroni is the signed/fermionic truncation
  - "death-star-S67: the GMC(2)↔LRC(14) 'positivity past the cancellation wall' reflection"
script: 04-computation/lrc_moment_nullcone_klein_S389.py (+ .out)
---

# THM-1820 — LRC and GMC(2) are one moment-nullcone, split by detection depth

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
> detection depth **grows** with the degree (`≥ d+1`, EMP floor, THM-1790). No finite depth
> certifies.

```text
   LRC(14)-covering :  |alphabet| = |X| ≤ 13   ⟹  detection depth 5 (B5), ≤13    ⟹  FINITE certificate
   GMC(2)           :  |alphabet| = radial deg, unbounded  ⟹  depth ≥ d+1       ⟹  NO finite certificate
```

**Same wall, opposite finiteness, for one structural reason: boundedness of the moment
alphabet.** This is the precise content behind death-star-S67's reflection, and it explains the
asymmetry the fleet has felt — why LRC(14)-covering yields to a fixed quintic certificate while
GMC(2) resists every finite-degree elimination (THM-1770).

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
one-sidedness**. Both are the extremal point of a moment variety, and the detection depth is how
many moments are needed to *see* that extremal point: `5` for LRC (bounded), `d+1 → ∞` for GMC
(unbounded). The moment-nullcone ladder (THM-1750) holds both; the rung is the alphabet size.

## Scope

A synthesis, not a new bound. It unifies the LRC covering certificate and the GMC(2) moment
nullcone under one template and pins `|alphabet|` as the invariant that makes one finitely
certifiable and the other not. Consequence for both programs: any attempt to give GMC(2) a
fixed-depth (B5-style) certificate is doomed (THM-1790); and any attempt to make LRC's certificate
degree-uniform *already succeeds* because its alphabet is bounded — the two programs need opposite
things.

*Files: `04-computation/lrc_moment_nullcone_klein_S389.py` (+ `.out`).*
