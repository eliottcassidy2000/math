---
id: THM-568
title: Apex-shell divisibility lemma — every tight LRC(14) crossing has 14 | D and D | active pair sum
status: PROVED (elementary) for the local apex-shell divisibility. CORRECTED by codex-S120/HYP-2929: the stronger primitive conclusion D=14 does not follow from the displayed arithmetic alone; shell-collapse remains a separate proof target.
author: kind-pasteur-2026-06-22-S31aa
depends_on:
  - THM-523    # q=14 witness (14-free ⟹ M_14 >= 1/14)
  - S31v       # comb-teeth union bound + M(R)>=1/13 margin
related:
  - THM-079    # the H=21 template (Move A reduce-to-atom + Move B forcing)
  - THM-560    # tight locus / dilation
  - HYP-2906   # mac-mini's (★)
  - HYP-2929   # correction: shell-collapse is not proved by this local lemma
---

# THM-568: the apex-shell divisibility lemma

## Statement
Let `S` be 13 distinct positive integers with `M(S) = max_t min_i ‖v_i t‖ = 1/14` (a tight set), and let
`t* = a/D` (lowest terms) be an optimum. Then:
1. **`14 | D`.**
2. The two binding runners `{v_a, v_b}` (those with `‖v_· t*‖ = 1/14`, at `±1/14`) satisfy `D | (v_a+v_b)`.

Equivalently, every tight crossing lies on an apex shell `D=14h`, and the
active pair sum is divisible by that shell.  The additional shell-collapse
claim `h=1` for primitive tight rows is not part of this theorem.

## Proof
At a tight optimum the origin lies in an empty open arc of length exactly `1/7`, with binding runners at
**both** endpoints `±1/14` (a single occupied endpoint would let the arc — hence `M` — grow, contradicting
optimality). For the `+1/14` binder, `v_a·(a/D) = m ± 1/14` for some `m ∈ ℤ`, so `14(v_a a − mD) = ±D`,
i.e. `D(1+14m) = 14 v_a a`. Since `1+14m ≡ 1 (mod 14)` is coprime to 14, `(1+14m) | v_a a` and
`D = 14·(v_a a/(1+14m))`, so **`14 | D`**. The two binders give `(v_a+v_b)t* ∈ ℤ`, so `D | (v_a+v_b)a`,
hence `D | (v_a+v_b)` because `(a,D)=1`. ∎

VERIFIED `lrc_apex_denominator_lemma_kps.py`: AP (`D=14`, binders `{1,13}`), GW (`D=14`, binders `{5,9}`),
`2·AP` (`D=28`), `3·GW` (`D=42`).

FORMALIZED `TournamentH7.LRCApexShell` (codex-S120): the integer equations
`14*(u*a-m*D)=D` and `14*(v*a-n*D)=-D` imply `14|D`, `D|(u+v)`, and
`14|(u+v)` with no extra assumptions beyond `IsCoprime D a`.

## Consequence — the reduction of (★)
`(★)` [`M(S)=1/14 ⟹ optimum at a denom-14 point`] reduces, after this local
apex-shell lemma, to two still-distinct parts (`S` 14-covering iff it contains a
multiple of 14):
- **14-free tight ⟹ a denom-14 optimum exists** (no equidistribution): THM-523 gives `M_14(S) ≥ 1/14`, so
  tightness forces `M_14 = 1/14` attained at `t=a/14`. This does not say every optimum has denominator `14`.
- **14-covering ⟹ `M(S) > 1/14`** (not tight): the 14-free part `R` has `M(R) ≥ 1/13 > 1/14` by proven
  LRC(≤13), giving an interval `I` with `min_R > 1/14`; the multiples of 14 cover `≤ |M14|/7` of `I`, so for
  `|M14| ≤ 6` a point survives with `M(S) > 1/14`. Residual: `|M14| ≥ 7` (the apex-localized second-moment).

So with Move A (peel, R1) and Move B (apex floor), THM-568 turns the THM-079-template proof of LRC(14)
into an apex-shell residual:

```text
primitive tight shell D=14h  =>  h=1
```

or, equivalently for the covering branch,

```text
>=7 multiples of 14 over a 14-free core  =>  not tight.
```

## Significance
The local arithmetic half is now PROVED and Lean-formalized: a tight set must
bind at `±1/14`, which pins `14 | D` and `D |` the active pair sum.  The
structural half of `(★)` still needs shell collapse or covering strictness.
KPS-S31ab gives supporting perturbation evidence (AP/GW-to-14-multiple swaps
have minimum `M=1/13`, no primitive tight row with `D!=14` found), but that is a
finite audit rather than the general `R union 14A` proof.

→ THM-523, THM-560, THM-079, S31v, HYP-2906, HYP-2929, `LRCApexShell.lean`,
`lrc_14covering_not_tight_kps.py`.
