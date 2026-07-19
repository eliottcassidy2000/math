---
id: THM-1013
title: Dilated empty-circle sieve
status: PROVED and kernel-formalized under its explicit distance-from-ndZ hypothesis; historical claims that it covers every dilated-AP-core compact family are CORRECTED
source: boxeph-2026-07-18-S86, scope corrected by codex-2026-07-18-S74
related: [THM-1099, THM-1149, MISTAKE-170]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCDilatedSieve.lean
---

# THM-1013 -- the dilated empty-circle sieve

> **Theorem.**  Fix positive integers `n,d`.  If every speed `v_i` satisfies
>
> ```text
> |v_i-ndm|>=d for every integer m,                       (1)
> ```
>
> then `t=1/(nd)` is `1/n`-lonely.

**Proof.**  For each speed, choose the nearest integer `m` to `v_i/(nd)`.
Then

```text
||v_i/(nd)||=|v_i-ndm|/(nd)>=d/(nd)=1/n.                (2)
```

∎

The theorem is formalized by `LonelyRunner.dilated_sieve` and
`dilated_sieve_lonely14` in `LRCDilatedSieve.lean`; the public axiom footprint
is `[propext, Classical.choice, Quot.sound]`, with no `sorry`.

For `n=13`, (1) says every speed lies at distance at least `d` from
`13d Z`, and (2) gives `M>=1/13`.

## Correct scope for an AP deletion core

Every member of the core

```text
d{1,...,12}
```

satisfies (1), but an additional speed need not.  Therefore the implication

```text
“contains a dilated twelve-term AP deletion”
  => “THM-1013 supplies a 1/13 witness”
```

is false without a separate condition on the extra speed.  The elementary
example `[12] union {13}` already shows the issue: the extra speed lies on
`13Z`, and the full set has value `1/14`, not `1/13`.

The historical S86/S113 prose promoted the conditional sieve to a statement
about all dilated-AP-core compact rows and then used that promotion in an
alleged LRC14 equivalence.  Those claims are withdrawn; the kernel theorem
itself is unchanged.

THM-1149 supplies the correct strict-cover statement:

```text
M(d[12] union {v})<1/13  =>  13d divides v.              (3)
```

Thus, after a tight AP deletion has been extracted, primitivity, one
14-carrier, and `rho<13` give a contradiction.  THM-1013 remains the direct
witness whenever the explicit distance hypothesis (1) holds; THM-1149
handles the complementary regeneration logic.  Neither theorem extracts a
tight deletion, proves the twelve-speed equality classification, or proves
LRC(14).

## Methodological lesson

The correct carrier is the distance of **every** speed from the lattice
`ndZ`.  Knowing only that a subset has a dilated-AP form forgets the extra
runner, precisely the coordinate on which the witness can fail.
