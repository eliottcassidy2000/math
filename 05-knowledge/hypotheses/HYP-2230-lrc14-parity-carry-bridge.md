---
id: HYP-2230
status: OPEN bridge hypothesis with exact congruence identities and S654 bounded evidence
source: codex-2026-06-05-S654
related:
  - HYP-2229
  - HYP-2218
  - HYP-2217
  - HYP-2167
  - HYP-2166
  - HYP-2165
  - HYP-2164
  - HYP-2049
  - THM-401
---

# HYP-2230: LRC14 Parity and Apex Obstructions Are One Carry Coordinate

## Claim

The even/odd number-theory carrier and the `n=14` LRC lift/CRT seam are the
same object when written over the `C=2n-1=27` residue tower.

For a lifted speed

```text
v_i = r_i + 27 k_i,
```

the carry coordinate has two exact projections:

```text
v_i mod 2  = r_i + k_i mod 2,
v_i mod 14 = r_i - k_i mod 14.
```

So `k_i` toggles the even/odd parity word and also decides the zero-divisor
condition:

```text
14 | v_i  iff  k_i == r_i mod 14.
```

Pair-sum pinches see the same carrier at pair level:

```text
14 | v_i + v_j  iff  k_i+k_j == r_i+r_j mod 14.
```

This turns the old "counterexamples must contain a multiple of 14" statement
into a carry-cocycle constraint over HYP-2166's quotient tower.

## Connection to Even/Odd Pair Work

HYP-2218 made the user's even/odd pair idea precise:

```text
E = p+q,
O = p+2q,
q = O-E,
p = 2E-O.
```

Even targets carry an unordered pair shadow, while odd targets carry an
ordered/doubled-prime bridge.  The duplicate diagonal sends

```text
p -> (2p, 3p),
```

with `p=7` giving `(14,21)`.

In LRC `n=14`, the doubled-prime bridge is not only a metaphor.  The required
apex obstruction is exactly a doubled object: the minimal carry making
`r+27k` divisible by `14` is `k == r mod 14`, and the resulting speed is even.
For AP residues `1..13`, the minimal apex speeds are

```text
28, 56, 84, ..., 364.
```

For `V*=(1..11,13,24)`, the same list appears except the residue `24` has
minimal carry `10`, producing `294=14*21`.  Thus the parity bridge, the
duplicate `(14,21)` shadow, and the LRC apex obstruction all live in the carry
fiber.

## Connection to Basel/Period Carrier

HYP-2229 says the Basel family should be read before scalar collapse:

```text
disjoint elementary packets <-> power-sum moments <-> Bernoulli local data.
```

The LRC14 analogue is:

```text
Res_27 shell packets <-> n-clock residues <-> carry/owner local data.
```

The scalar residue row `r` is like the zeta moment: useful, but incomplete.
The product/carry side channel is the object that keeps the equality honest.
Just as Newton identities recover power sums only after retaining elementary
packets, the n-clock obstruction is recovered from the `Res_27` row only after
retaining carries `k_i`.

## S654 Evidence

`04-computation/lrc14_parity_carry_bridge_s654.py` uses the verified S611
exact maximin oracle and records three checks.

### Single-Apex Carry Sweep

For each AP and `V*` coordinate, set exactly one carry to its minimal
apex value `k_i == r_i mod 14`, leaving all other carries zero.  This reaches
the high-carry residues `7..13` that S612's `L1<=6` search could not touch.

All `26` single-apex lifts are strict:

```text
AP best single-apex tax: M=28/365, margin=27/5110, via r=13 -> k=13.
V* best single-apex tax: M=2/25,  margin=3/350,  via r=13 -> k=13.
```

### Apex Plus Local Parity Toggles

S654 then adds all one- and two-coordinate `+27` toggles around every minimal
apex lift over AP and `V*`.

```text
rows screened: 2054
approximate near-floor rows exact-checked: 0
exact floor rows in checked set: 0
exact below-floor rows in checked set: 0
```

Every exact group minimum by base, apex residue, and extra-toggle count is
above `1/14`.

### Boolean Carry Lattice

The full Boolean carry lattice `{0,1}^13` is also profiled by parity data:

```text
AP rows=8192, groups=98, near-floor exact-checked=1, floor/below=1/0.
V* rows=8192, groups=98, near-floor exact-checked=1, floor/below=1/0.
```

In both cases the only exact floor row in the near-floor Boolean screen is the
zero-carry base row, matching S612 and strengthening the "no first carry"
picture.

## Proof Payoff

The no-multiple branch is already discharged by the `t=1/14` clock.  Therefore
any LRC14 counterexample must satisfy at least one carry congruence

```text
k_i == r_i mod 14.
```

S654 shows that the first forced-apex bridge over the two primitive floor
shadows is not dangerous by itself: it creates strict slack.  A new primitive
floor lift would need one of the following:

1. a larger globally coherent carry cocycle, not merely the first apex carry;
2. an owner/Cprime route from HYP-2165/HYP-2167;
3. a pair-sum carry alignment whose congruence survives the endpoint-owner
   certificates.

This does not prove LRC14, but it sharpens the next theorem target:

```text
apex carry tax theorem:
  minimal apex carries over AP/V* are strict,
  and any floor-preserving carry cocycle must be globally scalar/owner-routed.
```

## Tournament Analysis

The S654 proof-channel tournament uses proof channels as vertices:

```text
no-multiple t=1/14 exit
carry parity word
apex congruence k=r mod 14
pair-sum congruence k_i+k_j=r_i+r_j
minimal-apex tax sweep
owner/Cprime certificate reattachment
raw Res_27 row identity
```

The induced order is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
hamiltonian_paths=1
```

The chosen vertices are not runners.  They are carry obligations, because this
quotient preserves the exact predicate connecting parity, apex multiples, and
pair-sum denominators.  Raw runner vertices or raw `Res_27` rows discard that
side channel.

## Honest Status

The congruence identities are exact.  The computational evidence is bounded.
It extends S612 in the important high-carry direction for the AP/`V*` floor
shadows, but it does not cover arbitrary large coherent carries or arbitrary
least-positive strict shadows.

**See:** `04-computation/lrc14_parity_carry_bridge_s654.py`;
`05-knowledge/results/lrc14_parity_carry_bridge_s654.out`;
`07-reflections/lrc14-parity-carry-bridge-s654.md`; HYP-2167, HYP-2166,
HYP-2165, HYP-2164, HYP-2218, HYP-2229, THM-401.
