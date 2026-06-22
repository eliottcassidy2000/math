---
id: HYP-2910
status: PROOF-TARGET / tightness-star template; THM-568 proves apex-denominator half, THM-569 formalizes the q=14 unit-grid split, multiples-of-14 covering-strictness remains the formal target
source: codex-2026-06-22-S119
tags: [lrc14, tight-locus, thm079-template, goddyn-wong, apex-7, covering-core, tournaments, open-q-108]
depends_on:
  - HYP-2904
  - HYP-2905
  - HYP-2906
  - HYP-2907
  - HYP-2908
  - HYP-2893
  - HYP-2895
  - THM-079
  - THM-560
  - THM-568
  - THM-569
related:
  - HYP-2909
  - OPEN-Q-108
  - THM-523
  - THM-200
  - THM-201
  - THM-343
results:
  - 04-computation/lrc14_tightness_star_template_codex_s119.py
  - 05-knowledge/results/lrc14_tightness_star_template_codex_s119.out
  - 04-computation/lrc_14covering_not_tight_kps.py
  - 05-knowledge/results/lrc_14covering_not_tight_kps.out
---

# HYP-2910: exact tightness-star atlas supporting HYP-2909

The THM-079 analogy is now a useful proof template, but it must be stated with
one extra guardrail.

For `H=21`, THM-079 first reduces to an atom by multiplicativity over strong
components, then forbids the atom by a Moon/forcing step.  The LRC14 analogue is:

```text
Move A: reduce to a primitive bounded/top-balanced atom
Move B: prove the atom has M>1/14 unless it is the non-covering AP/GW boundary
```

HYP-2905/HYP-2906 supply theorem-level pieces of Move A: omit-prime direct
witnesses, dilation normalization, and the one-large-speed interval peeler
`v_max > 13 v_second`.  HYP-2904 supplies the positive-density finite-comb
version for scale-separated and multi-large branches.

The remaining target is the tightness-star theorem:

```text
(*)  primitive M(S)=1/14
     => S is in the denominator-14 AP/Goddyn-Wong tight boundary
     => S has no multiple of 14 and is non-covering.
```

Incoming THM-568 proves the apex-denominator half of this target:

```text
M(S)=1/14 at t=a/D  =>  14 | D, D | (v_i+v_j) for the binders,
and D=14*gcd(S).
```

So a primitive tight optimum is forced to denominator `14`.  The residual is
now the 14-covering branch: if `S` contains multiples of `14`, prove
`M(S)>1/14`.  THM-568 plus S31v closes the case with at most six multiples of
`14`; the live quantitative target is the `>=7` multiples-of-14
second-moment/equidistribution branch over the 14-free core's `1/13` margin.

Codex S120/THM-569 now formalizes the exact `q=14` unit-grid split in Lean:

```text
a in (Z/14Z)^*  =>  Lonely 14 v (a/14) <=> no speed is divisible by 14.
```

So the finite apex boundary is no longer just a script checksum.  The Lean
module `TournamentH7.LRCUnitGrid14` exposes named theorems for
`a=1,3,5,9,11,13` and the specialized corollary that any no-lonely-time family
must contain a multiple of `14`.

Post-rebase KPS-S31ab is strong incoming signal for the next obligation: it
claims the 14-covering residual is not tight and verifies the claim on AP/GW
replacement families, with minimum `M=1/13`.  Read conservatively, the script
does not by itself replace the still-needed general theorem; it proposes the
formal mechanism to extract next:

```text
S = R union M14, R 14-free
  -> the 1/13-margin supplied by the smaller runner theorem cannot be covered
     by the 14-multiple danger combs.
```

That is now the sharp formalization target after THM-569.

Equivalently, in the HYP-2908 language, a remaining bounded apex-7 over-cover
must state-lift to a tournament-conflict-realizable connected binary packet
graph with `I(.,2)=7`, hence to the forbidden `K_3` atom.

## Exact evidence from S119

The S119 audit verifies:

```text
AP {1,...,13}:          M=1/14
GW {1,...,11,13,24}:   M=1/14
argmaxes in both rows: k/14, k in {1,3,5,9,11,13}
```

The binding pairs match in both rows:

```text
t=1/14,13/14: binders {1,13}
t=3/14,11/14: binders {5,9}
t=5/14,9/14:  binders {3,11}
```

The exact denominator-14 grid obstruction is theorem-level:

```text
for k in {1,3,5,9,11,13}, 14 | v*k  <=>  14 | v.
```

Thus if `S` has no multiple of `14`, every listed `k/14` survives, so
`M(S)>=1/14`.  Any strict counterexample must contain a multiple of `14`.
This exact statement is now proved in Lean as THM-569, in the stronger
predicate form
`Lonely 14 v (k/14) <-> forall i, not (14 | v_i)`.

The exact bounded AP single-swap atlas through replacement `v<=80` has only
one non-AP tight row:

```text
{1,...,13} with 12 -> 24 = GW.
```

No single-swap row in that atlas has `M<1/14`.

Finally, the finite q-covering window `[1,18]` under the necessary condition
"for every q=2..14, some speed is a multiple of q" has strict slack:

```text
966 rows checked exactly
minimum M = 1/12
minimum row = (1,2,3,4,10,11,12,13,14,15,16,17,18)
0 tight rows, 0 below-threshold rows
```

This supports Move B and is now consistent with THM-568's sharper reduction:
the finite window is an exact checksum for 14-covering slack, not the source of
the apex-denominator theorem.

## The necessary correction

The implication `M(S)=1/14 => AP/GW` alone does not logically rule out a
strict bounded counterexample with `M(S)<1/14`.  The template closes only if
one of the following is also proved:

1. a variational/compression theorem: every bounded over-cover descends to a
   tight boundary row, contradicting `(*)`; or
2. a direct state-lift theorem: every bounded apex-7 over-cover produces the
   forbidden HYP-2908 `K_3` packet; or
3. an exact finite bounded-core atlas with positive minimum margin.

This is the same discipline as THM-079.  In THM-079, the atom reduction and
the atom-forcing theorem are both required.  For LRC14, `(*)` is the atom
forcing theorem only after the bounded-core reduction carries enough boundary
state to force equality or a forbidden packet.

THM-568 improves this status: the equality-to-apex-denominator part is proved
elementarily, so the remaining atom-forcing theorem can be stated as

```text
14-covering bounded atom  =>  not tight, in fact M>1/14.
```

The `<=6` multiples-of-14 subcase is closed by the 14-free core margin and the
comb-teeth union bound; the `>=7` subcase is the apex-localized Node-3
second-moment branch.

## Tournament Analysis / assumption challenge

Candidate vertices considered:

```text
raw runners
denominator-14 survivor points
AP/GW tight rows
single-swap atlas rows
q-covering windows
bounded covering cores
apex-7 over-cover packets
generic digraph shadows
tournament OCF conflict packets
```

Chosen quotient: proof carriers ordered by how much of the LRC predicate they
preserve.  The script's carrier tournament is transitive:

```text
denom14_grid_floor
  > AP_GW_exact_tight
  > single_swap_atlas
  > q_covering_window
  > K3_state_lift_target
  > generic_digraph_shadow
```

The challenged assumption is that "tight locus AP/GW" and "K3/H7" are
automatically the same statement.  They are not.  The bridge must be a
realizability theorem from apex-7 over-cover cells to the tournament/OCF packet
category.  HYP-2908 shows that loose digraph shadows are too broad.

## Current proof order

The sharp next theorem should be stated as:

```text
Every primitive top-balanced bounded covering LRC14 atom either
  (a) compresses to the AP/GW denominator-14 boundary, hence is non-covering;
  (b) has exact slack M>1/14; or
  (c) realizes the forbidden connected K3 packet, impossible by THM-201/343.
```

After THM-568, the sharpest version is:

```text
S = R union M14, where M14 are multiples of 14 and |M14|>=7.
R is 14-free and has a 1/13-margin interval by LRC(<=13).
Prove the multiples of 14 cannot cover that interval at threshold 1/14.
```

Proving this localized equidistribution/second-moment theorem would make the
THM-079 template genuinely close LRC14.
