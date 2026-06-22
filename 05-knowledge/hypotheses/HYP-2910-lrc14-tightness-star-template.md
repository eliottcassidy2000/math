---
id: HYP-2910
status: PROOF-TARGET / tightness-star template; THM-568 proves the local apex-shell arithmetic, THM-569 formalizes the q=14 unit-grid split, THM-571 closes the apex-majority branch, and HYP-2929 records the shell-collapse guardrail
source: codex-2026-06-22-S119
tags: [lrc14, tight-locus, thm079-template, goddyn-wong, apex-7, covering-core, tournaments, open-q-108]
depends_on:
  - HYP-2904
  - HYP-2905
  - HYP-2906
  - HYP-2907
  - HYP-2908
  - HYP-2911
  - HYP-2929
  - HYP-2893
  - HYP-2895
  - THM-079
  - THM-560
  - THM-568
  - THM-569
  - THM-570
  - THM-571
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
  - 04-computation/lrc14_apex_majority_shift_guard_codex_s121.py
  - 05-knowledge/results/lrc14_apex_majority_shift_guard_codex_s121.out
  - 04-computation/lrc14_apex_majority_gamma_descent_codex_s122.py
  - 05-knowledge/results/lrc14_apex_majority_gamma_descent_codex_s122.out
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

Incoming THM-568, corrected by HYP-2929, proves the apex-shell half of this
target:

```text
M(S)=1/14 at t=a/D  =>  14 | D and D | (v_i+v_j) for the binders.
```

So a primitive tight optimum is forced onto an apex shell `D=14h`, but the
collapse `h=1` is a separate theorem target.  The residual is now the
14-covering/shell branch: if `S` contains multiples of `14`, prove
`M(S)>1/14` for every shell height.  THM-571 closes the apex-majority
subbranch with at least seven multiples of `14`, modulo the accepted LRC<=13
input.  The remaining branch is the scale-separated / bounded-core case with
at most six multiples of `14`, where S31v supplies the comb estimate and a
finite-core compression/census theorem must supply the final uniform margin.

Incoming THM-569 now formalizes the exact `q=14` unit-grid split in Lean:

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

Codex S121/THM-570 narrowed the `>=7` multiples-of-14 branch.  For

```text
S = 14Q union R, |Q|>=7, |R|<=6,
```

the below-14 theorem gives a strict safe interval for `Q` in the scaled
coordinate `u`; all fourteen lifts `t=(u+k)/14` keep the `14Q` block safe.
The exact shift sieve proves the branch is safe if `R` has at most one speed
with `gcd(r,14)=7`: ordinary residual speeds forbid at most two shifts and at
most one per parity, while one half-step speed forbids at most one whole parity
class.  Thus the broad apex-majority residual has collapsed to HYP-2911:

```text
R contains at least two speeds divisible by 7 but not by 14,
and their u-dependent half-step phases cover both parities.
```

The guardrail is exact: at `u=2/49`, `r=7` forbids all even shifts and
`r=161` forbids all odd shifts.  So the next proof cannot be residue-only; it
must use the strict `Q` interval, actual phase labels, or an exact-period /
equidistribution estimate for half-step packets.

Codex S122/THM-571 supplies that missing move for the actual apex-majority
branch.  If two or more residual half-step speeds exist, then at least nine
speeds are multiples of `7`.  Scale those by `7`, use the proven below-14 LRC
input to choose a strict safe `7`-phase, and pigeonhole over the seven lifts.
At most four speeds are not divisible by `7`, and each forbids at most one of
the seven lifts.  Hence a lift survives.  Therefore every primitive 13-speed
row with `|M14|>=7` is LRC14-safe.

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

## S123 Steinhaus sequence addendum

HYP-2917 adds a useful exact sequence split for the three-gap/tight-locus route,
refining incoming HYP-2913/HYP-2914/HYP-2915 rather than replacing them.
For each denominator-14 unit `a`, retain the observer and all runner residues on
`Z/14Z`, with residue collisions recorded as zero cyclic gaps.  The AP row has
gap partition `((1,14))`; GW has the one-collision/one-missing partition
`((0,1),(1,12),(2,1))`.

The key correction is that this residue/gap partition is only a filter.  The row
`{1,...,11,13,36}` has the same coarse GW partition and is apex-locked, but
escapes off apex with

```text
M({1,...,11,13,36}) = 3/41 > 1/14.
```

Thus the missing Node-2 theorem should now be stated as:

```text
every non-AP/GW apex-locked reduced row has an off-apex escape M(S)>1/14.
```

S123 exact banks support this formulation: single AP replacements `v<=80` have
unique tight sets AP/GW and `805` apex-locked loose rows; local two-swaps
`<=18` have no non-AP tight row; bounded q-covering rows in `[1,19]` have
minimum `M=1/12`.  This does not close the global census, but it prevents a
false residue-only proof and gives a sharper finite sequence target.

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

THM-568/HYP-2929 improves this status: the equality-to-apex-shell part is
proved elementarily and formalized in Lean, so the remaining atom-forcing
theorem can be stated as

```text
14-covering bounded atom  =>  not tight, in fact M>1/14.
```

The `>=7` multiples-of-14 subcase is closed by THM-571.  The remaining
template risk is the `<=6` side only insofar as S31v's comb-teeth union bound
still has to be packaged with the bounded/intermediate finite-core census and
scale-separation reduction.

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

After THM-568/THM-571, the sharpest remaining version is:

```text
S = C union L, where L contains at most six scale-large multiples of 14.
Use S31v's comb-teeth union bound plus the bounded/intermediate finite-core
census to prove the danger combs cannot cover the C-safe interval.
```

Packaging this finite-core/scale-separated branch would make the
THM-079 template genuinely close LRC14.
