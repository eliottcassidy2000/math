---
id: HYP-2929
status: CORRECTION / PROOF-TARGET; Lean formalized safe arithmetic core
source: codex-2026-06-22-S120
tags: [lrc14, apex7, tightness, denominator, binding-pair, lean, thm568-correction]
depends_on:
  - HYP-+2909
  - HYP-2909
  - HYP-2910
  - THM-523
  - THM-568
related:
  - OPEN-Q-108
  - THM-079
  - THM-569
  - HYP-2908
  - KPS-S31ab
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCApexShell.lean
  - 04-computation/lean/TournamentH7/TournamentH7/LRCBindingPair.lean
results:
  - 05-knowledge/results/lrc_apex_shell_verify_codex_s120.out
  - 04-computation/lrc_14covering_not_tight_kps.py
  - 05-knowledge/results/lrc_14covering_not_tight_kps.out
---

# HYP-2929: the apex-shell correction for LRC14 tightness

THM-568's local binding-pair computation is valuable, but the currently written
strong conclusion

```text
primitive tight M(S)=1/14  =>  denominator D = 14
```

does not follow from the displayed arithmetic alone.

The theorem-level conclusion presently justified is the weaker apex-shell
statement.  If a lowest-term optimum `t=a/D` has opposite-side binders at
`+1/14` and `-1/14`, then

```text
14 | D,
D | (u+v),
14 | (u+v).
```

Equivalently, the tight point lies on an apex shell `D=14h`, and its active
pair has sum divisible by that shell.  The missing shell-collapse theorem is:

```text
primitive tight atom  =>  h=1.
```

This is exactly the "s_i+s_j = 14 exactly" refinement already listed as
remaining in HYP-+2909.  Primitivity of the full row does not by itself force
the active pair sum to be `14`; a primitive row can contain an active pair with
sum `28`, `42`, etc. unless an additional tight-locus or state-lift argument
rules that out.

## Formalized core

`LRCApexShell.lean` formalizes the safe part sorry-free.  In integer form, the
two active equations are

```text
14 * (u*a - m*D) =  D
14 * (v*a - n*D) = -D.
```

Lean proves:

```text
fourteen_dvd_den_of_pos_binding
fourteen_dvd_den_of_neg_binding
den_dvd_pair_sum_of_opposite_bindings
fourteen_dvd_pair_sum_of_opposite_bindings
```

The last theorem is the rigorous HYP-+2909 readout: tight opposite binders are
antipodal modulo `14`.

Incoming S50's `LRCBindingPair.lean` is the denominator-14 unit-grid companion:
at `t=a/14`, residues `si*a ≡ 1` and `sj*a ≡ -1 mod 14` imply
`14 | (si+sj)`.  `LRCApexShell.lean` extends the same arithmetic to arbitrary
apex shells `D=14h` by proving `D | (u+v)` from the two integer binding
equations.  These are compatible layers, not competing claims.

## Incoming Evidence

KPS-S31ab (`lrc_14covering_not_tight_kps.py`) is supporting evidence in the
right direction: in AP/GW-to-multiple-of-14 perturbations it finds minimum
`M=1/13>1/14`, no primitive tight 14-covering row, and no primitive tight row
with `D!=14` in its scan.  This should be read as a perturbation atlas, not a
general proof of shell collapse or of the full `R union 14A` residual.

## Corrected proof state

The THM-079 template should now be routed through one of these obligations:

```text
Shell collapse:
  primitive tight apex-shell row => D=14.

Covering strictness without shell collapse:
  if S contains a multiple of 14, then M(S)>1/14 even when the active shell is
  D=14h with h>1.

State lift:
  any h>1 apex-shell over-cover realizes the forbidden HYP-2908 K3 packet.
```

The third statement is the cleanest tournament-style target: shell height `h`
is a labelled packet, not a nuisance denominator.  The proof must show that an
attempt to keep `h>1` while covering the apex creates the same connected
`I(.,2)=7` atom that THM-201/THM-343 forbid.

## Assumption Challenge

Candidate carriers considered: runners, active binding pairs, denominators,
shell heights `h=D/14`, denom-14 residues, q-covering obligations, and HYP-2908
conflict packets.

The useful quotient for this correction is the active shell `(h, u+v)`.  It
preserves the local sawtooth arithmetic and the antipodal mod-14 predicate, but
it destroys the rest of the speed row.  That is precisely the challenged
assumption: local active-pair arithmetic alone cannot see row primitivity or
covering strictness.

## Impact

This does not weaken the LRC14 strategy; it prevents a false shortcut.  The
proof target is sharper:

```text
LRC14 = reductions + apex-shell arithmetic + shell-collapse/covering-strictness.
```

The `>=7` multiples-of-14 residual remains live, but should be read as an
apex-shell equidistribution problem rather than as a consequence of a proven
`D=14` denominator theorem.
