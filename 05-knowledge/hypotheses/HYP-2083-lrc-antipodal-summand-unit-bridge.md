---
id: HYP-2083
status: SUPPORTED by structural synthesis + S571 audit
source: codex-2026-06-03-S571
related:
  - HYP-2084
  - HYP-2082
  - HYP-2081
  - HYP-2075
  - HYP-2059
  - HYP-2055
  - HYP-2052
---

# HYP-2083: the `2n-1` antipodal witness reduction is the summand graph acted on by multiplicative units

The S553 witness family at modulus `C=2n-1` is the summand graph node `C`
equipped with the action of `(Z/CZ)^*`.

Dictionary:

```text
summand pair {a,C-a}      <-> antipodal residue shell P_a
missed pair               <-> candidate 2/C witness
a unit mod C              <-> invertible clock k=a^{-1}
nonunit a                 <-> composite-modulus blind spot
perfect transversal       <-> one speed in every summand shell
AP                        <-> all-lower transversal
```

**Claim.** A missed antipodal pair gives the S553 `2/(2n-1)` witness exactly
when the missed pair lies in a unit-visible shell.  If every missed pair is
nonunit, the instance lies in the composite-modulus hole and must be handled by
a second clock or gcd-stratum descent.

**Odd/even reading.**  Odd `C` has no midpoint, so the summand pairs are exactly
the cyclic antipodal pairs.  Even `C` has the fixed midpoint `C/2=-C/2`, and the
distinct-summand graph excludes the degenerate midpoint pair; this is the same
shape as the LRC apex/half-turn obstruction.

**Audit evidence (`lrc_antipodal_summand_units_s571.py`).**
- Prime `C=11,13`: all antipodal shells are unit-visible, one unit-action orbit.
- Composite `C=15`: shells split into unit, gcd-3, and gcd-5 strata; the n=8
  sporadics miss only nonunit shells.
- Composite `C=27`: shells split into unit, gcd-3, and gcd-9 strata; the n=14
  tight sporadic `V*={1,2,3,4,5,6,7,8,9,10,11,13,24}` misses `{12,15}` and
  doubles `{3,24}`, both in the gcd-3 stratum.
- Follow-on audit (`lrc_second_gap_transversal_audit_s572.py`): in bounded
  primitive boxes through `k=6`, every row below `2/(2n-1)` is already
  `n`-clock tight with `M(S)=1/n` and is a perfect antipodal transversal; the
  bounded flip-set menu is AP plus the known `{2}` sporadics.

**Proof-shape consequence.**  The reduced spectral-gap problem should split into:

```text
1. perfect transversals / flip-set structure;
2. nonunit gcd-stratum defects;
3. a second witness family or lift that restores invertibility for nonunit shells.
```

This reframes addition versus multiplication:

```text
addition       -> creates antipodal/summand shells
multiplication -> unit inverses select visible shells as clocks
```

**Honest scope.**  This is a structural bridge and bounded audit, not a proof of
the spectral gap.  The open theorem is to close the nonunit-hole branch, likely
by lifting the modulus or coupling the `n`-clock with the `2n-1` unit-action
clock.

**See:** `07-reflections/lrc-2n-minus-1-summand-unit-bridge-s571.md`,
`04-computation/lrc_antipodal_summand_units_s571.py` (+.out),
`04-computation/lrc_second_gap_transversal_audit_s572.py` (+.out);
S553, S553b, S559o, S560o, HYP-2084, HYP-2052, HYP-2055.
