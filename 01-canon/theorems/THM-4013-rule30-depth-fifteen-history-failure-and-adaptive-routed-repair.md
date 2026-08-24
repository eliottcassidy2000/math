---
id: THM-4013
title: "Rule 30 depth-fifteen projective-history failure and adaptive routed repair"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  canonical highest-shell universe n=2^m+t, 0<=m<=15, 0<=t<2^m, all 65,535
  records are reconstructed with the physical phase owner, exact carry,
  arithmetic shadow, projective histories, and routed chains. Three history
  levels leave 1,522 transition-mismatch fibres, including 65 cross-scale;
  the first cross-scale survivor is 1597/32967. Depths four, five, and six
  leave 197, 33, and 10 same-scale mismatches. The physical base chain leaves
  111; neither fixed off-ray chain is universal. Conditional on the base
  fibre, at most one selected full off-ray chain always resolves this finite
  universe, and each direction is independently required. This is a bounded
  lookup closure in one declared route family, not an all-scale observer,
  minimal finite-state theorem, gap bound, or Rule 30 prize consequence.
source: root + extended first-mismatch scout + independent no-import audit, 2026-08-24
audit: >
  PASS (independent standard-library reconstruction, 2026-08-24). The referee
  imports, executes, and reads neither producer nor THM-3511 helper. It derives
  valuations from the integer Rule 30 recurrence, rebuilds every finite-tree
  phase action, streams all 65,535 records, and exactly reproduces history,
  route, adaptive-selector, cross-scale, zero-base, carry, and direct-next
  controls. Normal and optimized streams are byte-identical. The exhaustive
  audit checks 131,070 tail images, 5,385,856 conjugacy entries, 196,605
  routes, 32,767 direct next-scale banks, and 98,304 same-scale top bits.
depends_on:
  - THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
  - THM-3824-rule30-fixed-division-tariff-and-physical-phase-separation
  - THM-4006-rule30-physical-cross-scale-observer-first-transition-mismatch
related:
  - THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary
  - THM-3804-rule30-all-period-amplitude-lattice-smith-law
script: 04-computation/rule30_depth15_history_adaptive_route_thm4013.py
output: 05-knowledge/results/rule30_depth15_history_adaptive_route_thm4013.out
script_sha256: 9d124486612ea62cd8c031dee4aad76975a1d630af41c998d1680cf8314306c7
output_sha256: a7fa61e33da8e9a3594ec3255e1b39597cfb7181c715390ce47fff68c0cbe16f
semantic_sha256: 6abf84f89c2d988a9c9eedf0b7d0e0a7f5a369e0faed6d467990d917b5369af6
independent_audit_script: 04-computation/rule30_depth15_history_adaptive_route_thm4013_independent_audit.py
independent_audit_output: 05-knowledge/results/rule30_depth15_history_adaptive_route_thm4013_independent_audit.out
independent_audit_script_sha256: 565e2d5e7cb1b1069dc65375d7dc2f4dad22fe8b08c596fb1257ea3ca999799c
independent_audit_output_sha256: 6630bc7e10c365c0b2a105526e8d46387d3e90ccf3211010c172c91335824a82
independent_audit_semantic_sha256: 836398e7aabefadaa1981d97bf2d0cba6c62c2c7b7a4c7c94e0badf22bbe0882
independent_record_stream_sha256: 5d203c998766dce0a06625767bc2e371138a4d2045e89403538e39736f31c17b
hash_basis: raw LF bytes
---

# THM-4013 -- depth-fifteen history failure and adaptive routed repair

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Retain the
Rule 30 orbit, signalizer gaps, carry decomposition, and physical phase
sections of THM-3511, THM-3516, THM-3824, and THM-4006. The exact universe is

```text
n=2^m+t,            0<=m<=15,            0<=t<2^m,     (1)
```

so it contains every integer `1<=n<=65535` exactly once. Every statement
below is exhaustive in `(1)` and makes no extrapolation to the next record.

## 1. The retained physical record

For each record, the observer reconstructs and retains:

1. the exact gap `d_m` and actual `d_m`-bit phase-tail extension `a`;
2. the physical phase owner on the complete depth-four tree above the gap;
3. four consecutive odd normalized units modulo `16`;
4. THM-3516's direct target parity and exact ordinary lower-residue carry;
5. up to six projective triples `(d,G mod 16,Z mod 16)`;
6. the physical routed chains

   ```text
   C_r=(a+r*2^d, s(a+r*2^d), s^2(a+r*2^d)), r=0,1,2; (2)
   ```

7. the center bit and next selected owner bank on rays `(0,1,2)`.

Here `s` is the actual phase owner obtained by finite-tree conjugation. The
route begins at the physical extension `a`; it is not replaced by the
marked-origin zero ray.

## 2. Finite projective history is not transition-compatible

Quotient by the retained arithmetic/owner record and the last `h` projective
levels. The exact transition censuses are

```text
h   fibres   nontrivial   mismatches   same-scale   cross-scale
0   10097       7854         7638         7578          1204
1   24302      15055        14231        14052          1467
2   54025       9551         8422         8132           420
3   63589       1873         1522         1458            65
4   65214        309          197          197             0
5   65442         86           33           33             0
6   65491         41           10           10             0. (3)
```

Thus additional history removes the tested cross-scale collisions after
depth four but does not make the observation operation-compatible. The first
depth-three mismatch remains THM-4006's `943/951`: both centers are one, with
next banks `(7,10,13)` and `(11,14,1)`.

The first genuine cross-scale depth-three survivor is

```text
(n,m,t)=(1597,10,573) and (32967,15,199),
gap=2, extension a=2,
owner bank=(9,12,11),
four-unit shadow=(9,11,9,3),
carry=(visible,0,0),
history=((1,7,5),(2,13,1),(2,13,9)),
centers=0,0,
next banks=(1,8,3),(9,4,15).                           (4)
```

The complete depth-four owner portrait also agrees. Equation `(4)` proves
that three histories plus exact ordinary carry still forget a load-bearing
transition coordinate.

## 3. Routed-chain repair and its sharp finite boundary

Start from the depth-three observer and append the physical base chain `C_0`.
The exact route census is

```text
added route                    fibres  nontrivial  mismatches
full C_0                        65351      183          111
C_0 + first address of C_1      65467       67           32
C_0 + first address of C_2      65469       65           30
C_0 + both first addresses      65502       32            8
C_0 + full C_1                  65490       44            9
C_0 + full C_2                  65485       49           14
C_0 + full C_1 + full C_2       65510       24            0. (5)
```

Both first off-ray addresses are insufficient, so the second routed
application is load-bearing. Neither fixed full direction works universally.

There is a sharper adaptive statement. First quotient by the `C_0` observer,
then select one full off-ray chain from that base fibre. The exact census is

```text
base already sufficient      65240 fibres / 65313 records
either C_1 or C_2 works          88 fibres /   176 records
C_1 specifically required       14 fibres /    28 records
C_2 specifically required        9 fibres /    18 records
both chains required             0 fibres /     0 records
unresolved                       0 fibres /     0 records.       (6)
```

The fibres containing `8556/10974` and `9922/12670` require `C_1` and `C_2`,
respectively. Hence no fixed direction can replace the selector. Within this
declared family, the sufficient adaptive addition is

```text
C_0 plus at most one base-fibre-selected full C_1 or C_2.        (7)
```

After the starting bank is known, `(7)` costs at most four literal owner
evaluations: two along `C_0`, then two along the selected off-ray chain. This
is a finite lookup closure, not a lower bound against arbitrary symbolic or
Mealy encodings and not an all-depth selector theorem.

## 4. Physical, carry, and same-scale controls

The two implementations verify

```text
physical owner-tail images       131070
finite-tree conjugacy entries   5385856
routed-chain images              196605
direct next-scale bank checks     32767
same-scale packed top bits        98304.                (8)
```

For every record, the direct target XOR and exact ordinary carry reproduce
bit `n` of `R_n`. Every available recovered upper bank agrees with a separately
constructed owner at scale `m+1`. The same-scale checks also verify the
positive margin `4c^2-13c+4` at each scale `0<=m<=15`; this explains why the
finite same-scale quotient collisions do not contradict THM-3824's separation
by the unbounded exact free defect.

The zero-base hostile is unchanged and remains decisive. At `n=6`, replacing
`a+r*2^d` by `r*2^d` predicts `(11,12,9)`, while the physical next bank is
`(15,8,5)`.

## 5. Independent reconstruction and scope

The producer uses the canonical THM-3511 gap helper. The standard-library
referee imports, executes, and reads neither it nor the producer. It generates
the integer Rule 30 orbit directly, derives all valuations, constructs every
marked finite-tree permutation and physical chart section, and hashes an
independently encoded stream of all `65,535` snapshots before running separate
fibre censuses. Its normal and optimized transcripts are byte-identical and
reproduce `(3)--(8)` without discrepancy.

This theorem stops at scale 15, owner depth four, four unit residues modulo
16, six projective residues modulo 16, and the declared three routed chains.
It does not cover `n>=65536`, deeper portraits, longer histories, greater
arithmetic precision, noncanonical phase charts, or arbitrary finite-state
encodings. In particular, it proves no finite signalizer graph, periodicity or
nonperiodicity, balance, density, complexity, center theorem, bounded gap, or
Rule 30 prize statement.

Reproduce from the repository root:

```text
python -B 04-computation/rule30_depth15_history_adaptive_route_thm4013.py
python -B -O 04-computation/rule30_depth15_history_adaptive_route_thm4013.py
python -B 04-computation/rule30_depth15_history_adaptive_route_thm4013_independent_audit.py
python -B -O 04-computation/rule30_depth15_history_adaptive_route_thm4013_independent_audit.py
```
