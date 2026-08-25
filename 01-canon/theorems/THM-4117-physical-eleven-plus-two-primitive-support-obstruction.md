---
id: THM-4117
title: "Physical eleven-plus-two primitive-support obstruction"
status: >
  PROVED RELATIVE TO THM-4049 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  No common dilation or primitive normalization of the canonical physical
  rank-eleven 11+2 row from THM-4049/MISTAKE-490 belongs to any of THM-4112's
  AP7+4 seam, AP8+5, or D0+6 supplier classes. The obstruction is
  the absolute primitive support and its unit origin. The row is itself
  1/14-lonely, so this is a method hostile rather than an LRC counterexample.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-4049-lrc14-d2-two-phase-residue-firewall
  - THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families
related:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-4105-primitive-reciprocal-phase-descent-and-quantitative-arrival
  - THM-4110-sparse-reciprocal-phase-graph-saturation-and-ap13-torsion-tariff
  - MISTAKE-490
script: 04-computation/lrc_physical_supplier_obstruction_thm4117.py
output: 05-knowledge/results/lrc_physical_supplier_obstruction_thm4117.out
independent_audit_script: 04-computation/lrc_physical_supplier_obstruction_thm4117_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_physical_supplier_obstruction_thm4117_independent_audit.out
script_sha256: c867edb1da4c0bea71f351722c32790774a3319c513d80102a02beaaea00ce3e
output_sha256: 56051358a6aa60142934c37fe8216ce81e4b42f6cfc190221c9cb347b3472649
independent_audit_script_sha256: da067a3f6145f889ca46ae06e2b9a1a71c71d7c5db3c659dacf25b01522cca71
independent_audit_output_sha256: 09979db4cd0c806899580f3b3900c193dbd982f2110328a314140350a0259637
independent_semantic_sha256: 73563bafba2bc5e8dd870aee18e32047a053eed8a9bce84ad33e7aaecceb0f88
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone Fraction implementation checks primitive normalization,
  all three supplier exclusions, the forced AP7 parity/gcd split, physical
  box membership, residues, phase clearances, six dilation controls, and a
  positive control for each supplier class. Normal, optimized, and frozen
  outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room implementation imports no primary code, enumerates
  every possible pair of AP7 odd tails and every admissible scale divisor,
  and evaluates phases by integer residues. It independently reproduces the
  unique failed AP7 attempt, both direct-core failures, all exact clearances,
  and all positive controls. Normal, optimized, and frozen outputs byte-match;
  the smallest theorem failure is none.
---

# THM-4117 -- physical 11+2 primitive-support obstruction

**PROVED RELATIVE TO THM-4049 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4112 supplies three explicit scale-separated families. This theorem
shows that the presently retained physical `11+2` type and residue data do
not force a row into any of them. The missing coordinate is not another
phase-sheet label: it is the literal primitive support, including its unit
origin.

## 1. The canonical physical row

For a finite positive integer set `A`, write

```text
g(A)=gcd(A),                 nu(A)=A/g(A).                 (1)
```

Let

```text
U=(1,4,6,8,10,12,14,15,16,18,22),
(p,q)=(1,3),                s=1,          t=2^45,
S=sU union t{p,q}
 =(1,4,6,8,10,12,14,15,16,18,22,2^45,3*2^45).            (2)
```

THM-4049 section 5, with the correction recorded in MISTAKE-490, proves
that `(2)` lies in THM-3818's finite physical box, has rank-eleven `11+2`
typing and no bounded crossing row. Its scale-one pair type is `(1,3)` since
`max U=22<2^45`. Those source facts are imported rather than reproved here.

Common dilation is compatible with primitive normalization and phase:

```text
nu(kA)=nu(A),
||(kv)(theta/k)||=||v theta||                 (k>=1).       (3)
```

Translation of speeds or separate dilation of the two components is not part
of `(3)` and is not asserted to preserve one common phase.

## 2. Supplier-exclusion theorem

The three THM-4112 targets are:

1. `AP8+5`: `{1,...,8}` plus five ratio-two outliers starting at `94`;
2. `D0+6`: `D0={3,4,5,6,8,10,12}` plus six ratio-two outliers starting at
   `16`;
3. an `AP7+4` seam `2rC union {x,y}`, where `r>=1`, `x,y` are distinct
   positive odd speeds, and `C` contains `{1,...,7}` plus four outliers
   satisfying either THM-4112's adaptive or parity-free gates.

> **Theorem.** For every positive integer `k`, neither `kS` nor
> `nu(kS)=S` belongs to any of these three classes.

### Direct cores

For the primitive row,

```text
{1,...,8} minus S = {2,3,5,7},
D0 minus S        = {3,5}.                                (4)
```

Thus `S` belongs to neither direct class. If `kS` were AP8+5, its required
speed `1` would force `k=1`, returning the first failure in `(4)`. If `kS`
were D0+6, its required speeds `3` and `4` would force `k|gcd(3,4)=1`, again
returning `(4)`.

### AP7 seam

The odd and even parts of `S` are

```text
O={1,15},
E=(4,6,8,10,12,14,16,18,22,2^45,3*2^45),
g(E)=2.                                                   (5)
```

If `k` is even, `kS` has no odd speeds and cannot have the two required odd
seam tails. If `k` is odd, parity forces those tails to be `kO`, leaving the
block `kE`. Since every allowed `C` contains `1`, it is primitive. Hence

```text
kE=2rC  =>  2r=g(kE)=2k  =>  r=k,
C=E/2=(2,3,4,5,6,7,8,9,11,2^44,3*2^44).                 (6)
```

The forced core in `(6)` misses `1`, so it cannot contain `{1,...,7}`. This
already fails before either four-outlier gate is consulted. Finally, since
`g(S)=1`, equation `(3)` gives `nu(kS)=S` for every `k`. **QED.**

## 3. Phase and residue boundary

The row is positively safe at

```text
theta=9/19,             min_(v in S)||v theta||=2/19>1/14. (7)
```

The phase `9/(19k)` gives the same clearance on `kS` by `(3)`. Thus the
normalization did not manufacture the obstruction by losing a physical safe
phase. Yet the same phase on `S` is not an antipodal THM-4112 certificate:

```text
min_(v in S) min(||v theta||,||v(theta+1/2)||)=1/38<1/14. (8)
```

The forced divided block in `(6)` has residues

```text
(2,3,4,5,6,7,8,9,11,32,40) modulo 56.                   (9)
```

At THM-4049's four firewall times the full-row clearances are exactly
`(1/28,1/28,1/56,1/56)`, whereas `(7)` is safe. This is a method hostile,
not an infeasibility witness or an LRC counterexample.

## 4. Exact bridge boundary

The connection exposed by `(2)--(6)` is

```text
source:     THM-4049's physical rank-eleven 11+2 row
target:     the union of THM-4112's three explicit supplier classes
map:        common dilation, primitive normalization, then the forced
            parity split for an AP7 seam
preserved:  support ratios, primitive support/parity, and physical clearance
destroyed:  raw common scale and the choice among phase lifts
missing:    absolute primitive support with its unit origin, plus a core-safe
            interval if the target is to be used
test:       direct missing cores (4), then the forced quotient (5)--(6)
outcome:    bridge false; the first AP7 failure is the missing unit 1.
```

In particular, the displayed rank/type and THM-4049 residue/parity
projection do not alone force supplier entry. This does not show that every
row of type `(1,3)` misses the suppliers, that failure of a sufficient
THM-4112 certificate is a cover, or that a future argument using actual
nonloneliness or a stronger survivor inequality cannot close a residual row.
The witness `(2)` is lonely.

## 5. Exact replay

The primary and clean-room audits respectively use a direct gcd recognizer
and exhaustive tail/scale enumeration with integer-residue phase arithmetic.
Run

```text
python3 -B 04-computation/lrc_physical_supplier_obstruction_thm4117.py
python3 -B -O 04-computation/lrc_physical_supplier_obstruction_thm4117.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_physical_supplier_obstruction_thm4117.py
python3 -B 04-computation/lrc_physical_supplier_obstruction_thm4117_independent_audit.py
python3 -B -O 04-computation/lrc_physical_supplier_obstruction_thm4117_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_physical_supplier_obstruction_thm4117_independent_audit.py
```

All six streams match their frozen outputs byte-for-byte. No claim of
LRC(14), a new lonely family, or completeness of the three THM-4112 supplier
classes is made. **QED.**
