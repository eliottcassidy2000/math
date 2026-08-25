---
id: THM-4119
title: "Infinite eleven-plus-two residue family outside three explicit ancestry suppliers"
status: >
  PROVED ELEMENTARY CONGRUENCE FAMILY + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. For 14 of the 19 residue classes of t, the thirteen-speed set
  U union {t,3t} is exactly 2/19-safe at the common phase 9/19. For every
  t>22, independently of that phase condition, no common dilation or
  primitive normalization belongs to any of THM-4112's three explicit
  thirteen-speed supplier classes. The family has natural density 14/19 in
  its parameter and includes THM-4117's canonical t=2^45 row. LRC(14) remains
  open.
source: codex-frontier-synthesis-creative-20260825g
depends_on: []
related:
  - THM-4049-lrc14-d2-two-phase-residue-firewall
  - THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families
  - THM-4117-physical-eleven-plus-two-primitive-support-obstruction
script: 04-computation/lrc_supplier_free_residue_family_thm4119.py
output: 05-knowledge/results/lrc_supplier_free_residue_family_thm4119.out
independent_audit_script: 04-computation/lrc_supplier_free_residue_family_thm4119_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_supplier_free_residue_family_thm4119_independent_audit.out
script_sha256: e246ac0d3fb07f07b9be930c8d547b14202b3a166ea15a0d0bed6865480a3697
output_sha256: 538ad0b9cbbcc000fb49e8867f8825db07595fd3ff9c3efd73d351055ca42584
semantic_sha256: b1d6101e71ed922a4f1a873dcd11d0f7c41be9bd1f2a3c1a6c53f8f07c61ac35
independent_audit_script_sha256: cd78f8fa07ab9e2fc3a16ab619def789d0f8d3cd5abab6143f5ed0716333bcc1
independent_audit_output_sha256: e777357295fb19006158fff2f5788aa6deec2574e7fab8b16b42af8109f66860
independent_semantic_sha256: b1d6101e71ed922a4f1a873dcd11d0f7c41be9bd1f2a3c1a6c53f8f07c61ac35
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone Fraction implementation classifies all nineteen phase
  residues, checks representative thirteen-speed rows, both direct-core
  exclusions, the parity/gcd AP7 failure, six dilation controls, and the
  canonical t=2^45 specialization. Normal, optimized, and frozen outputs
  byte-match.
independent_audit: >
  ACCEPT. A clean-room integer-residue implementation imports no primary code
  and exhausts all odd-tail choices and every possible AP7 scale divisor. It
  independently reproduces the sharp five-residue failure set, all fourteen
  safe classes, direct-core failures, dilation controls, and positive supplier
  controls. Normal, optimized, and frozen outputs byte-match; the semantic
  ledgers agree and the smallest theorem failure is none.
---

# THM-4119 -- an infinite supplier-free 11+2 residue family

**PROVED ELEMENTARY CONGRUENCE FAMILY + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.**

THM-4117's physical hostile is one member of a much larger positive family.
The fixed body below supports one rational phase on fourteen congruence
classes of the pair scale. At the same time, its literal primitive support
prevents every member—not only the safe ones—from entering THM-4112's three
explicit thirteen-speed supplier shapes.

## 1. The family and sharp phase classification

Put

```text
U=(1,4,6,8,10,12,14,15,16,18,22),
S_t=U union {t,3t},                         t>22.           (1)
```

The inequality makes the thirteen speeds distinct. For an integer `a`, write

```text
|a|_19=min(a mod 19,19-(a mod 19)).                       (2)
```

At the phase `theta=9/19`, exact residue arithmetic gives

```text
min_(u in U)||u theta||=2/19,
||t theta||=|9t|_19/19,
||3t theta||=|8t|_19/19.                                 (3)
```

> **Theorem 1 (sharp residue family).** Let
>
> ```text
> R={1,3,4,5,6,8,9,10,11,13,14,15,16,18} subset Z/19Z.  (4)
> ```
>
> If `t>22` and `t mod 19 in R`, then
>
> ```text
> min_(v in S_t)||9v/19||=2/19>1/14.                      (5)
> ```
>
> Hence every such `S_t` is `1/14`-lonely.

The inverse residues are `9^(-1)=17` and `8^(-1)=12` modulo `19`. Therefore
one of the two high speeds has residue distance at most one exactly when

```text
t mod 19 in {0,+/-17,+/-12}={0,2,7,12,17}.               (6)
```

At residue `0` the clearance in `(3)` is `0`; at the other four residues in
`(6)` it is `1/19`. Outside `(6)`, both high-speed numerators are at least
two, while the body already attains two. This proves `(4)--(5)` and shows the
classification is sharp for the displayed phase. Since `R` has fourteen
classes, the valid parameters have natural density `14/19`.

## 2. Uniform exclusion from the three explicit suppliers

For a finite positive set `A`, write `nu(A)=A/gcd(A)`. Since `1 in S_t`,

```text
nu(kS_t)=S_t                         for every k>=1.        (7)
```

> **Theorem 2 (supplier exclusion).** For every `t>22` and every positive
> integer `k`, neither `kS_t` nor `nu(kS_t)` belongs to THM-4112's AP8+5,
> D0+6, or AP7+4 seam class.

The direct cores fail uniformly:

```text
{1,...,8} minus S_t={2,3,5,7},
{3,4,5,6,8,10,12} minus S_t={3,5}.                       (8)
```

If `kS_t` contained the AP8 core, its required speed `1` would force `k=1`.
If it contained the D0 core, the required speeds `3,4` would force
`k|gcd(3,4)=1`. Both cases return the failures in `(8)`.

An AP7 seam has the form

```text
2qC union {x,y},                                         (9)
```

where `x,y` are odd and `C` contains `{1,...,7}`. If `k` is even, `kS_t`
has no odd speeds, so `(9)` is impossible. Suppose `k` is odd. If `t` is
odd, `kS_t` has the four odd speeds `k,15k,kt,3kt`, again impossible. If `t`
is even, the only odd speeds are `k,15k` and hence are forced to be the two
tails. The remaining even block is

```text
kE_t,  E_t=(4,6,8,10,12,14,16,18,22,t,3t),  gcd(E_t)=2. (10)
```

Every allowed `C` is primitive because it contains `1`; thus `(9)--(10)`
would force `2q=gcd(kE_t)=2k`, followed by

```text
C=E_t/2=(2,3,4,5,6,7,8,9,11,t/2,3t/2).                 (11)
```

The forced core misses `1`. This proves the AP7 exclusion, and `(7)` proves
the normalized exclusion. **QED.**

## 3. Relation to the physical hostile

The canonical THM-4117 parameter satisfies

```text
2^45=18 modulo 19.                                       (12)
```

It therefore lies in `R`, and `(5)` recovers its exact `2/19` clearance.
Only this special parameter is imported by THM-4117 from THM-4049's physical
finite box. The arbitrary parameters in `(4)` are an elementary LRC family;
no claim is made that all of them are physical residual rows of THM-3818.

The conclusion is also deliberately relative to THM-4112's **three explicit
thirteen-speed shapes**. It does not say that its general component-ancestry
lemma has no other sufficient inputs. Rather, the theorem proves that AP7,
AP8, and D0 cores are not necessary for `1/14` loneliness.

## 4. Exact replay and scope

The primary audit uses exact Fractions. The clean-room audit builds the phase
table with integer residues and brute-forces every possible odd-tail choice
and AP7 scale divisor. Run

```text
python3 -B 04-computation/lrc_supplier_free_residue_family_thm4119.py
python3 -B -O 04-computation/lrc_supplier_free_residue_family_thm4119.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_supplier_free_residue_family_thm4119.py
python3 -B 04-computation/lrc_supplier_free_residue_family_thm4119_independent_audit.py
python3 -B -O 04-computation/lrc_supplier_free_residue_family_thm4119_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_supplier_free_residue_family_thm4119_independent_audit.py
```

All six streams match their frozen outputs byte-for-byte. The five excluded
residues are failures only of the fixed phase `9/19`; they are not asserted
to be nonlonely. This family does not prove LRC(14), classify arbitrary
eleven-speed bodies, or exhaust component-ancestry suppliers. **QED.**
