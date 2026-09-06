---
id: THM-4405
title: "LRC14 one-zero relation norm sixteen through twenty atlas"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4398 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. All 26 one-zero-mod-three primitive relation-presentation sectors
  of coefficient norm 16, 18, or 20 have sharp maximum at most 6/77. Five
  sectors attain equality, always at the single physical comb (1,5,11).
  Combined with THM-4398, this closes the 40 such presentation sectors through
  norm 20, not arbitrary speed triples or LRC(14).
source: root + jc4385_elliptic / LRC14 continuation session, 2026-09-03
depends_on:
  - THM-4398-lrc14-one-zero-relation-residue-dichotomy-and-small-norm-atlas
related:
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
  - THM-4393-lrc14-minimal-ternary-unit-norm-eighteen-shell
  - THM-4394-lrc14-minimal-ternary-unit-norm-twenty-shell
primary_script: 04-computation/lrc14_one_zero_relation_norm16_to20_atlas_thm4405.py
primary_output: 05-knowledge/results/lrc14_one_zero_relation_norm16_to20_atlas_thm4405.out
primary_script_sha256: 66d53f2c466b7e3f41fb6df7fd8438e0b3222b98de0bdcac8551f3dbe3807d09
primary_output_sha256: 9c7f810da6a2f16fdb1b7676863f48e6abc0adbf7bb9a554c3bd0ea5bf77c835
independent_audit_script: 04-computation/lrc14_one_zero_relation_norm16_to20_atlas_thm4405_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_one_zero_relation_norm16_to20_atlas_thm4405_independent_audit.out
independent_audit_script_sha256: 62a758084d5a7a0a0c66788315e73a0fb4d0b48bbed74ba9d74ef1782edfdf0d
independent_audit_output_sha256: da0739060197f1546b5b344a7bf8e81767edbfdb493331e28ffb497773427bfb
hash_basis: raw LF bytes
audit: >
  PASS. A smaller standard-library-only referee independently generates the
  coefficient shapes, signed relation rays, speed universes, cube slices,
  physical interval dictionaries, affine carrier dictionaries, admissible
  cutoffs, and maxima. It checks 646 labelled charts carrier by carrier and
  records 16,090 optimization-live gates. Normal, optimized, and fixed-seed
  outputs byte-match; the optimized tripwire rejects.
---

# THM-4405 -- LRC14 one-zero relation norm sixteen through twenty atlas

**PROVED ELEMENTARY RELATIVE TO THM-4398 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** These are overlapping relation-presentation sectors, not a
partition into minimal-relation cases. The theorem proves no arbitrary
nonresonance, chart-entry, synchronization, or fourteen-runner result.
`LRC(14)` remains **OPEN**.

## 1. Complete coefficient-shape list

For each `N in {16,18,20}`, enumerate

```text
1<=a<=b<=c,  a+b+c=N,  gcd(a,b,c)=1,
exactly one of a,b,c is divisible by 3.                 (1)
```

There are exactly `4+12+10=26` shapes:

```text
16: (2,3,11) (2,5,9) (3,5,8) (5,5,6)

18: (1,2,15) (1,3,14) (1,5,12) (1,6,11) (1,8,9)
    (2,3,13) (2,7,9) (3,4,11) (3,5,10) (3,7,8)
    (4,5,9) (5,6,7)

20: (1,1,18) (1,3,16) (1,4,15) (1,6,13) (1,7,12)
    (1,9,10) (3,4,13) (3,7,10) (4,7,9) (6,7,7).       (2)
```

Every distinct coordinate permutation and mixed sign choice is retained,
modulo only overall sign. A row consists of sorted, distinct, primitive,
positive odd speeds, all nonzero modulo three, satisfying one such labelled
relation. No minimal-relation filter is imposed.

## 2. Inherited residue and tail mechanism

For a signed primitive relation `c dot w=0`, THM-4398 proves that exactly one
coefficient divisible by three forces

```text
delta=c dot n !=0 (mod 3),
C in C_delta+Zc,
one live k class modulo 3 in each defect fibre.         (3)
```

Strictness `|delta|<(3/14)||c||_1` gives the exact defect sets

```text
norm 16 or 18: {-2,-1,1,2},
norm 20:       {-4,-2,-1,1,2,4}.                       (4)
```

Every affine carrier is evaluated by the exact roof

```text
L_w(C)=max(0,min(
  3/(7w1), 3/(7w2), 3/(7w3),
  3/(14w1)+3/(14w2)-|C3|/(w1w2),
  3/(14w1)+3/(14w3)-|C2|/(w1w3),
  3/(14w2)+3/(14w3)-|C1|/(w2w3))).                    (5)
```

For a magnitude shape `p=(a,b,c)`, both audited implementations reconstruct
the cube slice

```text
I_delta(p)=1/(2abc) sum_(S subset {1,2,3}) (-1)^|S|
 (|delta|+(3/14)(a+b+c)-(3/7)sum_(i in S)a_i)_+^2.     (6)
```

One route uses signed box convolution; the other clips the corresponding
rational polygon. They agree exactly. The step-three unimodal sampling lemma
from THM-4398 gives, for `W=max(w)`,

```text
mu(F_w)<=B_p+E_p/W,
B_p=(1/3)sum_(delta in D_p) I_delta(p),
E_p=3|D_p|/7.                                          (7)
```

After the finite row maximum `M_p` is computed, set

```text
T_p=E_p/(M_p-B_p).                                     (8)
```

The proof universe ends at the largest positive odd ternary-unit height
`H_p<=T_p`. Each `H_p` is occupied, `(7)` at `H_p` is not yet strict, and the
next admissible height exceeds `T_p` and is strictly below `M_p`. Thus the
finite cutoff is exact for this tail argument rather than an arbitrary box.

## 3. Sharp atlas

The following table records every new sector. `H` is the admissible cutoff;
each maximum and winner is unique inside its presentation row.

| shape | `H` | sharp mass | physical winner |
|---|---:|---:|---:|
| `(2,3,11)` | 31 | `58/833` | `(5,7,17)` |
| `(2,5,9)` | 29 | `12/161` | `(1,11,23)` |
| `(3,5,8)` | 25 | `6/77` | `(1,5,11)` |
| `(5,5,6)` | 31 | `8/119` | `(5,11,17)` |
| `(1,2,15)` | 37 | `60/1001` | `(1,11,13)` |
| `(1,3,14)` | 31 | `12/175` | `(1,13,25)` |
| `(1,5,12)` | 31 | `64/931` | `(7,13,19)` |
| `(1,6,11)` | 35 | `6/91` | `(1,7,13)` |
| `(1,8,9)` | 37 | `118/1925` | `(7,11,25)` |
| `(2,3,13)` | 29 | `12/161` | `(1,11,23)` |
| `(2,7,9)` | 37 | `10/161` | `(5,13,23)` |
| `(3,4,11)` | 31 | `46/665` | `(5,7,19)` |
| `(3,5,10)` | 31 | `58/833` | `(5,7,17)` |
| `(3,7,8)` | 25 | `6/77` | `(1,5,11)` |
| `(4,5,9)` | 37 | `10/161` | `(5,13,23)` |
| `(5,6,7)` | 31 | `64/931` | `(7,13,19)` |
| `(1,1,18)` | 65 | `2/37` | `(1,19,37)` |
| `(1,3,16)` | 59 | `12/203` | `(5,17,29)` |
| `(1,4,15)` | 47 | `58/833` | `(5,7,17)` |
| `(1,6,13)` | 41 | `6/77` | `(1,5,11)` |
| `(1,7,12)` | 49 | `8/119` | `(5,11,17)` |
| `(1,9,10)` | 59 | `12/203` | `(7,11,29)` |
| `(3,4,13)` | 41 | `6/77` | `(1,5,11)` |
| `(3,7,10)` | 53 | `6/91` | `(1,7,13)` |
| `(4,7,9)` | 41 | `6/77` | `(1,5,11)` |
| `(6,7,7)` | 49 | `64/931` | `(7,13,19)` |

Consequently every new row is at most `6/77`, and equality occurs in exactly

```text
(3,5,8), (3,7,8), (1,6,13), (3,4,13), (4,7,9),        (10)
```

always at the same physical comb `(1,5,11)`. Combined with THM-4398, all 40
one-zero presentation shapes through norm 20 have maximum at most `6/77`;
ten presentation shapes attain equality, but there is still only one physical
equality comb.

## 4. Physical dictionaries and overlap

The primary and referee construct the finite universes independently by
relation-solving and brute triple generation. For each physical triple, a
definition-level sweep intersects the exact open nearest-integer intervals,
retains pairwise-distinct owner residues, and labels each component by
`C=w cross n`. Separately, every relation chart enumerates `(3)` and evaluates
`(5)`. The complete carrier-to-length dictionaries agree for all `646`
labelled relation incidences—not just their total masses.

The 26 proof boxes contain

```text
row triple incidences:                  636
labelled relation incidences:           646
positive row incidences:                559
physical carrier components:          1638
distinct physical triples in the union: 370
positive physical triples in the union: 340.           (11)
```

Every one of the 26 row winners also satisfies a signed `(1,1,2)` relation.
This is an exact finite compression result for these winners, proved by the
complete cross-product intersection lists; it is not extrapolated to future
norms.

At `(1,5,11)`, the five new charts all reconstruct the same two physical
components of length `3/77`, but their defects split in two ways:

```text
profile A: delta=-1,+1, mass 3/77 each
  (3,5,8), (1,6,13), (4,7,9);

profile B: delta=-2,+2, mass 3/77 each
  (3,7,8), (3,4,13).                                  (12)
```

This is a concrete warning that relation charts overlap and defect is not a
physical invariant. The raw carrier dictionary is the consumer-complete
address for the local mass.

## 5. Reproduction and boundary

Both verifiers are standard-library-only and use exact `Fraction` arithmetic.
All correctness gates raise explicitly under optimization.

```powershell
python -B 04-computation/lrc14_one_zero_relation_norm16_to20_atlas_thm4405.py
python -B -O 04-computation/lrc14_one_zero_relation_norm16_to20_atlas_thm4405.py
python -B 04-computation/lrc14_one_zero_relation_norm16_to20_atlas_thm4405_independent_audit.py
python -B -O 04-computation/lrc14_one_zero_relation_norm16_to20_atlas_thm4405_independent_audit.py
```

The primary records `67,051` live checks and the referee `16,090`. Normal,
optimized, and fixed-seed outputs byte-match the frozen LF artifacts; both
optimized tripwires reject when invoked.

The theorem closes 26 presentation sectors, not all relations of larger norm,
not all speed triples, and not the missing seam-entry or synchronization
steps. `LRC(14)` remains **OPEN**.
