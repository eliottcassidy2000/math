---
id: THM-4121
title: "Sharp projective tail-multiplier residue compiler"
status: >
  PROVED ELEMENTARY PROJECTIVE RESIDUE COMPILER + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. For the fixed eleven-speed body U and paired tails
  {t,ct}, the common phase 9/19 has density 0 when c=0 modulo 19, density
  16/19 when c=+/-1 modulo 19, and density 14/19 for every other nonzero
  multiplier class. The 16/19 value is sharp among all multiplier classes at
  this body and phase. Every member, safe or not at that phase, stays outside
  THM-4112's three explicit thirteen-speed supplier shapes after every common
  dilation and primitive normalization. LRC(14) remains open.
source: codex-frontier-synthesis-creative-20260825g
depends_on: []
related:
  - THM-4049-lrc14-d2-two-phase-residue-firewall
  - THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families
  - THM-4117-physical-eleven-plus-two-primitive-support-obstruction
  - THM-4119-infinite-supplier-free-eleven-plus-two-residue-family
  - THM-4125-arbitrary-multitail-projective-packing-and-sharp-density
script: 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121.py
output: 05-knowledge/results/lrc_projective_tail_multiplier_residue_compiler_thm4121.out
independent_audit_script: 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_projective_tail_multiplier_residue_compiler_thm4121_independent_audit.out
script_sha256: c187f19b711e0ec108bdf9157d9ae60a5504e1cd9fa485780d7473911f69b011
output_sha256: 99acb8e686bbf3df8bd2f6b10e4ea5060e9c5db5a5a084a1fef046c1a7fada67
semantic_sha256: 62712d11643e5e93ba830a327226530739cd5f2bbace569c472751516166a089
independent_audit_script_sha256: 3c3228528ec854e887ac2d2085dabdb79ac7cb56a5d96b91e7df40e1e0a65022
independent_audit_output_sha256: 26e134f1873b92d399f44df283aa7fbdb235027440470e3c98816ee55e858e73
independent_semantic_sha256: 62712d11643e5e93ba830a327226530739cd5f2bbace569c472751516166a089
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone Fraction implementation classifies all 361 pairs of
  multiplier and parameter residues, checks the full 0/1/2 clearance rule on literal
  thirteen-speed representatives, verifies the projective inverse formula
  and sharp upper bound, and checks four parity cases through six dilation
  controls. Normal, optimized, and frozen outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room integer-residue implementation imports no primary
  code and independently exhausts every odd-tail choice and every possible
  AP7 scale divisor on all 361 residue representatives and six dilations. It
  reproduces the complete projective and clearance tables, direct-core failures, positive
  supplier controls, and sharp 16/19 ceiling. Normal, optimized, and frozen
  outputs byte-match; the semantic ledgers agree and the smallest theorem
  failure is none.
---

# THM-4121 -- sharp projective tail-multiplier residue compiler

**PROVED ELEMENTARY PROJECTIVE RESIDUE COMPILER + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**

THM-4119 fixed the tail multiplier at `3`. Allowing the multiplier itself to
move in the projective residue line reveals a sharper family: the two tail
failure pairs can coincide. This raises the valid parameter density from
`14/19` to `16/19`, and the first tail proves that no multiplier can do
better at the same body and phase.

## 1. The two-parameter family

Put

```text
U=(1,4,6,8,10,12,14,15,16,18,22),
S_(c,t)=U union {t,ct},                  c>1, t>22.         (1)
```

Both inequalities are integral. They make all thirteen speeds positive and
distinct. For an integer `a`, write

```text
|a|_19=min(a mod 19,19-(a mod 19)).                       (2)
```

The fixed body has exact clearance

```text
min_(u in U)||9u/19||=2/19.                               (3)
```

Let bars denote residues modulo `19`.

> **Theorem 1 (complete projective classification).** Fix an integer `c>1`.
>
> - If `c bar=0`, no residue class of `t` is safe at phase `9/19`.
> - If `c bar!=0`, define
>
>   ```text
>   B_c={0,+/-(9^(-1)),+/-(9c)^(-1)} in Z/19Z.            (4)
>   ```
>
>   Then `S_(c,t)` has exact clearance `2/19` at phase `9/19` if and
>   only if `t bar` is outside `B_c`.
> - The set `B_c` has size `3` when `c bar=+/-1` and size `5` for every
>   other nonzero multiplier class.
>
> Consequently the safe parameters have natural density `16/19` for
> `c=+/-1 mod 19`, density `14/19` for all other nonzero classes, and
> density zero for `c=0 mod 19`.

### Proof

At the displayed phase the two tail numerators are

```text
|9t|_19,                 |9ct|_19.                        (5)
```

If `c bar=0`, the second numerator vanishes for every `t`. Suppose now that
`c bar!=0`. A residue has distance below `2` from zero exactly when it is
`0,+1`, or `-1`. Thus the first numerator in `(5)` fails exactly on

```text
{0,+/-9^(-1)},                                             (6)
```

and the second fails exactly on

```text
{0,+/-(9c)^(-1)}.                                         (7)
```

Their union is `(4)`. Outside it both tail numerators are at least two, and
the body equality `(3)` makes the full clearance exactly `2/19`; inside it
one tail has clearance zero or `1/19`.

The two nonzero sign-pairs in `(4)` coincide precisely when

```text
(9c)^(-1)=+/-(9^(-1)),
```

equivalently `c=+/-1 mod 19`. Otherwise they are disjoint. This gives the
three cardinalities and, by periodicity in `t`, the asserted natural
densities. Since `2/19>1/14`, every safe row is `1/14`-lonely. **QED.**

## 2. Sharpness at the fixed body and phase

> **Theorem 2 (sharp multiplier ceiling).** Among every integer multiplier
> `c>1` in the paired-tail family `(1)`, the largest possible density at the
> fixed phase `9/19` is `16/19`. It is attained exactly by the multiplier
> classes `c=+/-1 mod 19`.

For every nonzero multiplier class, the first tail alone forces the three
bad parameter residues

```text
{0,+/-9^(-1)}={0,2,17}.                                  (8)
```

Hence at most sixteen of nineteen classes can survive. Theorem 1 shows that
`c=+/-1` aligns the second tail's bad pair with `(8)` and attains all sixteen.
A zero multiplier class is worse, not exceptional: it has no safe parameter.
This proves the phase-specific optimum. **QED.**

For example, the admissible integer multipliers `c=18` and `c=20` attain the
optimum. The THM-4119 choice `c=3` instead has

```text
B_3={0,2,7,12,17},                                        (9)
```

so Theorem 1 recovers its density `14/19` classification exactly.

THM-4125 subsequently replaces the pair `{1,c}` by an arbitrary finite tail
multiplier bank. Its fixed `1/14` density depends only on signed projective
support, while the natural threshold changes regime at seven tails.

## 3. Uniform exclusion from three explicit suppliers

For a finite positive set `A`, write `nu(A)=A/gcd(A)`. Since `1 in S_(c,t)`,

```text
nu(k S_(c,t))=S_(c,t)                 for every k>=1.      (10)
```

> **Theorem 3 (supplier exclusion).** For all integers `c>1`, `t>22`, and
> `k>=1`, neither `k S_(c,t)` nor its primitive normalization belongs to
> THM-4112's AP8+5, D0+6, or AP7+4 seam class.

The two direct cores fail already in the primitive row:

```text
{1,...,8} minus S_(c,t)={2,3,5,7},
{3,4,5,6,8,10,12} minus S_(c,t)={3,5}.                   (11)
```

If `k S_(c,t)` contained the AP8 core, its required speed `1` would force
`k=1`. If it contained the D0 core, its required speeds `3,4` would force
`k` to divide `gcd(3,4)=1`. Both cases reduce to `(11)`.

It remains to exclude an AP7 seam

```text
2qC union {x,y},                                          (12)
```

where `x,y` are odd and the primitive eleven-speed core `C` contains
`{1,...,7}`. Such a row has exactly two odd speeds. If `k` is even,
`k S_(c,t)` has none. Suppose `k` is odd. If `t` is odd and `c` is odd, the
four speeds `k,15k,kt,kct` are odd. If `t` is odd and `c` is even, exactly
the first three are odd. Neither case can have the form `(12)`.

Finally suppose `t` is even. The only odd speeds are `k,15k`, so they are
forced to be the two tails in `(12)`. The remaining even block is

```text
k E_(c,t),
E_(c,t)=(4,6,8,10,12,14,16,18,22,t,ct),
gcd(E_(c,t))=2.                                           (13)
```

Because `C` contains `1`, it is primitive. Equations `(12)--(13)` therefore
force `2q=gcd(kE_(c,t))=2k` and

```text
C=E_(c,t)/2=(2,3,4,5,6,7,8,9,11,t/2,ct/2).              (14)
```

The forced core misses `1`, contradicting its AP7 core. This excludes all
three suppliers for every dilation. Equation `(10)` handles primitive
normalization. **QED.**

## 4. Scope and exact replay

The `16/19` optimum is for this fixed body, the paired-tail form `{t,ct}`,
and the single common phase `9/19`. It is not an upper bound for arbitrary
phases, bodies, or LRC constructions. A parameter in a bad class only fails
this displayed phase; it is not asserted to be nonlonely. Likewise the
supplier conclusion concerns THM-4112's **three explicit thirteen-speed
shapes**, not every possible input to its general component-ancestry lemma.
Unlike the `c=3,t=2^45` specialization of THM-4117/4119, the new optimal
multiplier rows are not asserted to be physical residual rows.

The primary audit uses exact Fractions and literal representatives for all
`19*19` residue pairs. The clean-room audit uses only integer residue
distances, independently exhausts every possible odd-tail choice and AP7
scale divisor, and includes positive controls for all three supplier
recognizers. Run

```text
python3 -B 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121.py
python3 -B -O 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121.py
python3 -B 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121_independent_audit.py
python3 -B -O 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_projective_tail_multiplier_residue_compiler_thm4121_independent_audit.py
```

All six streams match their frozen outputs byte-for-byte. The complete bad
sets for every multiplier residue are recorded there. This compiler enlarges
an explicit supplier-free family; it does not prove LRC(14). **QED.**
