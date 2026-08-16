---
id: THM-3477
title: "Moving-numerator Jordan--Pruefer response normal form for Jacobian effectivity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.  For the intermediate affine
  modification E_H=C[s,v,Y]/(sY-vH(v)) inside B=C[s,h], v=sh and
  Y=hH(v), the additive effectivity quotient B/E_H has an exact graded
  response normal form under multiplication by s.  A root of H at the
  transition multiplier v=0 produces an explicit reduced ladder of finite
  Jordan blocks; every nonzero root produces graded Pruefer arms.  The
  inverse limits are instead zero and formal power-series completions.
  THM-3383's two terminal orientations are exact specializations.  This is
  an effectivity-module refinement, not a new JC(2) exclusion.
audit: >
  The grading, CRT-stable transition, finite-interval and infinite-ray
  decompositions, kernel/cokernel counts, moving-numerator embedding into
  THM-3404, and direct-limit versus inverse-limit semantics were independently
  proof-audited.  A standard-library exact companion checks transition
  matrices, literal barcodes, both terminal orientations, inverse systems,
  moving numerators, and merged-root hostiles on a frozen finite universe.
  Normal and optimized replays are byte-identical.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3383-terminal-monomial-cone-polynomiality-fork
  - THM-3397-torsor-killing-versus-effective-boundary-valuations
  - THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction
related:
  - THM-3406-affine-modification-power-jets-and-principal-part-transgression
  - THM-3437-derived-boundary-jet-euler-conservation-and-prufer-recovery
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
script: 04-computation/jc_moving_numerator_jordan_prufer_response_thm3477.py
output: 05-knowledge/results/jc_moving_numerator_jordan_prufer_response_thm3477.out
script_sha256: 644cdb7508485bc808a5ddd237c3c68703f81661557e7b24574b73b405567db4
output_sha256: 4e681c77629ecbf24f83b07ed187e651ebc838e09818268fa53af541c03ead12
semantic_sha256: 9fca2dffaa87feb03b47505f2adb1ef197ef46978cc25530dc83c2a52dcd47d0
hash_basis: raw bytes
---

# THM-3477 -- moving numerators split finite Jordan and Pruefer responses

**PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.**

## 1. Statement

Let `A=C[v]`, let `0!=H(v) in A`, and define

```text
E_H=C[s,v,Y]/(sY-vH(v)),                                  (1)
B=C[s,h],                 v=sh,             Y=hH(v).      (2)
```

This map is injective.  Indeed, `E_H` embeds in `A[s,s^(-1)]` by
`Y=s^(-1)vH`; its image lies in `B`, and `s,v=sh` are algebraically
independent in `B`.  Moreover,

```text
B=E_H[v/s].                                                (3)
```

Let

```text
M_H=B/E_H                                                   (4)
```

be the **additive** quotient.  It is not a quotient ring.  Multiplication by
`s` preserves `E_H` and therefore acts on `M_H`.

Factor

```text
H(v)=c v^m product_(i=1)^b (v-alpha_i)^(e_i),              (5)

c!=0,       m>=0,       alpha_i!=0 distinct,       e_i>=1.
```

Put

```text
Pr_s=C[s,s^(-1)]/C[s].                                    (6)
```

Then the complete `C[s]`-module response normal form is

```text
M_H = M_0 direct-sum M_nz,                                (7)

M_0 ~= direct-sum_(j>=1) C[s]/(s^(ell_j))       if m>0,
M_0 =0                                           if m=0,  (8)

ell_j=j-floor(j/(m+1))=ceil(mj/(m+1)),                    (9)

M_nz ~= direct-sum_(i=1)^b direct-sum_(k>=0) Pr_s.       (10)
```

The grading carries more than the ungraded isomorphism:

```text
the j-th finite bar occupies
 [floor(j/(m+1))+1, j];                                  (11)

the (i,k)-th infinite ray is born in grade
 floor(k/e_i)+1.                                         (12)
```

When `m>0`, every positive finite length occurs once, and every multiple of
`m` occurs a second time.  At every nonzero root, exactly `e_i` new infinite
rays are born at each grade.  The maximal `s`-divisible submodule of `M_H`
is exactly `M_nz`; the zero-root arm `M_0` is reduced, although its finite
block lengths are unbounded.

This zero/nonzero-root dichotomy is the result.  It distinguishes a boundary
where the transition multiplier vanishes from one where it is a unit; root
multiplicity alone does not.

## 2. Exact graded quotient

Give `B=C[s,h]` the difference grading

```text
deg(s)=1,          deg(h)=-1,          deg(v)=0.           (13)
```

Every monomial has a unique normal form, so additively

```text
B = (direct-sum_(q>=0) s^q A)
    direct-sum (direct-sum_(r>=1) h^r A).                 (14)
```

Since `Y=hH(v)` has degree `-1`, the same calculation in `E_H` gives

```text
E_H = (direct-sum_(q>=0) s^q A)
      direct-sum (direct-sum_(r>=1) h^r H(v)^r A).        (15)
```

Distinct grades cannot cancel.  Therefore

```text
M_H ~= direct-sum_(r>=1) h^r Q_r,
Q_r=A/(H^r).                                               (16)
```

With `Q_0=0`, multiplication by `s` is the quiver transition

```text
T_r:Q_r -> Q_(r-1),             [f] |-> [v f].           (17)
```

Indeed `s h^r f=h^(r-1)vf` for `r>=2`, while `shf=vf` lies in the
degree-zero part of `E_H` when `r=1`.  Formula (17), not merely the list of
finite quotients, is the load-bearing operation sidecar.

## 3. CRT arms and the transition phase change

The factors in (5) are pairwise coprime.  CRT refines (16) into transition-
stable arms

```text
Q_r ~= A/(v^(mr))
       direct-sum (direct-sum_i A/((v-alpha_i)^(e_i r))).  (18)
```

Omit the first summand when `m=0`.

At a nonzero root `alpha_i`, multiplication by `v` is a unit.  If
`z_i=v-alpha_i`, then multiplying grade `r` by `v^r` conjugates (17) to the
ordinary truncation

```text
C[z_i]/(z_i^(e_i r)) -> C[z_i]/(z_i^(e_i(r-1))).          (19)
```

At the zero root, however, (17) is genuinely multiplication by the local
uniformizer:

```text
C[v]/(v^(mr)) --times v--> C[v]/(v^(m(r-1))).             (20)
```

This is the algebraic phase change behind (7)--(12).

For `r>=2`, the armwise transition dimensions are

```text
zero root:     dim ker=m+1,       dim coker=1;
alpha_i!=0:    dim ker=e_i,       dim coker=0.            (21)
```

At `r=1`, the target is zero and the kernel dimensions are `m` and `e_i`,
respectively.  In particular, (17) is surjective in every grade exactly when
`m=0`.

## 4. Finite intervals at v=0

Assume `m>0`.  In the zero arm use the basis

```text
v^k in Q_r,                    0<=k<mr.                   (22)
```

Under (20), the quantity

```text
j=k+r                                                        (23)
```

is invariant.  For fixed `j>=1`, the allowed grades are exactly

```text
floor(j/(m+1))+1 <= r <= j.                               (24)
```

Multiplication by `s` moves down this interval by one and kills its bottom.
It is therefore one Jordan block of length

```text
j-floor(j/(m+1))=ell_j.                                  (25)
```

The pairs `(r,k)` are partitioned by (23), proving (8), (9), and (11).

To read the multiplicities, write `j=(m+1)q+b`, `0<=b<=m`.  Then

```text
ell_j=mq                         if b=0,
ell_j=mq+b                       if 1<=b<=m.              (26)
```

Each positive integer therefore occurs once, with a second occurrence at
the positive multiples of `m`.

## 5. Infinite rays away from v=0

Fix a nonzero root `alpha_i` and use (19).  A monomial `z_i^k` occurs in
grade `r` exactly when

```text
0<=k<e_i r,
```

or

```text
r>=floor(k/e_i)+1.                                       (27)
```

For fixed `k`, the transition carries this basis element down through every
grade in (27) and kills it at the first missing grade.  The ray has no top.
As a `C[s]`-module it is exactly `Pr_s`: if its birth grade is `b`, send its
element in grade `b+n` to `s^(-(n+1))` modulo `C[s]`.  Distinct `k` give a
direct sum, proving (10) and (12).

Multiplication by `s` is surjective on every such ray, so `M_nz` is
`s`-divisible (indeed divisible as a `C[s]`-module).
The direct sum of finite cyclic blocks in (8) has no nonzero divisible
submodule.  This proves the maximality assertion.

## 6. Moving numerators inside the full principal-part tower

Set

```text
F(v)=vH(v).
```

THM-3404 gives the full localization module

```text
(E_H)_s/E_H
 =direct-sum_(r>=1) s^(-r) A/(F^r),                      (28)
```

whose transition is coefficient reduction.  The polynomial intermediate
ring `B` lands in it by the gradewise moving-numerator maps

```text
A/(H^r) -> A/(v^r H^r),             [f] |-> [v^r f].     (29)
```

These maps are injective.  At a nonzero root of `H`, `v^r` is a unit, so
(29) is the whole local principal-part arm and its Pruefer behavior survives.
At `v=0`, (29) is a proper moving sublattice.  When `m>0`, it removes exactly
`r` orders from the ambient `(m+1)r`-jet and changes the reduction tower into
the finite Jordan barcode (8); when `m=0`, it clears that ambient arm
completely.

Thus the new mechanism is a **moving-numerator Pruefer-to-Jordan
transmutation**, not a recopy of THM-3404's full principal parts.

## 7. Inverse limits are completions, not Pruefer modules

Keep the actual arrows (17); all inverse limits below use those bonding maps,
conjugated as in (19).  On the zero arm, the image in a fixed stage of
a sufficiently high stage is zero.  Hence the inverse system is
Mittag--Leffler and

```text
inverse-limit_r A/(v^(mr)) =0,              lim^1=0.       (30)
```

After the conjugation leading to (19), the nonzero-root arrows are ordinary
surjective jet truncations.  Therefore

```text
inverse-limit_r A/((v-alpha_i)^(e_i r))
 =C[[v-alpha_i]],                              lim^1=0.    (31)
```

The Pruefer modules in (10) are the actual direct unions of the finite
`s`-torsion chains inside the additive direct sum (16).  They are not the
inverse limits in (31).  Confusing those two transition systems would repeat
MISTAKE-397; a tower is its objects **and its arrows**.

## 8. THM-3383 specialization

Use THM-3383's notation `v=ut` and `L=1+cv`.

In its positive orientation `g=ae+1`, take

```text
s=u,        h=t,        H=v^(g-1)L^e,        Y=G_+,
m=g-1.                                                     (32)
```

In its negative orientation `g=ae-1`, take

```text
s=t,        h=u,        H=v^g L^e,            Y=G_-,
m=g.                                                       (33)
```

Equations (1)--(17) reproduce the exact intersections and response strings
of THM-3383.  They sharpen its statement “locally nilpotent with unbounded
response lengths” into the complete finite/transient versus
infinite/persistent barcode.  The legacy class `[h^q]` has exact length `q`;
its full-length component is already visible in the constant jet of the
nonzero `L` arm.

## 9. Hostiles, equality boundaries, and scope

- **Unit boundary.**  If `H` is a unit, then `Y` is a unit multiple of `h`,
  so `E_H=B` and `M_H=0`, as (5)--(10) prescribe.
- **No zero root.**  If `m=0`, the finite barcode is absent and the entire
  quotient is `s`-divisible.
- **The `g=1` positive terminal cell.**  It is exactly the preceding `m=0`
  boundary, not an omitted finite block.
- **Simple nonzero root.**  `e_i=1` gives one new infinite ray per grade; no
  multiplicity hypothesis is needed for Pruefer behavior.
- **Merged-root hostile.**  For

  ```text
  H_alpha=v^m(v-alpha)^e,
  ```

  every `alpha!=0` has the finite `m`-barcode plus `e` Pruefer births per
  grade.  At `alpha=0`, the two roots merge and the whole response becomes the
  pure finite barcode with parameter `m+e`.  Root collision changes operation
  type; it is not a continuous multiplicity-only specialization.
- **Additive type.**  `M_H` is an additive quotient and a `C[s]`-module, not
  a quotient algebra.  Multiplying two cosets is not defined.
- **No global Jacobian conclusion.**  The theorem classifies the effectivity
  debt of the displayed affine modification and THM-3383 cells.  It neither
  proves that every terminal Keller chart has this form nor excludes any
  remaining `C3`, `A4/S4`, `G1`, `JC(2)`, or `DC(2)` branch.

## 10. Exact companion

The standard-library companion checks `38,808` exact transition-matrix cells,
`1,764` literal barcode states, `168` nonzero-root unit conjugacies, `11,928`
moving-numerator incidences, `1,504` inverse-tower controls, `89` terminal
specialization cells with `712` legacy response strings, and `336` merged-root
cells.  Its universe is `m=0..7`, `e=1..6`, grades `1..8`, and block indices
`1..160`; the unit, `m=0`, `e=1`, grade-one, and merged-root boundaries are
explicit controls.  The all-parameter conclusions come from the proof above,
not from this finite bank.

Reproduce the byte-identical normal and optimized transcripts with

```bash
python3 04-computation/jc_moving_numerator_jordan_prufer_response_thm3477.py
python3 -O 04-computation/jc_moving_numerator_jordan_prufer_response_thm3477.py
```

## 11. Information contract

The information contract is exact:

```text
source:    E_H inside the polynomial intermediate ring B
target:    additive polynomiality debt B/E_H
operation: multiplication by the surviving polynomial coordinate s
preserved: grade, root label, jet multiplicity, transition arrows
lost:      ring multiplication, global chart entry, inverse effectivity
sidecar:   moving numerator v^r and the zero/nonzero location of each root
test:      H_alpha=v^m(v-alpha)^e at alpha!=0 and alpha=0.             (34)
```

**QED.**
