---
id: THM-2435
title: "Top-blocker essential parent and equivariant quotient-root section"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT
  HOSTILE AUDIT. In each of the three surviving deep-c_3 shapes, the
  common-91 quotient stalk has parent mass at least (4+k)/7, where k
  is the number of low quotient blockers. THM-2431 caps the six-unit
  exact-tiling locus at 3/7, so actual exact-depth top blockers fill
  gaps over parent mass at least (1+k)/7. One fixed blocker label
  supports a translation-equivariant marked quotient-root section of
  uniform mass at least 3/1274. Every one of the 91 quotient residue
  classes has exact energy rho/91, hence at least 3/115934; the 72
  unit classes have total energy at least 108/57967. For M=0 this is
  a physical all-91-class packet. For M>0, physical pullback retains
  a 7^M-fold ancestry kernel and no 91-unit current. An intrinsic
  nonempty top-blocker word has quotient mass at least 1/637. Every
  essential F_13 row has a nonzero integral zero-sum F_7 incidence
  defect, cyclic Dirichlet energy at least 4, and all six nontrivial
  row colours; its unnormalized parent-summed energy integrates to at
  least 8/7. The full defect either is independent of its F_13 row or
  has a genuine 91-unit mixed quotient character. The marked root is
  a derived 13-copy Boolean selector, not a canonical Abel endpoint
  probe. No physical positive-depth 91-unit current, terminal current,
  valuation-shape exclusion, row decrement, or proof of LRC(14)
  follows.
source: codex-2026-07-26-top-blocker-essential-parent
depends_on:
  - THM-2427-guard-top-thirteen-root-capacity-and-residual-types
  - THM-2431-repeated-step-rounding-exclusion-of-guard-top-zero-blocker-types
  - THM-2432-guard-top-pair-cage-and-low-blocker-residual-exclusion
  - THM-2424-coprime-common-root-crt-and-unit-residue-spectrum
related:
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2409-unfiltered-septimal-source-orbit-and-real-word-obstruction
  - THM-2419-valuation-normalized-homogenization-of-affine-sideband-shells
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - MISTAKE-273
script:
  - 04-computation/lrc14_top_blocker_essential_parent_thm2435.py
  - 04-computation/lrc14_top_blocker_essential_carrier_thm2435.py
output:
  - 05-knowledge/results/lrc14_top_blocker_essential_parent_thm2435.out
  - 05-knowledge/results/lrc14_top_blocker_essential_carrier_thm2435.out
script_sha256:
  - 6ae536f85fb7c53b7cd8da86e993bd582a94cb28080258737b7299c939fced69
  - 2b7c03f102fa53bc42ba44dc73194731bdd40c3a1eec42da9cfb4b56529c83cd
output_sha256:
  - 5be089f3a1068518263bf4254ac1db9a1d8e2f9a0655b56a115a130e99907e62
  - 6a02977e4d1ccf5489d6e463e940e507e7ced84013742ab546bf8e6e8703a585
hash_basis: working-tree bytes (LF)
---

# THM-2435 -- essential top blockers carry a flat quotient-root section

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT
HOSTILE AUDIT.**

The three shapes left by THM-2431 and THM-2432 have one common
mechanism:

```text
low quotient blockers are averaged out by T_(7^(M+1))
  -> the parent image is larger than 4/7
  -> the six unit words cannot tile all those parents
  -> a fixed actual top blocker fills a positive family of gaps
  -> one equivariantly marked gap has a flat C_91 quotient spectrum. (1)
```

The last arrow is deliberately a quotient statement. At positive
septimal depth, forgetting the `7^M` ancestors is exactly the
information loss that prevents a physical `91`-unit current.

## 1. The three typed shapes

Retain THM-2427's primitive scalar cover and put

```text
M=nu_7(H)<nu_7(c_3),                 t=5,

k=#{j in {1,2}:nu_7(c_j)<=M},

b=#{j in {1,2}:nu_7(c_j)=M}.                                  (2)
```

After THM-2431 and THM-2432 the residual is exactly

```text
M=0:   (k,t,b,W)=(1,5,1,8), (2,5,2,9),

M>0:   (k,t,b,W)=(2,5,1,8).                                  (3)
```

Set

```text
N=7^(M+1),                    d=7^M,

u_0=H/d,                      u_i=q_i/d,

J={j in {1,2}:nu_7(c_j)=M},   v_j=c_j/d  for j in J.           (4)
```

The six `u_i` are units modulo `91`. Every exact-depth top blocker has

```text
v_j=13w_j,                    7 not|w_j.                        (5)
```

All identities below are almost everywhere. Endpoints and the finite
family of inherited exceptional pullbacks are removed once and for
all.

## 2. The sharpened parent image

Write `c_j=13C_j` and put

```text
A=D_(C_1)^c intersection D_(C_2)^c intersection D_(C_3)^c,

P=T_N(A).                                                       (6)
```

On a generic `N`-preimage fibre:

- if `nu_7(C_j)>M`, then `C_j=N e_j`, so the corresponding quotient
  blocker is constant on the whole fibre;
- if `nu_7(C_j)<=M`, that blocker occupies exactly `N/7` preimages.

There are only `k<=2` low quotient blockers. Their union therefore
cannot fill a fibre. Conversely, one dangerous high quotient blocker
removes the whole fibre. Hence the image is typed exactly by

```text
P
 =intersection_(nu_7(C_j)>M) D_(C_j/N)^c                 a.e.    (7)
```

and the union bound gives

```text
mu(P)>=1-(3-k)/7=(4+k)/7.                                      (8)
```

This is the first gain over the older `4/7` invoice.

## 3. The quotient stalk and its forced gaps

Use the quotient coordinate

```text
z=T_d x.                                                        (9)
```

For `Y in P`, the composed seven- and thirteen-root operation is the
common quotient stalk

```text
z_s=(Y+s)/91,                    s in Z/91Z.                     (10)
```

At `M=0`, this is the physical stalk. At `M>0`, it is the
`T_(7^M)` quotient stalk from THM-2430; this distinction is
load-bearing.

Put

```text
U=E_(u_0) union D_(u_1) union ... union D_(u_5),                 (11)

mathcal E={Y: every root in (10) belongs to U}.                  (12)
```

The guard contributes `26` incidences and the five ordinary words
contribute `13` each. Thus every cover in (12) is an exact one-fold
tiling. THM-2431's fixed repeated-pair argument gives

```text
mu(mathcal E)<=3/7.                                             (13)
```

Equations (8) and (13) imply

```text
mu(P minus mathcal E)>=(1+k)/7.                                 (14)
```

THM-2427 has `W>k`, so its top words alone cover every relevant
orbit. Consequently each quotient root missing from the six-unit
union (11) is covered by at least one **actual exact-depth top
blocker** `D_(v_j)`, `j in J`. Low blockers may shadow it but are not
needed for existence.

Define

```text
F^q={z:T_91 z in P, z notin U},

F_j^q=F^q intersection D_(v_j),

P_j={Y: the fibre over Y meets F_j^q}.                           (15)
```

The `P_j` cover `P\mathcal E`. Since `|J|=b`, one fixed label `j`
satisfies

```text
mu(P_j)>=(1+k)/(7b).                                            (16)
```

The three resulting invoices are

| shape | `mu(P_j)` | one marked quotient root |
|---|---:|---:|
| `M=0,(1,5,1,8)` | `>=2/7` | `>=2/637` |
| `M=0,(2,5,2,9)` | `>=3/14` | `>=3/1274` |
| `M>0,(2,5,1,8)` | `>=3/7` | `>=3/637` |

The fixed label is chosen before any Fourier character.

There is also an order-free blocker observable. For `z in F^q`, set

```text
omega(z)={j in J:z in D_(v_j)}.                                (16a)
```

This top-blocker word is nonempty. Its exact nonempty word strata
partition `F^q`; for `b=2` they are `{1}`, `{2}`, and `{1,2}`.
Moreover,

```text
mu(F^q)
 >=mu(P\mathcal E)/91
 >=(1+k)/637.                                                   (16b)
```

For either one-blocker shape there is only one word. For the
two-blocker shape, (16b) is at least `3/637`, so one of its three
words has mass at least `1/637`. Uniformly, therefore, one intrinsic
word `w` satisfies

```text
mu({z in F^q:omega(z)=w})>=1/637.                               (16c)
```

A singleton word is an edge candidate after a prescribed expiration;
the double word is a genuine fork/hyperedge in the sense of THM-2305.
No arbitrary orientation is imposed on the fork. At `M=0` these are
physical blocker words; at `M>0` they are quotient words until an
ancestry section is supplied.

## 4. The full labelled packet already sees every quotient character

Use CRT coordinates

```text
Z/91Z = F_7 x F_13.                                             (17)
```

Because `v_j=13w_j`, the blocker `D_(v_j)` selects exactly one
`F_7` sheet on every active generic fibre. On that sheet define

```text
S_j(Y)
 ={r in F_13: the CRT root (source_sheet_j(Y),r) lies in F_j^q}.
                                                                    (18)
```

The set `S_j(Y)` is nonempty for `Y in P_j`. It is proper: the
length-`2/7` guard meets three or four of the thirteen target roots
on the same source sheet, and those roots lie in `U`. Thus

```text
1<=|S_j(Y)|<=10.                                                 (19)
```

A `C_91` character restricts on the sheet to a phase times a
`C_13` character sum over `S_j(Y)`. For the trivial target character
the sum is `|S_j(Y)|`. For a nontrivial target character it cannot
vanish: by irreducibility of `Phi_13`, a rational `0/1` polynomial of
degree at most twelve vanishing at a primitive thirteenth root would
have to be `0` or `Phi_13`, contrary to (19). Therefore the direct
fixed-label packet `1_(F_j^q)` has positive polyphase energy in
**every one of the 91 quotient residue classes**.

This conclusion uses no chosen root and no signed termwise
factorization. The next section supplies a sharp uniform floor by
selecting one root equivariantly.

## 5. A translation-equivariant marked root

For a nonempty proper `S subset F_13`, define the cyclic binary word

```text
w_(S,a)=(1_S(a+r))_(r=0)^12.                                    (20)
```

The thirteen rotations are distinct. Indeed, a nontrivial
translation stabilizer in the prime cyclic group would make `S`
empty or full. Let

```text
sigma(S)=the unique a for which w_(S,a) is lexicographically maximal,

with 1>0.                                                        (21)
```

The maximal word begins with `1`, so `sigma(S) in S`, and

```text
sigma(S+t)=sigma(S)+t.                                          (22)
```

Pair `sigma(S_j(Y))` with the unique blocker source sheet in the CRT
coordinates (17). This selects one quotient root on every parent in
`P_j`. Call the resulting quotient packet `G^q`.

The construction is rational finite-step: away from finitely many
rational endpoints, all 91 root truth values are locally constant.
Fixing half-open representatives makes the selector literal without
changing any measure. It is equivariant under root translations, but
not asserted to be reflection- or affine-unit-equivariant.

Put

```text
rho_q=mu(G^q)=mu(P_j)/91.                                       (23)
```

For the normalized `C_91` root transform and every
`m in Z/91Z`,

```text
C_m(Y)
 =1/91 1_(P_j)(Y) exp(-2 pi i m sigma_91(Y)/91),                 (24)

|C_m(Y)|^2=1/91^2 1_(P_j)(Y).                                  (25)
```

The polyphase identity gives the exact flat bank

```text
sum_(n congruent m mod 91)|Ghat^q(n)|^2
 =mu(P_j)/91^2
 =rho_q/91                                                     (26)
```

for all 91 characters, including zero. Uniformly across the three
shapes,

```text
rho_q>=3/1274,

energy in each C_91 class >=3/115934,                            (27)

total energy in the 72 unit classes >=108/57967.                 (28)
```

The exact companion checks all `8190` nonempty proper `C_13` masks,
their `630` translation orbits, and all `106470` marker-equivariance
instances.

## 6. Physical pullback and the ancestry kernel

When `M=0`, `d=1` and `G^q` is a genuine physical packet. Thus both
`M=0` shapes have positive physical Fourier energy in every residue
class modulo `91`, with the floors (26)--(28). THM-2424's
endpoint-Prony lemma gives a finite surviving physical lift in every
nonzero class.

When `M>0`, define the physical pullback

```text
G^phys=(T_d)^(-1)(G^q),                    d=7^M.                 (29)
```

Each selected quotient root has `d` physical ancestors. On the
natural `C_(91d)` fibre, the normalized transform is

```text
C_q^phys(Y)=0                         unless d|q,

C_(dm)^phys(Y)
 =1/91 1_(P_j)(Y) exp(-2 pi i m sigma_91(Y)/91).                 (30)
```

Hence the 91 annihilator characters retain exactly the energy (26),
and all other characters vanish. Equivalently,

```text
Ghat^phys(n)=0 unless d|n,

Ghat^phys(dn)=Ghat^q(n).                                        (31)
```

Because `gcd(d m,91)=7` whenever `m` is a `91`-unit, the
positive-depth packet does **not** supply a physical `91`-unit
frequency. Selecting one physical ancestor would require a new
ancestry section and would break the present quotient-gauge theorem.

## 7. The linear row-defect carrier

The same quotient stalk has a complementary signed carrier. Use CRT
coordinates

```text
(h,r) in F_13 times F_7
```

for the root index `s` in (10). For a parent `Y` and target row `h`,
let `m_(Y,h)(r)` be the multiplicity at `(h,r)` among the guard and
five ordinary unit masks in (11). Every ordinary word meets each
`F_13` row once and the guard meets it twice. Hence

```text
sum_(r in F_7)m_(Y,h)(r)=7.                                   (31a)
```

Define the full two-dimensional incidence defect by

```text
d_Y(h,r)=1-m_(Y,h)(r).                                        (31b)
```

Every row is integral and has sum zero. Its positive support is
exactly the set of unit-cover holes in that row, so (by the
top-blocker cover preceding (15)) it lies in the one or two
top-blocker source sheets. If the row contains a point of `F^q`, its
defect is nonzero.

For an integral zero-sum vector `e` on `F_7`, put

```text
E_7(e)=sum_(r in F_7)(e(r+1)-e(r))^2.                          (31c)
```

Every nonzero such vector has the sharp bound

```text
E_7(e)>=4.                                                     (31d)
```

To prove it, note that the seven edge differences are integral and
sum to zero. Energy one is impossible. Energy two would have one
`+1` and one `-1`; then `e` has two consecutive integer levels `a`
and `a+1`, occurring `7-m` and `m` times for some `1<=m<=6`.
The zero-sum law `7a+m=0` is impossible. Energy three would require
three signed unit differences to sum to zero, also impossible.
Since the energy is integral, (31d) follows. A separated
source-to-duplicate arrow attains equality.

Every nontrivial discrete row colour is also nonzero. If a primitive
seventh root `zeta` and `beta=1,...,6` satisfied

```text
sum_(r in F_7)e(r)zeta^(beta r)=0,                             (31e)
```

then the rational coefficient polynomial of degree at most six would
be divisible by `Phi_7=1+X+...+X^6`. All its coefficients would be
equal, and their zero sum would make `e=0`, a contradiction.

The normalization of the integrated energy is load-bearing. Define
the **unnormalized** row sum

```text
mathscr D(Y)=sum_(h in F_13)E_7(d_Y(h,-)).                      (31f)
```

Each `Y in P\mathcal E` has at least one essential row. Equations
(14) and (31d) give the typed and uniform floors

```text
integral_P mathscr D(Y)dY
 >=4(1+k)/7
 >=8/7.                                                        (31g)
```

If the thirteen rows were averaged rather than summed, the uniform
floor would be `8/91`.

When `b=1`, an essential row has exactly one hole at the blocker
source `r_0`; all six other sites are unit-covered. Equation (31a)
then forces exactly one duplicate `r_1`, and

```text
d_Y(h,-)=e_(r_0)-e_(r_1),                  r_0!=r_1.             (31h)
```

This is a canonical row arrow `r_0 -> r_1`. Its energy is six for
adjacent residues and four otherwise; among all `42` arrows the
histogram is `6:14, 4:28`. For `b=2`, holes remain supported on the
one or two top sources (which may coincide), retaining the
two-source hole-to-duplicate flow.

Finally, the **full** defect detects exactly the target-unit fallback.
With any fixed CRT character convention, set

```text
dtilde_Y(alpha,beta)
 =sum_(h in F_13)sum_(r in F_7)
   d_Y(h,r)zeta_13^(-alpha h)zeta_7^(-beta r).                  (31i)
```

All mixed coefficients with `alpha!=0` and `beta!=0` vanish iff
two-way Fourier decomposition gives `d_Y(h,r)=a(h)+b(r)`. The
row-sum law then says `7a(h)+sum_r b(r)=0`, so `a` is constant and
can be absorbed into `b`. Conversely, an `h`-independent defect has
no nonzero-`alpha` coefficient. Therefore

```text
d_Y varies with h
 iff dtilde_Y(alpha,beta)!=0
     for some alpha!=0, beta!=0.                               (31j)
```

Such a mixed pair is nonzero modulo both `13` and `7`, hence is a
unit character modulo `91`. Thus the only quotient-stalk fallback
without a `91`-unit mixed defect mode is a vertical,
row-independent defect. This signed row DFT differs from the Boolean
packets of Sections 4--5: it is linear in the six unit incidences and
gives pointwise nonvanishing, but it is not by itself an integrated
physical coefficient or endpoint current.

## 8. What is canonical, and what is not

The direct packet `F_j^q` is a Boolean factor-repair packet carrying
one fixed actual top-blocker label. The flat quantitative section
`G^q` is more selective: it is a finite Boolean function of the
thirteen root translates of `F_j^q`. Current THM-2334/THM-2410
single-copy probe algebra is not proved closed under this nonlinear
thirteen-copy operation. Therefore (26) is a derived root-section
Fourier bank, not a canonical Abel endpoint or relation current.

The owner typing is:

- in `M=0,(1,5,1,8)`, the unique top label is an exclusive blocker
  owner on the gap;
- in `M=0,(2,5,2,9)`, one has only exclusive owner versus
  double-top collision;
- in `M>0,(2,5,1,8)`, one has only exclusive top owner versus
  top/low-blocker collision.

All coverage is punctured and almost everywhere. MISTAKE-273's exact
open-mask seam shows that no literal endpoint owner follows.
Integrated complex coefficients can also cancel while the energy
bank (26) remains flat.

Thus THM-2435 does not yet produce a signed handoff, terminal word,
common endpoint phase, physical positive-depth `91`-unit relation,
valuation-shape exclusion, scalar-row decrement, or proof of
LRC(14). It isolates two orthogonal missing coordinates:

```text
M=0:  thirteen-copy Boolean section -> lawful common-endpoint probe,

M>0:  quotient root -> physical ancestry section,

both: energy -> retained complex current phase.                  (32)
```

## 9. Exact companion

Run

```text
python3 04-computation/lrc14_top_blocker_essential_parent_thm2435.py
python3 -O 04-computation/lrc14_top_blocker_essential_parent_thm2435.py
python3 04-computation/lrc14_top_blocker_essential_carrier_thm2435.py
python3 -O 04-computation/lrc14_top_blocker_essential_carrier_thm2435.py
```

The two dependency-free companions:

- exhausts all `352947` labelled seven-bin unit profiles;
- verifies the three sharpened `(k,b)` measure invoices;
- exhausts all `8190` nonconstant `C_13` masks, `630` translation
  orbits, and `106470` marker-equivariance checks;
- performs `745290` exact cyclotomic/character checks without
  floating point;
- verifies the exact `C_91` flat-energy floors; and
- checks the positive-depth ancestry annihilator and the loss of
  physical `91`-unit residues;
- checks generic ordinary, guard, and blocker row-incidence laws;
- exhausts all `1716` seven-incidence row multiplicity vectors,
  obtaining the sharp energy floor four and the one-source arrow
  histogram `4:28, 6:14`;
- verifies the sharpened `2/7`, `2/637`, `3/1274`, `1/637`, and
  `8/7` carrier invoices and the exact `91`-branch Jacobian; and
- checks explicit one-/two-source punctured repairs while retaining
  THM-2427's local exact tiling as a zero-essential hostile.

Normal and optimized runs must reproduce

```text
05-knowledge/results/lrc14_top_blocker_essential_parent_thm2435.out
05-knowledge/results/lrc14_top_blocker_essential_carrier_thm2435.out
```

byte-for-byte, respectively.

## 10. Independent audit

Independent audits have accepted the sharpened parent image, the
fixed-label measure invoices, the lexicographic selector, and the
quotient/physical ancestry distinction. Final replay and immutable
hash audit are pending. QED conditional on that final audit.
