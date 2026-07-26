---
id: THM-2383
title: "Polarized complete-subcube Gram tomography"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For maps
  from a labelled Boolean cube to a complex Hilbert space, the
  polarized complete-subcube Dirichlet bank has the exact 0/2 Walsh
  multiplier. Labelled Mobius inversion recovers the cross-Gram
  coefficient on every nonempty Walsh support. Two squared-energy
  quadratures recover each complex cross bank. A reference family
  recovers a coefficient exactly on the span of its same-support
  reference vectors, and spans are necessary for universal recovery.
  Constants, coordinate labels, a complex quadrature, and an absolute
  reference are all sharp: explicit hostiles lose each one. THM-2370's
  complemented clone cube is the no-reference terminal-orientation
  hostile. THM-2380 remains the separate physical cross-word twist
  service. No canonical LRC current, knot realization, Gaussian scalar
  moment consequence, ledger decrement, or conjecture closure follows.
source: codex-2026-07-25-polarized-subcube-gram-tomography
depends_on: []
related:
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2374-binary-allocation-complete-subcube-dirichlet-spectrum
  - THM-2375-gaussian-angular-complete-line-charge-tomography
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
script: 04-computation/polarized_complete_subcube_gram_tomography_thm2383.py
output: 05-knowledge/results/polarized_complete_subcube_gram_tomography_thm2383.out
script_sha256: 81809764c5fdd273783be9688f3c724f8fd4b4a52ce39f13988306c5042add33
output_sha256: 1dafb6feff2a3450af748c5338376b1ece3b026242b90dc43159707371ed8fef
hash_basis: working-tree bytes (LF)
---

# THM-2383 -- polarized complete-subcube Gram tomography

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2374 recovers the squared norm of every Boolean ANOVA component,
but its sharp THM-2370 clone hostile loses the component signs and the
absolute terminal orientation. The missing abstract datum is exactly a
polarized reference:

```text
squared complete-subcube bank
  -> support-labelled coefficient norms;

polarized complete-subcube bank against spanning references
  -> support-labelled complex coefficients.                    (1)
```

This theorem proves (1), identifies the minimal energy measurements
which produce the polarization, and separates that algebra from the
physical reference/twist service still missing in LRC(14).

## 1. Boolean and Hilbert conventions

Let `S` be a finite labelled set and put

```text
X=F_2^S.
```

For `A subset S`, write

```text
chi_A(x)=(-1)^|A intersect x|.                        (2)
```

Let `H` be a complex Hilbert space whose inner product is linear in the
first slot. For `F:X->H`, use normalized Walsh coefficients

```text
Fhat(A)=2^(-|S|) sum_(x in X) F(x)chi_A(x),

F(x)=sum_(A subset S) Fhat(A)chi_A(x).               (3)
```

For `U subset S`, identify `F_2^U` with the coordinate subgroup of
`X`. For maps `F,R:X->H`, define the cross-Dirichlet form

```text
B_U(F,R)
 =E_(x in X,h in F_2^U)
   <F(x+h)-F(x),R(x+h)-R(x)>.                        (4)
```

The squared bank is

```text
D_U(F)=B_U(F,F)>=0.                                  (5)
```

All expectations in (3)--(5) are uniform on the displayed **labelled**
sets.

## 2. The polarized 0/2 multiplier

For every `U subset S`,

```text
B_U(F,R)
 =2 sum_(A:A intersect U!=empty)
      <Fhat(A),Rhat(A)>.                             (6)
```

### Proof

Expand both differences in Walsh characters. Orthogonality in `x`
kills unequal supports. On the surviving support `A`, the average over
`h in F_2^U` is

```text
E_h (chi_A(h)-1)^2
 =0                              if A intersect U=empty,
 =2                              otherwise.          (7)
```

This proves (6). The argument is a finite sum and works unchanged for
any real or complex Hilbert space. QED.

Equation (6) is not merely a positivity statement. It retains the full
complex Gram coefficient

```text
g_A(F,R)=<Fhat(A),Rhat(A)>                           (8)
```

inside a triangular labelled support transform.

## 3. Exact labelled Mobius inversion

Put `SminusV=S\V` and define

```text
Q(V)
 =(B_S(F,R)-B_(SminusV)(F,R))/2

 =sum_(empty!=A subset V) g_A(F,R).                  (9)
```

For every nonempty `A subset S`, Boolean Mobius inversion gives

```text
g_A(F,R)
 =sum_(V subset A)(-1)^(|A|-|V|) Q(V).              (10)
```

Equivalently,

```text
g_A(F,R)
 =1/2 sum_(V subset A)(-1)^(|A|-|V|)
      (B_S(F,R)-B_(SminusV)(F,R)).                  (11)
```

Thus the complete labelled bank `{B_U:U subset S}` recovers the
cross-Gram coefficient on every exact nonempty Walsh support.

The empty support is absent from every difference in (4). No choice of
subcubes or references recovers `Fhat(empty)` without a separate mean
measurement.

## 4. Two squared-energy quadratures recover the complex bank

The cross form need not be given directly. The four real squared banks

```text
D_U(F),  D_U(R),  D_U(F+R),  D_U(F+iR)              (12)
```

give

```text
Re B_U(F,R)
 =[D_U(F+R)-D_U(F)-D_U(R)]/2,

Im B_U(F,R)
 =[D_U(F+iR)-D_U(F)-D_U(R)]/2.                      (13)
```

The sign in the second line uses the convention that the inner product
is linear in its first slot. Combining (13) with (9)--(11) reconstructs
every complex cross-Gram coefficient using squared energies only.

The singleton banks `D_U(F)` and `D_U(R)` need be measured only once.
For each reference, one ordinary union and one genuine complex
quadrature then suffice. This is the Boolean-support analogue of the
target-current polarization in THM-2380; it does not construct a
physical phase twist. In particular, the Hilbert-space vector
`F+iR` need not be a lawful indicator packet or positive physical
operation.

## 5. Spanning references are exactly the missing algebraic sidecar

Let

```text
R^1,...,R^m:X->H
```

be labelled references, and fix a nonempty support `A`. Put

```text
r_j=Rhat^j(A),

V_A=span_C{r_1,...,r_m}.                             (14)
```

Equations (10) and (13) recover

```text
g_j=<Fhat(A),r_j>                  for every j.       (15)
```

These numbers determine the orthogonal projection of `Fhat(A)` onto
`V_A`. In particular, if `Fhat(A) in V_A`, they determine `Fhat(A)`
exactly.

For an explicit formula, choose a basis `r_1,...,r_d` of `V_A` and
write

```text
Fhat(A)=sum_(k=1)^d c_k r_k.
```

With

```text
G_(jk)=<r_k,r_j>,

g_j=<Fhat(A),r_j>,                                  (16)
```

the Hermitian Gram matrix is invertible and

```text
c=G^(-1)g.                                          (17)
```

Consequently:

- references spanning `H` recover every nonconstant coefficient of
  every `H`-valued cube;
- supportwise references spanning a known coefficient subspace recover
  every coefficient in that subspace; and
- without the span hypothesis, the cross bank recovers only a
  projection. The self bank also gives the residual norm, but not its
  direction or phase.

Quantitatively, if the supportwise analysis map obeys the lower frame
bound

```text
a_A ||v||^2
 <=sum_j |<v,r_j>|^2                 for v in V_A,   (17a)
```

then

```text
||v||<=a_A^(-1/2)(sum_j |<v,r_j>|^2)^(1/2).         (17b)
```

The constant is sharp already for one reference
`r=sqrt(a_A)e`. Thus reference conditioning, not merely formal span,
controls stable tomography.

The span condition is necessary for universal recovery. If
`V_A!=H`, choose `0!=u in V_A^perp`. The two maps having only Walsh
coefficient `+u` or `-u` at support `A` have identical self-Dirichlet
banks and zero cross banks against every reference, but their
coefficients differ. Over a complex Hilbert space, the whole phase
orbit `e^(i theta)u` is invisible to the reference bank.

The reconstruction presumes that the reference vectors are known in
the ambient Hilbert space. If only their Gram matrix is known, the
coefficient is determined only up to the corresponding common unitary.
No bank here measures a cross-support product
`<Fhat(A),Rhat(B)>` with `A!=B`.

## 6. Four sharp information boundaries

### 6.1 The constant mode

Replacing `F(x)` by `F(x)+c` for a fixed `c in H` changes neither (4)
nor (5). A separate mean/anchor is necessary and sufficient for the
empty support.

### 6.2 One real union

Take `H=C`, `R(x)=chi_A(x)`, and

```text
F_+(x)= i chi_A(x),

F_-(x)=-i chi_A(x).                                 (18)
```

The two maps have the same self bank and the same ordinary-union bank
`D_U(F+R)` for every `U`, because both real Gram parts vanish. Their
imaginary Gram parts are opposite. The `F+iR` quadrature in (13)
separates them. Thus one real union is insufficient over `C`.

### 6.3 Forgotten coordinate labels

For distinct `s,t in S`, the cubes

```text
F_s(x)=u chi_{s}(x),

F_t(x)=u chi_{t}(x)                                 (19)
```

have the same bank after coordinate relabelling and the same profile
after grouping subcubes only by `|U|`. The labelled bank distinguishes
them immediately. Exact support recovery therefore requires labelled
subcubes.

### 6.4 No absolute reference: the clone terminal hostile

Let `|S|=n`, fix `0!=v in H`, and put

```text
F(x)       =(1-|x|/n)v,

Ftilde(x)  =(|x|/n)v=v-F(x).                        (20)
```

Then

```text
D_U(F)=D_U(Ftilde)                    for every U,   (21)
```

because every difference changes sign. Yet all nonconstant Walsh
coefficients change sign and

```text
F(1,...,1)=0,

Ftilde(1,...,1)=v.                                  (22)
```

This is the Hilbert-valued core of THM-2370's exact Boolean clone
hostile. The complete squared spectrum cannot decide which terminal is
null. A spanning **oriented** reference does. Without one, any
homogeneous quadratic **difference** bank invariant under
`F -> v-F` remains unable to orient the terminal. The uncentered scalar
endpoint table `x -> ||F(x)||^2` is not such a bank and can distinguish
the two displayed terminal values.

Finally, the factor `2` in (6) uses the complete subgroup, including
`h=0`. If the zero shift is deleted, a nontrivial character instead has
multiplier

```text
2^(|U|+1)/(2^|U|-1),
```

so the displayed complete-subcube normalization is load-bearing.

## 7. Scope and cross-frontier use

The theorem supplies an algebraic target for several live debts:

- For THM-2374, it identifies the exact additional measurements needed
  to turn ANOVA norm recovery into signed/complex coefficient recovery.
  It does not realize an abstract allocation table by knots or identify
  a minimizing owner.
- For THM-2380, it says why one endpoint-matched cross-word quadrature
  is sufficient once a reference current exists. It does not produce
  that physical LRC observable, align endpoint gauges, or transfer a
  derived pair current into the canonical word current.
- For THM-2375, it is the finite Boolean counterpart of labelled
  isotypic Gram tomography. These Hermitian measurements are richer
  than scalar Gaussian moments and give no new NC2 implication.

No scalar LRC profile is excluded; the ledger remains `165`.
LRC(14), the knot realization problems, and the Gaussian moment
conjecture do not follow from this theorem.

## 8. Exact companion

The dependency-free `Fraction`/Gaussian-rational companion:

- verifies the direct cross-Dirichlet multiplier on `48` labelled
  reference/subcube cells in rank four;
- checks `96` real and imaginary polarization identities;
- recovers `45` exact support-labelled cross-Gram coefficients by
  Mobius inversion;
- reconstructs `15` three-dimensional coefficients from nonorthogonal
  spanning references by exact Gram elimination;
- checks the proper-span, one-quadrature, constant, and unlabelled
  hostiles; and
- verifies all `32` squared subcube identities in the rank-five clone
  terminal hostile.

Run

```bash
python3 04-computation/polarized_complete_subcube_gram_tomography_thm2383.py
python3 -O 04-computation/polarized_complete_subcube_gram_tomography_thm2383.py
```

Both transcripts must match

```text
05-knowledge/results/polarized_complete_subcube_gram_tomography_thm2383.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Two independent hostile audits rederived the `0/2` multiplier, labelled
Mobius inversion, both linear-first polarization signs, Gram solve,
supportwise span criterion, sharp frame bound, and every information
boundary. The first audit caught and repaired an overbroad clone claim:
only homogeneous quadratic difference banks invariant under
`F -> v-F` share that hostile, not the uncentered endpoint norm table.
Both audits replayed normal and optimized execution against the stored
transcript and independently reproduced the recorded LF hashes. QED.
