---
id: THM-2265
title: "Balanced-cut spectral coupling and bounded resonance"
status: >
  SUPERSEDED / WEAKER DUPLICATE OF THM-2145. The proof and exact checks are
  valid, but THM-2145 already proves the same abstract positive crossing
  mechanism with universal height 584 and defect-six height 298, while
  THM-2166 gives the stronger low-carry/core-sparse form. Retained only for
  its independent derivation and anisotropic bandwidth formula; do not use
  this theorem as the canonical crossing dependency.
source: codex-2026-07-25-balanced-cut-spectral-coupling
depends_on:
  - THM-1166
  - THM-1221
related:
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-2085-explicit-height-57-rank-seven-selberg-gate
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
  - THM-2166-hybrid-core-smoothing-low-carry-crossing
script: 04-computation/lrc14_balanced_cut_spectral_coupling_referee_codex_20260725.py
output: 05-knowledge/results/lrc14_balanced_cut_spectral_coupling_referee_codex_20260725.out
script_sha256: 8d8f634c68da777d20e870758d0d9980e43c71457129c5dad40cd16f76fca20d
output_sha256: bbc663b973c2f8ce75a4e8923953506e6db1a48102e44560441a3b4ede0c9bf1
hash_basis: working-tree bytes (LF)
---

> **SUPERSEDED / WEAKER DUPLICATE.** THM-2145, which predates this file,
> already proves the abstract common-frequency crossing lemma, universal
> height `584`, and defect-six height `298` with carry at most `20860`.
> THM-2166 further improves the defect-six carrier to `|nu|<=708`, far
> height `298`, and a core representation of support at most two and height
> `57`. The derivation below remains correct and supplies a convenient
> anisotropic bandwidth formula, but no current theorem should depend on its
> weaker headline bounds.

# THM-2265 -- balanced-cut spectral coupling and bounded resonance

Work on `T=R/Z`. For `0<h<1/2`, let

```text
chi_h(x)=1_[h,1-h](x),                                  (1)
G_A={t:chi_h(at)=1 for every a in A}.                   (2)
```

Endpoint choices do not change any integral below.

## 1. A positive band-limited approximation

For `N>=2`, let

```text
F_N(t)=(1/N)*(sin(pi N t)/sin(pi t))^2,
J_N(t)=F_N(t)^2 / integral_T F_N^2.                     (3)
```

Then `J_N>=0`, `integral J_N=1`, and its Fourier support is contained in
`[-2N+2,2N-2]`. Put

```text
q_N=J_N*chi_h.                                          (4)
```

Thus `0<=q_N<=1`, it has the same Fourier support, and

```text
||q_N-chi_h||_1<3/(2N).                                 (5)
```

### Proof of (5)

Represent `T` by `[-1/2,1/2]` and write `|t|` for centered distance. Parseval
gives

```text
integral F_N^2=(2N^2+1)/(3N)>2N/3.                      (6)
```

The elementary bounds

```text
F_N(t)<=min(N,1/(4N|t|^2))                              (7)
```

follow from `|sin(pi t)|>=2|t|`. Splitting at `1/(2N)` gives

```text
integral |t|F_N(t)^2 dt
 <=1/2-1/(4N^2)<1/2.                                   (8)
```

Hence `integral |t|J_N(t)dt<3/(4N)`. Translation of one circle interval by
`t` changes its indicator in `L1` by at most `2|t|`; for a long or wrapping
interval this is the same assertion applied to its complement. Averaging this
inequality against `J_N` proves (5). QED.

## 2. Positive spectral coupling

Partition a finite labelled speed set into nonempty blocks `F` and `E`.
Choose one integer bandwidth `N_i>=2` for each speed and set

```text
eta_i=3/(2N_i),
P_F(t)=product_(i in F) q_(N_i)(v_i t),
P_E(t)=product_(i in E) q_(N_i)(v_i t).                 (9)
```

Product telescoping, using that every factor lies in `[0,1]`, gives

```text
integral P_F >=mu(G_F)-eta_F,
integral P_E >=mu(G_E)-eta_E,
integral P_F P_E <=mu(G_(F union E))+eta_F+eta_E,       (10)
```

where `eta_F=sum_(i in F)eta_i` and similarly for `E`.

Suppose the two lower bounds in (10) are positive and

```text
(mu(G_F)-eta_F)(mu(G_E)-eta_E)
  >mu(G_(F union E))+eta_F+eta_E.                       (11)
```

Then there are integers `c_i`, with

```text
|c_i|<=2N_i-2,                                         (12)
```

such that

```text
sum_(i in F)c_i v_i + sum_(j in E)c_j v_j=0,           (13)
sum_(i in F)c_i v_i=-sum_(j in E)c_j v_j!=0.           (14)
```

### Proof

Write `supp_hat` for Fourier support. If

```text
supp_hat(P_F) intersect (-supp_hat(P_E))={0},           (15a)
```

Fourier orthogonality would give

```text
integral P_F P_E=(integral P_F)(integral P_E).          (15)
```

Equations (10)--(11) contradict (15). Therefore some `k!=0` occurs in the
support of `P_F` and `-k` occurs in the support of `P_E`. Expanding each block
product expresses these two frequencies as

```text
k=sum_(i in F)c_i v_i,       -k=sum_(j in E)c_j v_j,    (16)
```

with (12), proving (13)--(14). Notice that this argument does not assume
independence of the speeds: failure of independence is exactly the conclusion.
QED.

## 3. Every 6+7 cut of a strict LRC(14) counterexample

Now set `h=1/14`, take thirteen distinct positive integer speeds, and assume

```text
G_all is empty, equivalently M(v)<1/14.                 (17)
```

For any six-speed block `F`, the universal THM-1166 pair floor `1/91` and
the spanning-tree Hunter inequality give

```text
mu(G_F)>=1-6/7+5/91=18/91.                              (18)
```

For its seven-speed complement `E`, THM-1221 gives

```text
mu(G_E)>=15/154.                                        (19)
```

Take every `N_i=1162`. The strict margin in (11) is

```text
(18/91-9/1162)*(15/154-21/(2*1162))
  -39/(2*1162)
 =465/35106344>0.                                       (20)
```

The spectral support bound is `2N-2=2322`. Consequently, for **every**
labelled partition into six and seven coordinates, there is a relation

```text
sum_(f in F)a_f f+sum_(e in E)b_e e=0,
|a_f|,|b_e|<=2322,                                     (21)
```

whose two displayed block sums are both nonzero. This is stronger than merely
finding a bounded relation somewhere: no chosen balanced cut can separate all
bounded resonances.

Under the estimates (18)--(19), `1162` is the first successful common
bandwidth: the `N=1161` margin is

```text
-2293/699620922<0.                                      (22)
```

This is optimal only for these safe-mass and squared-Fejer error bounds.

## 4. An anisotropic relation budget

If bandwidth `N_F` is used on all six coordinates of `F` and `N_E` on all
seven coordinates of `E`, the exact sufficient condition is

```text
(18/91-9/N_F)*(15/154-21/(2N_E))
  >9/N_F+21/(2N_E).                                    (23)
```

It yields coefficient bounds `2N_F-2` on `F` and `2N_E-2` on `E`. This is an
independently derived relation budget, superseded in its headline consequences
by THM-2145/2166 and parallel to THM-2144's Selberg--Kraft budget:
one may trade a smaller coefficient atlas on one side of a cut for a larger
carry range on the other. For example,

```text
(N_F,N_E)=(600,4427)
```

has positive exact margin `9/326526200` and gives height bounds `(1198,8852)`.
No claim of optimal bandwidth allocation is made.

## 5. Defect-six core carry

In a defect-six AP-core row, let `E` be the seven retained elements of
`{1,...,13}` and `F` the six replacements. Equation (21) gives the nonzero
core carry

```text
C=sum_(e in E)b_e e=-sum_(f in F)a_f f,                 (24)
0<|C|<=2322*sum_(e in E)e<=2322*70=162540.              (25)
```

Thus the six replacement values lie on a finite union of bounded-coefficient
affine hyperplanes with bounded right side. If all replacement speeds exceed
`162540`, the nonzero far coefficients cannot all have the same sign, so the
far block contains an internal cancellation as well.

This reduces a raw six-parameter branch by one affine equation. It does
**not** prove that enumerating all coefficient templates is practical: the
height-`2322` atlas is enormous, coefficients may vanish, and (24) does not
choose a canonical relation. A relation-carry dynamic program must still
prove compression or pruning beyond the existence theorem.

## 6. Verification and transfer boundary

The dependency-free referee checks the Fejer normalization and moment bound,
both safe-mass floors, the exact margins (20)--(23), and (25). Normal and
optimized runs are intended to be byte-identical:

```text
python 04-computation/lrc14_balanced_cut_spectral_coupling_referee_codex_20260725.py
python -O 04-computation/lrc14_balanced_cut_spectral_coupling_referee_codex_20260725.py
```

The mechanism is positivity plus finite Fourier support. It does not transfer
unchanged to signed Wick channels, infinite-memory carry languages, or
infinite fibres: in those settings a shared support frequency can cancel, or
the finite bandwidth may forget the continuation sidecar.
