---
id: THM-3961
title: "Arbitrary-q hidden repetition exactly detects monogenic conductor debt"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For an irreducible cubic
  T^3-3PT-q(P,t), remove only the forced h^2 factor from
  q(h^2,t)-2h^3 when q is exactly once divisible by P. The resulting
  hidden polynomial is squarefree if and only if the cubic hypersurface is
  normal. If P^2 divides q, the zero section is generically singular and
  the hypersurface is nonnormal. Thus every normal member is excluded from
  a same-function-field planar Keller chart by its global monogenic
  different. The only unresolved members of this arbitrary-q monogenic
  grammar are the two explicitly typed nonnormal conductor-debt loci; this
  does not close JC(2).
source: jc-degree6-one-place / post-THM-3960 arbitrary-q extension, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift, 2026-08-24). The audit
  reconstructed the ramification/Jacobian dictionary, the height-one prime
  attached to every repeated adjusted factor, the h=0 forced-square seam,
  the P^2 zero-section failure, the vertical/common-factor boundary cases,
  and the R1+S2 normality converse. It also checked the irreducible-domain
  scope and the normalization-form Zariski Main/global-different unit
  obstruction. Normal and optimized 51-gate runs match the frozen transcript
  after canonical LF normalization on Windows; hashes and documentation
  checks pass.
depends_on: []
related:
  - THM-3960-natural-one-parameter-cubic-normal-monogenic-closure
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
script: 04-computation/jc2_arbitrary_q_hidden_repetition_normality_thm3961.py
output: 05-knowledge/results/jc2_arbitrary_q_hidden_repetition_normality_thm3961.out
script_sha256: e23bc41290cc397df5de5caa7bfc2cf24d911fd21862d455f6b6343b03b34064
output_sha256: 965d720d954670a8c4662adfe6c3940c1202199aa011dfbdc66fe3ce40c23764
semantic_sha256: 1993c77cd4dc081abb018ace8f177bc80ada708f19622a7a4c15c29fb385564b
hash_basis: raw LF bytes
---

# THM-3961 -- adjusted hidden squarefreeness is the exact normality gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of
characteristic zero. Let

```text
R=k[P,t],
F=T^3-3PT-q(P,t),                     q in k[P,t],        (1)
A=R[T]/(F),                           X=Spec A.
```

Assume throughout that `F` is irreducible over `k(P,t)`. Thus `A` is a
domain, finite free of rank three over `R`, with basis `(1,T,T^2)`. The
irreducibility hypothesis is a genuine first gate: this theorem does not
classify all polynomial roots of an arbitrary `q`. It classifies normality
and the Keller obstruction exactly after that domain gate has passed.

Write the unique `P`-adic expansion

```text
q=q0(t)+P q1(t)+P^2 q2(P,t)                            (2)
```

and define the hidden ramification polynomial

```text
K(h,t)=q(h^2,t)-2h^3.                                  (3)
```

There are exactly three cases.

```text
(I)   q0 != 0:                 L=K;
(II)  q0=0, q1 != 0:          K=h^2 L,
                              L=q1(t)-2h+h^2q2(h^2,t);
(III) q0=q1=0:                P^2 divides q.            (4)
```

In `(I)` and `(II)`, `X` is normal if and only if `L` is squarefree in
`k[h,t]`. In `(III)`, the height-one zero section `(P,T)` is singular, so
`X` is not normal. The forced factor `h^2` in `(II)` is therefore harmless;
every other repeated hidden factor is genuine codimension-one normalization
debt.

## 1. The ramification scheme and its exact Jacobian dictionary

The relative different generator and the other two partial derivatives are

```text
F_T=3(T^2-P),
F_P=-3T-q_P,
F_t=-q_t.                                                (5)
```

On the ramification scheme put

```text
P=h^2,                           T=-h.                   (6)
```

Then `F=-K`, and direct differentiation gives the identities

```text
K_h=-2h F_P|_(6),                K_t=-F_t|_(6).          (7)
```

Equivalently,

```text
A/(T^2-P) = k[h,t]/(K),          h=-T.                   (8)
```

Thus, away from `h=0`, the singular scheme of `X` on its necessarily
ramified locus is exactly the singular scheme of the plane curve `K=0`.
In case `(II)`, after restricting away from `h=0`, this is equivalently the
singular scheme of `L=0`, because there `K=h^2L` and `h` is a unit.

At `h=0`, or equivalently `P=T=0`, the four hypersurface equations
`F=F_T=F_P=F_t=0` reduce exactly to

```text
q0(t)=q1(t)=q0'(t)=0.                                  (9)
```

This separate row is why the forced `h^2` in `(II)` must not be mistaken
for a nonnormal divisor.

## 2. Every repeated nonzero hidden factor is codimension-one singular

In cases `(I)` and `(II)`, `h` does not divide `L`: in `(I)`,
`L(0,t)=q0(t)` is not the zero polynomial, and in `(II)`,
`L(0,t)=q1(t)` is not the zero polynomial. Suppose an irreducible
`M(h,t) != h` occurs in `L` with multiplicity at least two. The ideal

```text
p_M=(P-T^2, M(-T,t)) in A                              (10)
```

is prime: after imposing `P=T^2`, its quotient is the domain
`k[T,t]/(M(-T,t))`, since `M(-T,t)` divides the resulting equation
`K(-T,t)`. It has height one. At its generic point `h=-T` is nonzero;
the repeated-factor identities `L=L_h=L_t=0`, together with `(5)--(7)`,
make every partial derivative of `F` vanish. Hence `A_(p_M)` is not regular.
The hypersurface violates `R1` and is not normal.

This is scheme-theoretic, not merely a reduced-support count. A repeated
nonzero hidden component supplies an actual height-one component of the
Jacobian scheme and hence an actual conductor divisor after normalization.

## 3. Squarefree adjusted hidden curve implies normality

Conversely, suppose `L` is squarefree. A reduced plane curve in
characteristic zero has only finitely many singular points: otherwise
`L,L_h,L_t` would share an irreducible curve component, forcing that
component to divide `L` twice. By `(7)`, all singular points of `X` with
`h != 0` therefore form a finite set.

At `h=0`, equation `(9)` again gives only finitely many points. In case
`(I)`, they lie among the multiple roots of the nonzero one-variable
polynomial `q0`; the additional equation `q1=0` can only remove points. In
case `(II)`, they lie among the roots of the nonzero polynomial `q1`.
Consequently the entire singular locus of the two-dimensional domain `X`
is finite. The hypersurface `A` is Cohen--Macaulay and hence `S2`; the
finite singular locus gives `R1`. Serre's criterion proves

```text
in cases (I),(II):       A normal  iff  L squarefree.     (11)
```

There is no hidden vertical exception. In case `(I)`, the coefficient of
`h^3` in `L=K` is `-2`; in case `(II)`, the coefficient of `h` in `L` is
`-2`. Hence `L` is primitive over `k[t]` and has no vertical factor.
Common factors or common zeros of `q0` and `q1` can create only the finite
exceptional set `(9)`, unless both rows vanish identically, which is the
next case.

## 4. Exact double-`P` debt

If `q0=q1=0`, then `q=P^2q2(P,t)` and

```text
K=h^3(h q2(h^2,t)-2).                                  (12)
```

The quotient `A/(P,T)=k[t]`, so

```text
E0=(P,T)                                                 (13)
```

is a height-one prime. Every equation in `(5)` vanishes along `E0`, as do
`F` and `F_t`, because `P^2` divides `q`. Thus the whole generic point of
`E0` is singular. The domain `A` fails `R1` and is nonnormal.

This is the second, and only other, nonnormality mechanism in the stated
universe. Notice that `(12)` has a forced `h^3`, rather than the harmless
`h^2` of case `(II)`; its extra order is exactly the lost linear `P` row.

## 5. Every normal case is excluded by the global different

Suppose `(I)` or `(II)` holds and `L` is squarefree. By `(11)`, `X` is the
finite normalization of `A2_(P,t)` in its cubic function field. If a
same-function-field planar Keller map existed, normalization-form Zariski
Main would give an open immersion

```text
j:U isomorphic to A2 -> X                                (14)
```

on which `X -> A2_(P,t)` is etale. The global monogenic presentation gives

```text
Omega_(A/R)=A dT/(delta dT),
delta=F_T=3(T^2-P).                                      (15)
```

Therefore `delta|_U` is a unit. But every unit on `A2` is a scalar, while
`delta` cannot equal a scalar in `A`: its unique expression in the free
`R`-basis `(1,T,T^2)` has `T^2` coefficient `3`. This contradiction is
direct; it does not invoke THM-3862 circularly. Hence

```text
L squarefree  =>  no same-function-field planar Keller chart. (16)
```

## 6. Exact hostile controls

The forced-square distinction and both debt mechanisms occur in small
irreducible examples.

1. For `q=Pt^2`, Eisenstein at the prime `P` proves `F` irreducible, while

   ```text
   K=h^2(t^2-2h),                 L=t^2-2h.               (17)
   ```

   The adjusted curve is smooth, so the surface is normal even though
   `q1=t^2` has a double zero. Its only candidate at `h=0` is the isolated
   origin.

2. For `q=P^2+P`, viewing `-F` as a quadratic in `P` gives discriminant

   ```text
   (T+1)^2(4T+1),                                        (18)
   ```

   which is not a square in `k(T)`. Thus `F` is irreducible in
   `k[P,T]`, and Gauss applied to the primitive monic `T`-polynomial gives
   the required irreducibility over `k(P,t)`. Here

   ```text
   K=h^2(h-1)^2,                  L=(h-1)^2,              (19)
   ```

   and the entire line `P=1,T=-1` is singular. The forced `h=0` component
   itself remains regular.

3. For `q=P^2`, the same quadratic-in-`P` argument has discriminant
   `T^2(4T+9)`, again nonsquare, and the same Gauss step proves the domain
   gate. Meanwhile

   ```text
   K=h^3(h-2).                                            (20)
   ```

   The zero section `P=T=0` is singular, exactly as `(12)--(13)` predict.

4. For `q=t`, the surface is actually `A2` after solving
   `t=T^3-3PT`; it is smooth, and `K=t-2h^3` is squarefree. This is a
   positive normality control, not a Keller counterexample: the finite map
   still ramifies where `T^2=P`.

## 7. What is closed, and the exact remaining invoice

Within the irreducible arbitrary-`q` globally monogenic grammar `(1)`, the
normal locus is completely classified and excluded by `(11)` and `(16)`.
Only two conductor-debt loci escape that argument:

```text
(A) L has a repeated irreducible factor M != h;
(B) P^2 divides q, so the zero section E0 is singular.    (21)
```

Their normalizations are proper overorders of `A`; they need not remain
monogenic, and the global element `F_T` need not generate the normalized
different. No Keller chart is asserted to exist there. Subsequent THM-3962
closes every `t`-constant `q(P)` cylinder, including both debt types, and
THM-3963 closes the moving scalar family `q=c(t)P^2`. General moving repeated
factors and arbitrary `P^2q2(P,t)` remain open normalization problems.

This theorem does **not** prove JC(2). It proves that any counterexample in
the arbitrary polynomial-`q` depressed-cubic lane must exploit one of the
two nonnormal conductor debts in `(21)`, or leave this globally monogenic
function-field presentation altogether.

## Reproduction

```bash
python3 04-computation/jc2_arbitrary_q_hidden_repetition_normality_thm3961.py
python3 -O 04-computation/jc2_arbitrary_q_hidden_repetition_normality_thm3961.py
sha256sum 04-computation/jc2_arbitrary_q_hidden_repetition_normality_thm3961.py \
  05-knowledge/results/jc2_arbitrary_q_hidden_repetition_normality_thm3961.out
python3 agents/check_docs.py
```
