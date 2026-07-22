---
id: THM-2129
title: "Quartic Faber three-coefficient boundary classification and square-face carrier"
status: >
  PROVED. For every reduced quartic Faber degree n not divisible by four, the
  centered boundary value, first flux, and second flux are three consecutive
  coefficients of (1+4u+beta u^2+gamma u^3)^(n/4). Their common zero set is
  empty when n is odd and is the single point (beta,gamma)=(4,0) when n is
  twice odd. In depressed coordinates this point is
  z^4-2z^2+1=(z^2-1)^2. More generally, every nontrivial proper-power
  positive-weight face of a depressed quartic that contains z^4 is a square,
  with canonical lower sidecar (q,r-p^2/4). This classifies the balanced
  quartic boundary collision but does not close the quartic source-fiber
  stratum or JC(2).
source: codex-2026-07-22-JC2-quartic-three-coefficient-boundary
depends_on:
  - THM-2102
related:
  - THM-2084
  - THM-2118
  - THM-2127
  - THM-2132
  - HYP-3132
  - HYP-3138
  - HYP-3161
script: 04-computation/jc2_quartic_faber_three_coefficient_boundary_codex_20260722.py
output: 05-knowledge/results/jc2_quartic_faber_three_coefficient_boundary_codex_20260722.out
script_sha256: 6204284fc5e6e32dbf395b03a897b907ca524d66c1c186f41ce5ad190bf21cbb
output_sha256: 09b078a76512ca7d756c9e0f5321b7c684e21f3dd45e9ae4bef81a2e1204e136
hash_basis: repository blobs with LF line endings
---

# THM-2129 -- quartic boundary triple and the square carrier

## 1. The quartic tail has two independent fluxes

Work over a characteristic-zero differential field containing `C(x)`. Put

```text
P_0=z^4+p z^2+q z+r,                 alpha=n/4,        (1)
```

where `n>=1` and `4` does not divide `n`. Let

```text
P_0^alpha=E_n+c_1 z^-1+c_2 z^-2+c_3 z^-3+...,         (2)
```

with `E_n=Pol_z(P_0^alpha)`, and define

```text
Phi_n=4c_1,
Psi_n=4c_2,
R_n=4c_3+p c_1.                                      (3)
```

The condition `4 not|n` is exactly the reduced Faber condition: if `n=4d`,
then `E_n=P_0^d`, so a target shear removes that term from a reduced mate.

For the Hamiltonian derivation

```text
L=P_(0,x) partial_z-P_(0,z) partial_x,                 (4)
```

one has the exact polynomial identity

```text
L(E_n)=(z^2+p/4) Phi_n' + z Psi_n' + R_n'.            (5)
```

Indeed `L(P_0^alpha)=0`. If `S=P_0^alpha-E_n`, then the polynomial part of
`-L(S)` receives contributions only from `c_1,c_2,c_3` and is

```text
4c_1' z^2+4c_2' z+4c_3'+2p c_1'+p'c_1,               (6)
```

which is (5). Thus a constant-coefficient Faber combination satisfying a
nonzero constant Keller equation has two constant flux combinations,
`Phi'=Psi'=0`, while `R'` carries the Jacobian. The middle flux is a genuine
sidecar; a scalar boundary/flux pair would lose it.

## 2. Translation makes the boundary triple consecutive

Impose the centered boundary

```text
r=-1-p-q,                                             (7)
```

and set

```text
beta=6+p,                    gamma=4+2p+q,
T(u)=1+4u+beta u^2+gamma u^3,
A_k=[u^k]T(u)^(n/4).                                  (8)
```

With `t=z-1`, equation (7) gives

```text
P_0(t+1)=t^4+4t^3+beta t^2+gamma t
          =t^4 T(t^-1).                               (9)
```

Translation preserves the polynomial part at infinity. Its constant term is
therefore

```text
E_n(1)=A_n.                                           (10)
```

It also preserves the residue, so `c_1=A_(n+1)`. Expanding
`z^-1=(t+1)^-1=t^-1-t^-2+...` and
`z^-2=t^-2+...` shows

```text
A_(n+2)=c_2-c_1.                                      (11)
```

Consequently

```text
E_n(1)=Phi_n=Psi_n=0
  iff A_n=A_(n+1)=A_(n+2)=0.                          (12)
```

This is the quartic analogue of THM-2118's adjacent-coefficient carrier, but
the extra depressed coefficient requires three consecutive entries.

## 3. The coefficient recurrence

Write `F(u)=T(u)^(n/4)=sum_(k>=0) A_k u^k`. The identity

```text
T F'=(n/4)T'F                                         (13)
```

gives, with `A_j=0` for `j<0`,

```text
(k+1)A_(k+1)
 +4(k-n/4)A_k
 +beta(k-1-n/2)A_(k-1)
 +gamma(k-2-3n/4)A_(k-2)=0.                          (14)
```

All divisions used below are by displayed nonzero rational numbers, so there
is no analytic continuation or generic-parameter assumption.

## 4. A cubic term makes three zeros impossible

Assume the three zeros in (12). More generally, if

```text
A_(j+1)=A_(j+2)=A_(j+3)=0,                            (15)
```

then (14) at `k=j+2` gives

```text
gamma(j-3n/4)A_j=0.                                   (16)
```

Because `4` does not divide `n`, the number `3n/4` is not an integer. If
`gamma!=0`, start with `j=n-1` and descend (16) to `A_0=0`, contradicting
`A_0=1`. Hence every common boundary zero satisfies

```text
gamma=0.                                              (17)
```

## 5. Odd degrees are empty

Now `T=1+4u+beta u^2`. If `beta=0`, then

```text
A_n=4^n binom(n/4,n)!=0,                              (18)
```

since `n/4` is not an integer in `{0,...,n-1}`. Thus `beta!=0` at a putative
common zero. If two consecutive coefficients vanish, (14) becomes

```text
A_(j+1)=A_(j+2)=0  implies  beta(j-n/2)A_j=0.          (19)
```

For odd `n`, `n/2` is not an integer. Starting from `A_n=A_(n+1)=0`, equation
(19) again descends to the contradiction `A_0=0`. Therefore the common zero
set in (12) is empty for every odd reduced degree.

## 6. Twice-odd degrees have one square point

Let `n=2s` with `s` odd. For `s>1`, equation (19) descends from
`A_n=A_(n+1)=0` until it gives `A_(s+1)=0`; recall from (18) that `beta!=0`.
For `s=1`, this is already the assumed equality `A_2=0`. In either case
`A_(s+1)=A_(s+2)=0`. Solving the
order-two recurrence (14) upward now makes every `A_j` with `j>s` vanish.
Thus

```text
F(u)=T(u)^(s/2)                                       (20)
```

is a polynomial of degree at most `s`. Squaring in `C[u]` gives

```text
F(u)^2=T(u)^s.                                        (21)
```

Since `s` is odd and the right side is a square in the UFD `C[u]`, `T` itself
is a square. Its constant and linear coefficients force

```text
T(u)=(1+2u)^2,
(beta,gamma)=(4,0).                                   (22)
```

Conversely, (22) gives `F=(1+2u)^s`, whose coefficients above degree `s`
vanish, so (12) holds. Returning through (8) and (7), the unique point is

```text
(p,q,r)=(-2,0,1),
P_0=z^4-2z^2+1=(z^2-1)^2.                             (23)
```

This proves the all-degree classification.

## 7. Every nontrivial proper-power quartic face is a square

The exceptional point is structural, not an accidental resultant zero. Let

```text
P=z^4+p(x)z^2+q(x)z+r(x)                              (24)
```

and take a positive weight whose leading face contains `z^4` and at least one
other term. If that face is `lambda H^m`, `m>=2`, comparison of `z`-degrees
gives `m in {2,4}`. For `m=4`, write the degree-one polynomial `H=az+b`.
The leading `z^4` coefficient makes `a` constant and nonzero, while the absent
`z^3` coefficient gives `b=0`, leaving only `z^4`, a contradiction. For
`m=2`, write `H=az^2+bz+c`. Again `a` is constant and the absent `z^3` term
gives `b=0`. Therefore the face is a square, its `qz` term is absent, and

```text
r_face=p_face^2/4.                                    (25)
```

Globally complete the square:

```text
H_0=z^2+p/2,                delta=r-p^2/4,
P=H_0^2+qz+delta.                                     (26)
```

The exact square-completion remainder is the ordered pair `(q,delta)`. Since
its `z`-degree is below two,

```text
H_0 divides qz+delta  iff  q=delta=0.                  (27)
```

The latter would make `P=H_0^2`, which cannot be a Keller component because
`{P,Q}=2H_0{H_0,Q}` cannot be a nonzero constant. Therefore the first nonzero
weighted face of the remainder is the nonliftable approximate-root quotient
anticipated by THM-2102. Terminal, resonant, and all-order handling still
requires THM-2127's hypotheses or a new argument.

The claim is scoped to a positive-weight face containing `z^4`; a coefficient-
only face omitting `z^4` is not classified here.

An unrelated LRC reflection-fold packet suggested the correct carrier.
HYP-3132 isolates an even biquadratic resolvent, while HYP-3138/3161 warn that
the discarded odd coordinate must be resurrected. Under the typed map

```text
even fold -> H_0^2,             odd leakage -> qz,
remaining even offset -> delta,                                  (28)
```

the preserved lesson is exactly the need for the ordered `(q,delta)` sidecar.
No LRC inequality or tournament invariant is transferred to the Jacobian
problem; only the loss ledger survives.

## 8. Frontier and exact referee

The theorem eliminates every odd balanced quartic centering collision and
identifies the twice-odd collision with the unique proper-square boundary.
It does not prove that all centering poles are balanced, remove the remaining
quartic coefficient poles, or terminate arbitrary square-face resonances.
Therefore quartic source fibers, JC(2), and DC(2) remain open.

For a nonmonic quartic `A(x)y^4+...` with reduced mate leading coefficient
`q_n(x)`, the top Keller equation is

```text
n A' q_n-4A q_n'=0,             hence q_n^4=cA^n.      (29)
```

If `n` is odd, UFD valuations make `A` a fourth power. If `n` is twice odd,
they force only `A=V^2`; monic depression may then require the quadratic
extension adjoining `U` with `U^2=V`. The local theorem remains valid after
that extension, but a global quartic closure must retain and descend this deck
parity. No polynomial-coordinate conclusion is hidden here.

The companion works over exact rational polynomial rings. It independently
checks (5), compares (14) with the direct multinomial formula, verifies the
translation triangle (11), and computes Groebner certificates through `n=14`:
odd ideals are the unit ideal, while the
twice-odd ideals have radical `(beta-4,gamma)`. It uses explicit exceptions,
so normal and optimized runs execute the same checks. The recurrence proof,
not the finite sweep, supplies every all-degree quantifier.

Tournament Analysis is not faithful here. There is no intrinsic binary
relation or tie Hamiltonian path: the ordered coefficient triple, recurrence
index, parity, and square-defect sidecar are all essential. QED.
