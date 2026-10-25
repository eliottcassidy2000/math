# The unordered Frobenius repair has an effective torus-count model in every dimension

Status: **PROVED**, with independent finite-field controls. Session59, 2026-10-25. This completes the coefficient-positivity gate left open in `RESULT.md` §§6/9. It does **not** prove that the actual unordered configuration quotient admits a torification.

## 1. A positive surprise after the ordered obstruction

The ordered primitive carrier has negative first nonzero `(q-1)` coefficient for every `k>=3`. Restoring the full geometric unordered count, including nonrational Frobenius orbits, reverses that obstruction: its entire `(q-1)` polynomial has **nonnegative integer coefficients** in every dimension. Consequently it has an explicit count-equivalent finite disjoint union of split tori, hence a standard monoid/F1 count model. The distinction between *count-equivalent model* and *descent of the original quotient* is essential.

The parent independently derived the primitive Frobenius formula and proposed the positivity test after checking `k<=8`; the factorization below proves the all-`k` result. No literature-priority claim is made.

## 2. Exact primitive Frobenius formula

Keep `C_k,D_k,U_k` from `RESULT.md`, and let `V_k(q)` count Frobenius-stable unordered triples with **every coordinate active**. Put `V_0=V_1=0`. Let

```text
H_k(q)=1+q+...+q^(k-1).
```

For `k>=1`, the type-`12` condition now forces the conjugate pair's common support to be all coordinates; the rational third point is arbitrary. A type-`3` orbit must likewise have full support. Thus the three cycle types give

```text
6V_k(q) = D_k(q)
          +3[(q^2-1)^(k-1)-(q-1)^(k-1)] H_k(q)
          +2[(q^3-1)^(k-1)-(q-1)^(k-1)].       (V)
```

At `k=1` both brackets vanish and `D_1=0`. This is an exact count of geometric unordered configurations; it is **not** the arithmetic primitive count `M_k`, whose split part is `D_k(3)/6`.

## 3. Positive integral factorization

For `k>=2` and `t=q-1`, substitute the primitive formula (P) from `RESULT.md`, `H_k(1+t)=((1+t)^k-1)/t`, and `q^3-1=t(t^2+3t+3)` into (V). Collecting the two diagonal-subtraction terms gives

```text
6 V_k(1+t) = t^(k-2) B_k(t),

B_k(t) = t^(k-1)[(t+3)^k-3]
         +3[(t+2)^(k-1)-2][(t+1)^k-1]
         +2t(t^2+3t+3)^(k-1).                 (B)
```

Every coefficient of every summand of `B_k` is nonnegative:

- `(t+3)^k-3` has positive constant term for `k>=2` and nonnegative higher coefficients;
- `(t+2)^(k-1)-2` has constant `2^(k-1)-2>=0` and nonnegative higher coefficients;
- `(t+1)^k-1` and `t^2+3t+3` have nonnegative coefficients.

It remains to show coefficientwise divisibility by six—not merely integer values at field sizes. In `F_3[t]`,

```text
B_k(t) = t^(2k-1)+2t^(2k-1)=0.
```

In `F_2[t]`, the third summand vanishes and the first two are equal:

```text
B_k(t) = t^(k-1)[(t+1)^k-1]
         +t^(k-1)[(t+1)^k-1]=0.
```

Therefore every coefficient is divisible by both two and three, hence by six. Formula (B) proves

```text
V_k(1+t) belongs to Z_(>=0)[t]   for every k>=0.
```

This is stronger than finite-field integer-valuedness, which alone would not imply integer polynomial coefficients.

## 4. Explicit effective carrier and its limits

Write `V_k(1+t)=sum_r v_(k,r)t^r`. Define

```text
E_k^prim = disjoint_union_r [v_(k,r) copies of G_m^r over Z].
```

This is a finite disjoint union, explicitly determined by (B). It satisfies

```text
#E_k^prim(F_q)=V_k(q)
```

for every finite prime power `q`, and so has the same extension-field counts and Hasse–Weil zeta as the primitive unordered configuration model whenever that geometric quotient is used. Each torus has its standard monoid-scheme model; their finite disjoint union is an honest effective F1/monoid *count carrier*.

No isomorphism, stratification map, common moduli functor, Grothendieck-ring identity, or arithmetic gap bijection to this torus union is asserted. In particular a counting model is not a torification of the actual quotient. The torus decomposition is canonical only after choosing the coefficient recipe and labels for the repeated copies, not canonically attached geometric charts on the quotient.

The active-coordinate union is Frobenius invariant. Its deletion/reinsertion partition therefore gives

```text
U_k(q)=sum_(j=0)^k binom(k,j)V_j(q).             (UA)
```

Thus `U_k(1+t)` also lies in `Z_(>=0)[t]` for every `k`, and the same construction provides an effective total torus-count carrier `E_k`. This repairs the numerical positivity obstruction without identifying the models.

## 5. First surviving rank and controls

For `k>=3`, the first summand in (B) begins at degree at least two. The coefficient of `t` in the other two summands is

```text
3k(2^(k-1)-2)+2*3^(k-1).
```

Consequently the lowest torus rank of `E_k^prim` is `k-1`, with multiplicity

```text
v_(k,k-1)=k*2^(k-2)-k+3^(k-2),  k>=3.
```

The `k=2` boundary has an additional linear contribution from the first summand: `V_2(1+t)=t^3+3t^2+2t`, so its minimum-rank multiplicity is **2**, not the value1 obtained by incorrectly applying the `k>=3` formula. The top coefficient is1 at rank `3k-3`.

The minimum-rank multiplicities for `k=2,...,8` are

```text
2, 6, 21, 62, 171, 460, 1233.
```

Examples:

```text
V_3(1+t)=t^6+7t^5+18t^4+19t^3+6t^2,
V_4(1+t)=t^9+10t^8+42t^7+92t^6+109t^5+69t^4+21t^3.
```

For total unordered configurations, (UA) gives

```text
[t] U_k(1+t)=k(k-1),
[t^2] U_k(1+t)=3binom(k,2)+6binom(k,3).
```

Compare the ordered coefficient `3binom(k,2)-7binom(k,3)`: the cycle correction changes its sign pattern, while the interpretations of the two counts remain different. At `q=3,k=3`, the primitive split arithmetic part is138, type12 contributes390, and type3 contributes224, giving `V_3(3)=752`.

## 6. Verification and connection ledger

The shared `verify.py` now checks (V) against (B) by two polynomial computations through `k=20`, checks divisibility/nonnegativity coefficientwise, checks (UA)'s two lowest coefficients, and independently enumerates primitive Frobenius orbit types over explicit extensions `F_(q^2),F_(q^3)` for prime `q=2,3,5`, `k=1..3`. Normal and optimized replays use explicit exception gates. Read the matching current hashes before reproducing.

Connection: ordered primitive filter → geometric unordered filter → effective torus-count model. The first arrow forgets order **and restores Frobenius descent types**; the second sends a counting polynomial to copies of standard tori. Preserved by the second arrow: every finite-field point count, extension zeta, and polynomial absolute-zeta transform. Lost: incidence, moduli, topology, canonical points, and any numerical gap meaning. Decisive hostile to reversing it: positivity supplies no map back to the quotient. The strongest survivor is the explicit count-equivalent monoid model.

The new question is now sharply geometric: does the actual quotient admit a torification, or can an invariant invisible to point counts obstruct it? This note does not settle that question. No new external theorem is imported beyond the torus-count and zeta definitions already pinned in `RESULT.md`.
