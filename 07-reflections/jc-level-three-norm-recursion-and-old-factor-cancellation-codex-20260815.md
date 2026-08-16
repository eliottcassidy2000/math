# Level-three Keller norm recursion: what level two forces, and what it cannot

**Status: CONDITIONAL ALGEBRAIC SYNTHESIS + FINITE-EXACT THREE-SLICE
EVIDENCE.**  This note computes no global next image divisor and proves no
higher-iterate component theorem.  Its exact companion works only on the
three declared one-parameter slices.

## Inheritance pass

The closest proved mechanism is THM-2582's discriminant-of-a-norm theorem.
The canonical hostile is the old prediction that the grade-two odd divisor
was `LH`: THM-2582 proves that `L` cancels and only `H` survives.  The
corrected near miss is therefore to assume that the same cancellation happens
again merely because it happened once.  The least-used relevant sidecar is
not the set tower but the **valuation of the function-field norm along the
old nonproperness divisor**.

Write `N` for the degree-three function-field norm induced by the fixed
sporadic map `F`, and write `N^j` for the transitive norm from a `j`-step
generic inverse fibre.  Let `Delta_r` be the discriminant square class of the
degree-`3^r` `x`-eliminant for `F^r`, in the normalization-insensitive group
`K^*/K^{*2}`, where `K=Q(a,b,c)`.

## 1. The all-level formal recursion

Suppose at levels through `r+1` that:

1. the generic inverse algebra is finite and separable of degree three;
2. the level-`r` eliminant over each first inverse sheet is separable;
3. those three degree-`3^r` blocks are pairwise coprime;
4. the product of the blocks is the intended level-`r+1` eliminant, up to a
   nonzero element of `K`.

THM-2582 then applies with odd block degree `3^r` and gives

```text
[Delta_(r+1)] = [N(Delta_r)] [Delta_1]^(3^r)
                = [N(Delta_r)] [Delta_1].                 (1)
```

Induction yields the exact conditional recursion

```text
[Delta_r] = product_(j=0)^(r-1) [N^j(Delta_1)]
          = [(-1)^r product_(j=0)^(r-1) N^j(L)],          (2)
```

because THM-2473 gives `[Delta_1]=[-L]` and every norm degree `3^j` is odd.
This is a genuine general recursion for the square class of the recursively
defined norm-product eliminants.  It does **not** identify the irreducible
divisors of the iterated norms.

At level two, THM-2582 proves the exact identity

```text
N(L)=H/(64L),                                             (3)
```

so `(2)` recovers `[Delta_2]=[H]`.  At level three, equations `(1)` and `(3)`
already force two equivalent formulas without constructing the next image
divisor:

```text
[Delta_3]=[-L N(H)]                                      (4)
         =[-H N^2(L)].                                   (5)
```

Indeed `H=64 L N(L)`, and after applying the cubic norm the factor `64^3`
is a square.  Equation `(5)` is useful bookkeeping, but it is not new
geometric information: `N^2(L)` is exactly the unknown that `(4)` hides in
`N(H)`.

## 2. The exact valuation boundary

For every irreducible polynomial `R` in `Q[a,b,c]`, denominator clearing does
not alter discriminant square class because a degree-27 polynomial scales its
discriminant by the even power `52`.  Equations `(4)` and `(5)` therefore give

```text
v_R(Delta_3) = 1_(R=L) + v_R(N(H))        (mod 2),        (6)
             = 1_(R=H) + v_R(N^2(L))      (mod 2).        (7)
```

Thus the level-two identity determines the level-three class as an element of
`K^*/K^{*2}`, but it does **not** determine its prime-divisor parity until the
valuations of `N(H)` are known.

To state the geometric conditional precisely, suppose additionally that:

- `closure(F(V(H)))` has dimension two and prime equation `J`;
- the norm divisor has no numerator component other than `J` and no pole
  component other than the old boundary `L` after honest saturation;
- `J` is distinct from `L` and `H`; and
- for some `u in Q^*` and integers `e,d>0`,

```text
N(H)=u J^e/L^d times a square in K^*.                     (8)
```

Then

```text
[Delta_3]=[-u L^(1-d) J^e],                               (9)
```

and the valuation ledger is

```text
v_J(Delta_3) = e,       v_L(Delta_3) = 1-d,
v_H(Delta_3) = 0,       v_R(Delta_3) = 0 otherwise        (mod 2).  (10)
```

The “newest factor only” conclusion requires **both** `e` and `d` odd, plus
control of the constant class `[-u]`.  Separability and pairwise coprimality
justify the norm-product discriminant formula; they do not prove `(8)`, odd
pole order, image distinctness, or image multiplicity.  Irreducibility of the
closure follows once `V(H)` is irreducible, but being a hypersurface and being
distinct from earlier closures remain separate obligations.

## 3. Finite-exact slice probe

The companion evaluates the transported 361-term polynomial `H` inside the
cubic algebra

```text
Q(A)[w]/(L(A,b_0,c_0)w^3+(4-3b_0c_0)w-2c_0)
```

using the exact inverse graph from THM-2576.  The determinant of
multiplication by `H(q(w))` is `N(H)`; no degree-27 eliminant or discriminant
is expanded.  On every slice the same quotient algebra first rederives the
proved control `N(L)=H/(64L)`, and the multiplication determinant for `N(H)`
is checked against the independent resultant of the reduced element with the
monic cubic.

On all three slices it finds

```text
(b,c)=(1,2):  N(H)=K_(1,2)/(2^21 L^7),
(b,c)=(3,1):  N(H)=K_(3,1)/(2^35 L^7),
(b,c)=(1,3):  N(H)=K_(1,3)/(2^35 L^7).                    (11)
```

Each `K_(b,c)` is primitive and irreducible of degree `86` over `Q`, and is
coprime to both specialized old factors `L` and `H`.  Each quadratic `L` is
separable (`disc=-32,-500,-500`).  Therefore on every tested slice,

```text
[Delta_3]=[-L N(H)]=[-2K_(b,c)],                          (12)
```

and the old factor has exponent `1-7=-6`, which is even.  There is no hidden
`H` factor and no hidden cancellation between the new numerator and `L`.

This is strong evidence for a global pole order seven and a new degree-86
image slice, but it is not a proof: specialization can merge components,
change content, or miss a component supported on the frozen parameter locus.
It also does not prove generic separability or pairwise coprimality of the
three degree-nine blocks needed to turn `(4)` into a statement about the
actual degree-27 eliminant.

## 4. Cheapest decisive global computation

The next computation should continue to avoid the degree-27 discriminant.
Work directly in `Q(a,b,c)[w]/(E)` and compute the `3 by 3` multiplication
determinant of `H(q(w))`.  Then:

1. prove that `L^7 N(H)` is polynomial up to an explicit rational constant;
2. primitively normalize and factor/saturate its numerator `J`;
3. prove `gcd(J,LH)=1`, exhibit a rank-two point on `F(V(H))`, and separate
   `V(J)` from both earlier image components by exact points;
4. determine the generic divisor multiplicity `e`; and
5. independently produce one squarefree degree-27 specialization to certify
   generic block separability and pairwise coprimality.

Steps 1--4 decide the missing valuation in `(6)` and the component question;
step 5 licenses the discriminant recursion.  Computing the full degree-27
discriminant first would spend vastly more algebra while obscuring the one
load-bearing datum: the odd pole order of `N(H)` along `L`.

## Scope

The proved content here is the abstract conditional recursion `(1)--(2)` and
the three finite exact norm calculations `(11)`.  The global factor `J`, the
global identity suggested by the slices, generic degree-27 separability,
higher image distinctness, exact positive discriminant multiplicities, and
all Jacobian-conjecture classification consequences remain open.

Reproduce with

```bash
python3 04-computation/keller_level_three_norm_slice_probe_20260815.py
python3 -O 04-computation/keller_level_three_norm_slice_probe_20260815.py
```
