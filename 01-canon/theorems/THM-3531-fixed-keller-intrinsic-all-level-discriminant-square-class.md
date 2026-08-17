---
id: THM-3531
title: "Fixed Keller intrinsic all-level discriminant square class"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the degree
  3^n function-field extension induced by the nth iterate of the fixed
  sporadic Keller map, the basis-independent trace-discriminant square class
  is [(-1)^n P_(n-1)] in the raw cleared-norm normalization.  Hence the
  newest prime is the sole odd horizontal discriminant divisor at every
  level.  This theorem alone is not an all-level primitivity theorem;
  THM-3535 subsequently supplies that separate monodromy gate.
source: codex/intrinsic-discriminant-tower/2026-08-16
audit: >
  The independent audit checked tower orientation, trace-discriminant
  transitivity, raw scalar classes, odd clearing exponents, fixed-level
  generic primitive forms, and the distinction between intrinsic and
  coordinate discriminants.  It retained the mandatory constant class [2]
  and supplied even-block, even-clearing, rescaling, and nonprimitive-
  coordinate hostiles.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
  - THM-3528-fixed-keller-all-level-cleared-norm-polynomiality-and-finite-sheet-defect
  - THM-3530-fixed-keller-all-level-image-prime-and-component-tower
related:
  - MISTAKE-413
  - THM-3508-level-two-sporadic-keller-three-coordinate-primitive-discriminant-square-class
  - THM-3519-level-three-sporadic-keller-three-coordinate-primitivity-and-common-discriminant-class
  - THM-3525-level-five-degree243-separability-and-discriminant-square-class
  - THM-3526-level-six-degree729-separability-and-discriminant-square-class
  - THM-3533-fixed-keller-newest-prime-reduced-different-and-index-square
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_intrinsic_discriminant_raw_tower_audit_20260816.py
output: 05-knowledge/results/keller_intrinsic_discriminant_raw_tower_audit_20260816.out
script_sha256: 54a9d6977275feb97759b33be8ff90a677cc86222f44548f9c5f0952ec78e553
output_sha256: 16a1950ff176f08c7983f3b20429a8dede9ea1802d9e166081fd1dffb9ea63e0
semantic_sha256: 9ccc7adc713ebf30e69af6456e212b94f4a43da4cec6215f730c58cc00ecfe46
hash_basis: LF-normalized bytes
---

# THM-3531 -- the intrinsic discriminant follows the newest raw prime

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Retain the fixed sporadic Keller map `F`, and let

```text
K_0=Q(a,b,c),
K_n=Q(x_n,y_n,z_n),             (a,b,c)=F^n(x_n,y_n,z_n). (1)
```

Then `K_n/K_0` is a finite separable function-field extension of degree
`3^n`.  Let

```text
d_n in K_0^*/K_0^*2                                  (2)
```

be its intrinsic discriminant square class: the determinant of the trace
pairing in any `K_0`-basis.  A basis change multiplies this determinant by a
square, so `(2)` is canonical.

Use THM-3528's exact raw tower

```text
P_0=L,
P_(j+1)=L^e_j N(P_j),
(e_(j+1),m_(j+1))=(7e_j-2m_j,3e_j-2m_j).              (3)
```

## 1. The theorem

For every `n>=1`,

```text
d_n=[(-1)^n P_(n-1)] in K_0^*/K_0^*2.                 (4)
```

By THM-3530, `P_(n-1)` is absolutely irreducible and the raw primes are
pairwise distinct.  Therefore the sole odd horizontal prime divisor of the
intrinsic discriminant class at level `n` is

```text
V(P_(n-1)).                                             (5)
```

Equation `(5)` is a parity statement.  It neither gives the exact positive
discriminant multiplicity nor identifies a discriminant ideal on the
nonfinite affine coordinate-ring map.

## 2. One step has class `[-L]`

THM-2473 gives the primitive one-step cubic

```text
E(w)=Lw^3+(4-3bc)w-2c                                 (6)
```

with

```text
disc_w E=-4(27ac^2-9bc+8)^2 L.                        (7)
```

Passing to the monic cubic changes its discriminant by `L^4`, a square.
Thus

```text
d_1=[-L]=[-P_0].                                       (8)
```

This starts `(4)` and retains the sign.  The affine map need not be finite:
the dominant generically finite map still induces the finite separable field
extension used in `(1)`--`(2)`.

## 3. Trace-discriminant transitivity has the odd-block orientation

Choose a generic inverse chain

```text
F(q_i)=q_(i-1),               K_i=Q(q_i).              (9)
```

For the tower

```text
K_0 subset K_1 subset K_(n+1),                         (10)
```

the extension `K_(n+1)/K_1` is the same `n`-level construction after base
change to the first inverse point.  Trace-discriminant transitivity gives

```text
d(K_(n+1)/K_0)
 =Norm_(K_1/K_0)(d(K_(n+1)/K_1))
  d(K_1/K_0)^[K_(n+1):K_1].                            (11)
```

The exponent in `(11)` is `3^n`, hence odd.  Evaluating the `n`-level class
at `q_1` and taking the one-step cubic norm therefore yields exactly

```text
d_(n+1)=[N(d_n)] [d_1].                                (12)
```

This is the field-theoretic form of THM-2582's odd-block discriminant law.
No coordinate polynomial or block-coprimality computation is needed for
`(12)`.

## 4. Raw norm normalization closes the induction

Suppose `(4)` holds at level `n`.  The cubic norm of the constant sign keeps
that sign because its degree is odd.  Equations `(3)`, `(8)`, and `(12)` give

```text
d_(n+1)
 =[N((-1)^nP_(n-1))] [-L]
 =[(-1)^(n+1) L N(P_(n-1))]
 =[(-1)^(n+1) P_n L^(1-e_(n-1))].                    (13)
```

The first packet coordinate is always odd:

```text
e_(j+1)=7e_j-2m_j=e_j mod 2,              e_0=1.      (14)
```

Hence `1-e_(n-1)` is even, its power of `L` disappears in square classes,
and `(13)` becomes

```text
d_(n+1)=[(-1)^(n+1)P_n].                               (15)
```

This proves `(4)` for all `n`.

## 5. The raw scalar class is load-bearing

The first six raw-to-named normalizations are

| rung | raw polynomial |
|---:|---|
| `P_0` | `L` |
| `P_1` | `H/2^6` |
| `P_2` | `J/2^53` |
| `P_3` | `G/2^159` |
| `P_4` | `R_5/2^477` |
| `P_5` | `R_6/2^1431` |

Here `53=3*6+35` uses the exact identity
`N(H)=J/(2^35L^7)`; after `P_2`, each displayed denominator exponent cubes
under the cubic norm.  Substitution into `(4)` gives

```text
n=1: [-L],
n=2: [ H],
n=3: [-2J],
n=4: [ 2G],
n=5: [-2R_5],
n=6: [ 2R_6].                                         (16)
```

These are exactly the previously audited classes.  In particular, the
constant `[2]` at levels three through six is not disposable.  Dropping it
repeats MISTAKE-413.

If a named normalization is `Q_j=c_jP_j`, then `(4)` reads

```text
d_(j+1)=[(-1)^(j+1)c_j^(-1)Q_j].                      (17)
```

An arbitrary rational rescaling can therefore change the displayed square
class even though it leaves the zero divisor unchanged.

## 6. Generic primitive forms exist separately at every depth

For each fixed `n`, there is a Zariski-open set of coefficient triples
`(alpha,beta,gamma)` for which

```text
theta=alpha x_n+beta y_n+gamma z_n                    (18)
```

is primitive for `K_n/K_0`.  Indeed, a finite separable extension has only
finitely many intermediate fields.  Their intersections with the
three-dimensional span in `(18)` are proper linear subspaces: if one
contained all three source coordinates, it would contain `K_n`.  An infinite
base field is not the union of finitely many proper subspaces.

The minimal polynomial of a primitive `theta` is separable of degree `3^n`,
and its discriminant has class `(4)`.  Total primitivity also makes the
conjugate block polynomials separable and pairwise coprime in any chosen
tower decomposition.

This is only a fixed-level generic argument.  THM-3535 subsequently avoids
the invalid countable-intersection inference and proves, by full wreath
monodromy plus one-step non-descent, that every nonzero constant linear form
works simultaneously for every `n`.

## 7. Intrinsic is not literal-coordinate separability

THM-3508 supplies the sharp hostile already at depth two.  Its intermediate
coordinate `t` has a separable cubic minimal polynomial, but multiplication
by `t` on the rank-nine algebra has characteristic polynomial equal to the
cube of that cubic.  The degree-nine power matrix has rank three and the
naive degree-nine eliminant has discriminant zero.

Thus a nonzero intrinsic field discriminant does not make a chosen
nonprimitive coordinate eliminant squarefree or full degree.  The literal
`x/y/z` claims through the audited levels remain stronger, separate theorems.

The other sharp boundaries are:

1. if a clearing exponent were even, `L^(1-e)` in `(13)` would retain the old
   boundary class;
2. if the block degree were even, the outer class in `(11)` would disappear;
3. rescaling a raw rung by `2` changes the constant square class; and
4. the theorem gives parity, not exact positive multiplicity.

No arbitrary Keller-map law, exact discriminant multiplicity, `JC(2)`,
`DC(2)`, LRC, or general Jacobian-conjecture classification follows from
this theorem.  THM-3533 separately proves newest-prime normalization
multiplicity one, and THM-3535 supplies the simultaneous-coordinate result.

## Reproduction

```text
python -B 04-computation/keller_intrinsic_discriminant_raw_tower_audit_20260816.py
python -B -O 04-computation/keller_intrinsic_discriminant_raw_tower_audit_20260816.py
```

Normal and optimized transcripts match the stored output.  The companion
checks odd clearing through depth 31 plus its symbolic parity recursion, the
trace-discriminant sign ledger, all six raw scalar invoices, and every hostile
listed above.
