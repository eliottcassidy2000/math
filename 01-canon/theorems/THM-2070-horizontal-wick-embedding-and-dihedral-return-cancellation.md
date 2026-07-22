---
id: THM-2070
title: "Horizontal Wick embedding and dihedral return cancellation"
status: >
  PROVED. Every finite one-variable Laurent polynomial embeds exactly as a
  horizontal-support two-variable Gaussian polynomial: for B large enough,
  P_f(Z,W)=sum_q c_q Z^(B+q)W^B satisfies
  E[P_f^m]=(Bm)! CT(f^m) for every m. Thus the face Laurent polynomial in the
  GMC(2) proof is not a restricted subclass; replacing its DvdK input requires
  reproving the arbitrary one-variable existence theorem. Moreover
  f=u^2+u+u^(-1)-u^(-2) has support return set {m>=2} and gcd-one support gaps,
  but CT(f^m)=0 for every odd m and also for m=2, while CT(f^4)=-12. The
  mechanism is the coefficient-level involution f(-u^(-1))=-f(u). Hence the
  support numerical semigroup controls feasibility only, not complex-weighted
  noncancellation; S222's aperiodic/eventual-nonzero claim and S223's
  finite-sporadic-cancellation/completed-bypass claim are refuted. Their
  positive-coefficient and two-charge conclusions survive. The independently
  proved Galois orbit theorem in THM-2067 supplies the correct algebraic
  replacement for bare existence.
source: codex-2026-07-21-dvdk-referee
depends_on:
  - THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2
  - THM-1840-the-functional-agnostic-both-signs-single-character-nonvanishing
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2067-galois-orbit-product-closes-one-variable-dvdk
  - HYP-8890
  - HYP-8895
---

# THM-2070 -- horizontal Wick embedding and dihedral return cancellation

This theorem separates three objects that the proposed S222/S223 DvdK bypass
had blended:

1. the support return set, which records whether a balanced word exists;
2. the coefficient-weighted constant term, which may cancel on that return
   set; and
3. the exact lowest-face Laurent polynomial consumed by the GMC(2) proof.

The first does not determine the second, and the third is as general as an
arbitrary one-variable Laurent polynomial.

## 1. Every Laurent polynomial is a horizontal Wick face

Let

```text
f(u)=sum_(q in S) c_q u^q in C[u,u^(-1)]
```

have finite exact support `S`. Choose an integer

```text
B >= max(0,-min S)
```

and define

```text
P_f(Z,W)=sum_(q in S) c_q Z^(B+q) W^B in C[Z,W].          (1)
```

If `Z` is a circular complex Gaussian and `W=conj(Z)`, then for every
`m>=1`,

```text
E[P_f(Z,W)^m]=(Bm)! CT_u(f(u)^m).                         (2)
```

### Proof

An ordered word `(q_1,...,q_m)` in the expansion of (1) contributes the
monomial

```text
Z^(Bm+q_1+...+q_m) W^(Bm).
```

Its Gaussian expectation is zero unless `sum_i q_i=0`. In that balanced case
both exponents equal `Bm`, so its Wick weight is the same scalar `(Bm)!`.
After factoring out that scalar, the sum of the balanced coefficient words is
exactly `CT(f^m)`. This proves (2).

The charge of the monomial indexed by `q` is

```text
(B+q)-B=q,
```

and for the exposing slope `lambda=1` its tilted height is

```text
(B+q)-lambda*q=B.                                        (3)
```

Thus the entire exact support of `P_f` is one equality face, and its face
Laurent polynomial is literally `f`. There is no coefficient, support, or
charge restriction hidden in the GMC(2) lowest-face construction.

### Consequence for the dependency

The implication

```text
all Gaussian moments of every horizontal P_f vanish
  => the charges of P_f have one strict sign
```

already implies the one-variable statement

```text
CT(f^m)=0 for every m>=1
  => supp(f) is strictly one-sided.                       (4)
```

Therefore a proof of NC2 cannot bypass (4) by claiming that its face Laurent
polynomials form an easier subclass. It must either cite or reprove the bare
one-variable theorem. The Galois orbit-product proof in
`THM-2067-galois-orbit-product-closes-one-variable-dvdk.md` now does reprove
exactly (4); it replaces the DvdK citation without asserting a stronger
eventual-nonvanishing theorem.

## 2. A cofinite return set with infinitely many cancellations

Consider the real signed Laurent polynomial

```text
f(u)=u^2+u+u^(-1)-u^(-2),
S={-2,-1,1,2}.                                             (5)
```

Define the support return set

```text
R(S)={m>=1 : there exist q_1,...,q_m in S with
                   q_1+...+q_m=0}.
```

Then

```text
R(S)={m>=2}.                                               (6)
```

Indeed, every even `m=2k` is represented by `k` copies of the pair `(1,-1)`,
and every odd `m=3+2k` is represented by `(-2,1,1)` followed by `k` copies of
`(1,-1)`. There is no return of length one because `0` is not in `S`.
The support gaps have gcd one, so this is an aperiodic/cofinite support in the
precise sense used by S222/S223.

Nevertheless

```text
CT(f^m)=0 for every odd m.                                (7)
```

The reason is the coefficient-level Laurent involution

```text
tau(u)=-u^(-1),             f(tau(u))=-f(u).              (8)
```

Constant term is invariant under `u -> -u^(-1)`. Hence

```text
CT(f^m)=CT(f(tau(u))^m)=(-1)^m CT(f^m),
```

which proves (7) in characteristic zero.

The first even return cancels too:

```text
CT(f^2)=2(c_1 c_(-1)+c_2 c_(-2))=2(1-1)=0.               (9)
```

But this is not a null Laurent polynomial. Squaring (5) gives coefficients

```text
exponent:  -4 -3 -2 -1  0  1  2  3  4
f^2:        1 -2  1 -2  0  2  1  2  1,
```

so

```text
CT(f^4)
 =2[(2)(-2)+(1)(1)+(2)(-2)+(1)(1)]
 =-12 != 0.                                               (10)
```

Thus `R(S)` has conductor two, but its coefficient-weighted constant terms
vanish at every odd return length, an infinite set, and also at length two.
The return semigroup does not even predict the first nonzero power in this
four-term example.

### The same witness is an actual GMC(2) face

Taking `B=2` in (1) gives

```text
P_f(Z,W)=W^2(Z^4+Z^3+Z-1).                               (11)
```

All four monomials lie on the horizontal face `b=2`, and

```text
E[P_f^m]=(2m)! CT(f^m).                                  (12)
```

Consequently the first, second, and every odd Gaussian moment in this family
vanish, while

```text
E[P_f^4]=-12*8! != 0.
```

The hostile control is therefore internal to the exact face setting of
THM-2022, not an irrelevant Laurent-polynomial pathology.

## 3. The general missing sidecar is dihedral coefficient phase

The parity mechanism is not isolated. For `alpha in C^*` and
`epsilon in {1,-1}`, suppose

```text
f(alpha/u)=epsilon f(u).                                  (13)
```

Then constant-term invariance gives

```text
CT(f^m)=epsilon^m CT(f^m).                                (14)
```

For `epsilon=-1`, every odd constant term vanishes. In coefficients, (13)
means

```text
c_(-q)=epsilon alpha^q c_q.                               (15)
```

On a balanced return word, the factors `alpha^q` multiply to one; reversal
therefore changes its coefficient weight by `epsilon^m`. Odd return words
pair with equal magnitude and opposite sign even though their existence is
fully visible in `R(S)`.

This identifies the lost coordinate in the semigroup quotient:

```text
weighted Laurent polynomial  ->  support return semigroup
preserves: mass and existence of balanced words
destroys: coefficient phase on reversal/orbit classes
needed sidecar: the coefficient character (15), or the full orbit sum.
```

Minimal return words and numerical-semigroup conductors are therefore useful
selectors, but they cannot be noncancellation certificates without this
sidecar.

## 4. Audit of S222 and S223

### Refuted claims

The witness (5)--(10) refutes all of the following general mixed-coefficient
claims:

- aperiodic support implies a unique uncancelled dominant saddle;
- `CT(f^m)` is nonzero for all sufficiently large `m`;
- mixed signs remove only finitely many, or merely sporadic, elements from the
  support return set;
- the numerical semigroup or its Frobenius number completes the DvdK bypass or
  supplies the open effective first-nonzero bound.

There is an even earlier obstruction to S222's displayed positive-radius
equation for arbitrary complex coefficients: for

```text
g(u)=u+i u^(-1),
```

the proposed equation `sum_k k c_k r^k=0` becomes `r-i/r=0`, which has no
positive real solution. A complex critical point may exist, but it does not
restore the positive-coefficient triangle inequality or rule out several
equal-modulus contributions. Those are precisely the steps that need proof
in a genuine complex saddle argument.

### Strongest surviving statements

The following parts remain correct.

1. `R(S)` is an additive submonoid of the positive integers. When its gcd is
   one it is cofinite; otherwise it is periodic after dividing by its gcd.
2. For positive real coefficients,
   `CT(f^m)>0` exactly when `m in R(S)`. The support problem and the weighted
   problem coincide because cancellation is impossible.
3. For a two-charge polynomial

   ```text
   a u^p+b u^(-q),       a b!=0,
   ```

   the balanced channel is unique at every return length, and the exact
   THM-1840 binomial formula proves nonvanishing when
   `(p+q)/gcd(p,q)` divides `m`.
4. A saddle analysis may still prove asymptotics under additional positivity,
   dominance, or nonresonance hypotheses. It is not a coefficient-uniform
   replacement theorem as stated in S222.

The full arbitrary-complex bare existence statement is now supplied by the
Galois orbit norm in the Galois THM-2067 artifact. The support semigroup is a
useful feasibility skeleton; the orbit product, not the conductor, controls
the cancellations required by GMC(2).

## 5. A separate scope repair for THM-1840

The core two-charge constant-term formula in
`THM-1840-the-functional-agnostic-both-signs-single-character-nonvanishing.md`
is correct. Its stated definition of a `charge-graded` functional, however,
is too weak for the displayed identity

```text
F(Lambda^m)=F(1) CT(Lambda^m).
```

Merely saying that `F(u^k)` depends on the charge `k` and that `F(1)!=0` does
not make `F` annihilate nonzero charges. For example, let

```text
Lambda=u+u^(-1),
F(1)=1,  F(u^2)=F(u^(-2))=-1,
F(u^k)=0 otherwise.
```

Then `F` is diagonal by charge in the stated sense, but

```text
F(Lambda^2)=-1+2-1=0,
F(1) CT(Lambda^2)=2.
```

The exact repair is to require the charge-selection rule

```text
F(u^k)=0 for every k!=0,        F(1)!=0.                  (16)
```

Under (16), `F(h)=F(1)CT(h)` for every Laurent polynomial `h`, so the
functional corollary is immediate. This correction does not affect the true
two-charge constant-term result used as a positive control above.

## 6. Scope

- **PROVED:** the horizontal Wick embedding (1)--(3), the exact return set
  (6), the infinite dihedral cancellations (7), and the nonzero witness (10).
- **REFUTED:** the mixed-complex eventual-nonzero and semigroup-completeness
  claims in S222/HYP-8890 and S223/HYP-8895.
- **SURVIVES:** positive coefficients, unique balanced channels, and the
  two-charge THM-1840 formula.
- **PROVED ELSEWHERE:** the Galois THM-2067 gives a self-contained algebraic
  proof of the weaker existence theorem actually needed by GMC(2).
- **OPEN:** a uniform effective bound such as the Sturmfels/ESV width bound.
  Neither a support conductor nor the Galois orbit proof supplies it.

No numerical experiment is used in the proof; every displayed identity is an
exact finite expansion.
