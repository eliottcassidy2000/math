---
id: THM-2067
title: "A Galois orbit-product proof of one-variable constant-term nonvanishing"
status: >
  PROVED. Let f in C[z,z^-1] have a nonzero positive exponent and a
  nonzero negative exponent. Then CT(f^m) is nonzero for some m>=1. More
  generally, if CT(f^m)=0 for every m>=1, the nonzero exponents of f have one
  strict sign. THM-1550 turns total vanishing into Pi(t)=c t for the product
  of the M small roots of X^M=tR(X). The polynomial X^M-tR(X) is irreducible
  over C(t), so its Galois group is transitive. Multiplying Pi=c t over the
  orbit of that M-subset gives (ct)^r=C^eta, while Vieta makes C the nonzero
  constant (-1)^d r_0/r_d. The t-adic valuations disagree. This is a
  project-internal proof of exactly the bare existence input used by THM-2022;
  it does not prove DvdK's stronger critical-value/limsup theorem. THM-2087
  later refines the same root object to an unconditional effective, generally
  non-sharp first-return bound.
source: codex-2026-07-21-dvdk-galois-orbit
depends_on:
  - THM-1550-an-exact-criterion-for-the-toral-nullcone
related:
  - THM-1605-tnc-proved-monodromy-transitivity
  - THM-1630-tnc-is-duistermaat-van-der-kallen-theorem-2
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2058-primitive-phase-packets-and-deck-fan-intervals
  - THM-2087-effective-compound-root-bound-for-one-variable-constant-terms
---

# THM-2067 -- Galois orbit-products close the one-variable input

The one-variable theorem needed by the NC2/GMC(2) proof is weaker than the
full Duistermaat--van der Kallen theorem. It asks only:

> If a complex Laurent polynomial has exponents of both signs, must one of
> its positive powers have nonzero constant term?

The answer follows from THM-1550 and a finite Galois orbit count. This is the
algebraic version of THM-1605's monodromy orbit-product proof; it removes
analytic continuation from the closing step and isolates a short norm lemma.

## 1. The orbit-product lemma

Let `K=C(t)`, let `Phi in K[X]` be separable and irreducible of degree `d`,
and let `Omega` be its roots in a splitting field `L`. Write

```text
C_Phi=product_(alpha in Omega) alpha.
```

Assume `C_Phi in C^*`. Let `S` be a nonempty proper subset of `Omega`, and
suppose its product

```text
p_S=product_(alpha in S) alpha
```

belongs to `K`. Then `p_S` is constant.

### Proof

Let `G=Gal(L/K)`. Irreducibility makes `G` transitive on `Omega`. Let

```text
O={sigma(S):sigma in G},             r=|O|.
```

Since `p_S in K`, every orbit subset has the same product:

```text
product_(alpha in T) alpha=p_S       for every T in O.       (1)
```

For a root `alpha`, let `eta_alpha` count the members of `O` containing it.
Transitivity makes `eta_alpha` independent of `alpha`; call the common value
`eta`. Double counting gives `d eta=r|S|`, so `eta>0`. Multiplying (1) over
the orbit gives

```text
p_S^r
 =product_(T in O) product_(alpha in T) alpha
 =product_(alpha in Omega) alpha^eta
 =C_Phi^eta in C^*.                                      (2)
```

A rational function whose positive power is a nonzero constant has no zero
or pole on `P^1`, hence is constant. This proves the lemma. Equivalently, (2)
is the norm of the subset product from its stabilizer field down to `K`.

For the application one needs even less: `p_S=c t` with `c!=0`. Taking the
`t`-adic valuation in (2) gives `r=0`, an immediate contradiction.

## 2. The root polynomial is irreducible and has constant total product

Let

```text
R(X)=r_0+r_1X+...+r_dX^d in C[X],
r_0 r_d!=0,
d=M+N,                 M,N>=1,
Phi(X)=X^M-tR(X).
```

Then `Phi` is irreducible over `C(t)`.

Indeed, first work in `C[X,t]`. Since `Phi` has degree one in `t`, any
nontrivial factorization has a factor `A(X)` independent of `t`. Comparing
the constant and linear coefficients in `t` shows that `A` divides both
`X^M` and `R(X)`. But

```text
gcd(X^M,R)=1                                             (3)
```

because `R(0)!=0`; hence `A` is a unit. Thus `Phi` is irreducible in
`C[X,t]`. It is primitive as a polynomial in `X` over `C[t]`, so Gauss's
lemma gives irreducibility in `C(t)[X]`. Characteristic zero gives
separability.

For `t!=0`, `Phi` has degree `d` in `X`, leading coefficient `-t r_d`, and
constant coefficient `-t r_0`. Vieta therefore gives

```text
product_(alpha in Omega) alpha=(-1)^d r_0/r_d in C^*.    (4)
```

Equations (3)--(4) are the two structural facts needed by the orbit-product
lemma: transitivity and a constant total norm.

## 3. Closing the constant-term theorem

Put

```text
Lambda(u)=u^(-M)R(u).
```

For small nonzero `t`, the equation `Phi(u)=0` has `M` roots tending to zero
and `N` roots tending to infinity. Let `S_0` be the small-root subset and

```text
Pi(t)=product_(u in S_0)u.
```

THM-1550 proves, by exact Wiener--Hopf factorization, that

```text
CT(Lambda^m)=0 for every m>=1
        iff
Pi(t)=c t,       c=(-1)^(N+d+1)r_0!=0.                 (5)
```

For completeness, its forward direction is short. Factor

```text
1-t Lambda(u)
 =(-t r_d)
   product_(u_i in S_0)(1-u_i/u)
   product_(a_j notin S_0)(u-a_j).                     (6)
```

The logarithms of the small-root factors contain only negative powers of
`u`; after writing `u-a_j=-a_j(1-u/a_j)`, the nonconstant part of each
large-root logarithm contains only positive powers. Hence

```text
CT_u log(1-t Lambda)
 =log((-t r_d)(-1)^N product_j a_j).                   (7)
```

The left side is
`-sum_(m>=1) CT(Lambda^m)t^m/m`. If all those coefficients vanish, the
argument of the logarithm in (7) is one. Substituting (4) and separating the
small-root product gives (5).

The small and large algebraic germs in (6) generate a splitting field of
`Phi` over `C(t)` inside a common ramified germ field. Thus the germ identity
`Pi=c t` is an equality in a splitting field, and `S_0` is a genuine
`M`-subset of `Omega`. Since

```text
0<M<d,                                                    (8)
```

the orbit-product lemma applies. Equations (4)--(5) would give

```text
(c t)^r=((-1)^d r_0/r_d)^eta,                            (9)
```

whose two sides have respective `t`-adic valuations `r>0` and `0`. This is
impossible. Therefore some `CT(Lambda^m)` is nonzero.

Now let `f in C[z,z^-1]` be arbitrary. If its constant coefficient is
nonzero, already `CT(f)!=0`. Otherwise, failure of strict one-sidedness means
that the support contains both a negative and a positive exponent; shifting
by the lowest exponent puts it in the form above. Consequently

```text
CT(f^m)=0 for every m>=1
  implies
supp(f) subset Z_(>0) or supp(f) subset Z_(<0).          (10)
```

This is exactly the contrapositive seed-existence statement consumed by
THM-2022.

## 4. What has and has not been removed

The THM-2022 paper proof may now cite THM-2067 instead of importing DvdK:
the lowest balanced face gives a non-one-sided one-variable Laurent
polynomial, and (10) supplies the nonzero face constant term used before
finite-field reduction.

This does **not** reproduce DvdK's stronger conclusion about a nonzero
critical value and

```text
limsup |CT(f^m)|^(1/m).
```

By itself this orbit contradiction gives no bound on the first nonzero `m`.
THM-2087 later adds the compound polynomial of every same-size root-subset
product and obtains the explicit bound

```text
m <= binom(M+N,a)+binom(M+N-1,a-1)-1,
a=min(M,N).
```

That estimate is generally non-sharp; the sharp Sturmfels/ESV target
`m<=M+N` remains open when `a>=2`. In particular, S222's general mixed-complex saddle
asymptotic and S223's assertion that mixed-sign cancellations are merely
finite are not needed here and remain unproved. The numerical-semigroup
return set controls which balanced words exist, but arbitrary complex
coefficients can cancel their total weight. The orbit norm handles those
cancellations without claiming eventual nonvanishing.

The Lean interface named `DvdK1` is not automatically discharged: THM-1550's
root-factorization criterion and this splitting-field argument would still
need formalization. The mathematical dependency, however, is now internal
and algebraic rather than a citation to the stronger analytic theorem.

## 5. The underlying object and the tournament audit

The useful vertices are algebraic root branches, not Laurent monomials or
balanced channels. The real carrier is the transitive `G`-set `Omega`
together with the orbit hypergraph `O` of the small-root subset. Uniform
incidence is what turns the product over hyperedges into the Vieta norm.

There is no faithful tournament reduction. Orienting two roots by valuation,
modulus, or label destroys the subset product and the orbit incidence count;
ties at ramified places are also structural. One may orient orbit subsets by
size or lexicographic labels to schedule checks, but the score sequence
preserves none of (2), (4), or (9). The challenged assumption is therefore
that pairwise root comparisons should carry noncancellation. They do not: the
proof is hypergraphic and multiplicative.

THM-2058 has the exact same uniform-incidence orbit norm for transported LRC
phase packets. The difference is now sharp rather than metaphorical: here
Vieta supplies a constant total product while the hypothesized subset product
has positive `t`-valuation; the LRC packet orbit currently has no analogous
nonconstant seed to oppose its uniform norm. That missing pair of quantities,
not a missing tournament orientation, is why the transfer stops.

## 6. Audit

The proof was checked independently against the following possible failure
points: primitivity in Gauss's lemma, separability, transitivity of the Galois
action, identification of the analytic small germs inside a splitting field,
nonvanishing of `c`, properness of the subset, and the incidence exponent.
All reduce to (3)--(8). No computation, genericity assumption, saddle
asymptotic, or coefficient sign condition enters the proof.
