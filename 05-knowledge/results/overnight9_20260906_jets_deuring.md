# A divided reflected solution fixes the Deuring precision loss at one digit

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[independent proof and exact audit](overnight9_20260906_jets_independent_audit.md)
passes without requested repairs, including the adjacent-unit refinement.
This upgrades the eighth-round three-node Deuring ceiling test to an exact
largest-loss formula at every depth and every unit lift. It does not give
the full Smith partition or extend the general n-node coefficient packet.

## 1. Inheritance and exact statement

Use complete uniform Hasse jets at three distinct integer nodes, on the
integer coefficient module of degree below three times the multiplicity.
Let p be an odd prime and set

```
k=(p-1)/2,       m=k+1=(p+1)/2.
```

Assume all three pairwise p-adic depths equal e. After a translation and
a p-adic unit change of variable, the nodes are `0,p^e,p^e a`, with
`a,a-1` units in Z_p. Put

```
H_p(a)=sum_(j=0)^k binom(k,j)^2 a^j,
sigma=[H_p(a)=0 modulo p].
```

**Theorem.** If e=0, all Smith exponents are zero. For every e>=1 the
largest p-Smith exponent is exactly

```
L=(3m-1)e-sigma.                                      (1)
```

Thus the single Deuring residue bit determines the exact precision loss
through all higher unit lifts. No second digit is required for this target.
The determinant also gives the penultimate determinantal valuation

```
D_(3m-1)=(3m^2-3m+1)e+sigma,       D_(3m)=3m^2e.       (2)
```

These two ideals do not determine all intermediate factors. The full p=7
equilateral four-jet list remains the separate eighth-round theorem.

The closest proved mechanism is the
[eighth Deuring coefficient theorem](overnight8_20260906_jets_residue.md),
which gave `L=(3m-1)e` precisely when H_p(a) is nonzero modulo p and only
an upper bound when it vanishes. The exact inverse-denominator equality
is the proved
[THM-4443](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md).
The canonical hostile is to retain one top numerator alone: its simple
residue roots lift to make that numerator divisible by arbitrarily high
p-powers. The corrected near miss is to treat this individual cancellation
as a simultaneous precision loss. The least-used sidecar is the divided
difference between a polynomial solution and its reflected solution.

The concept board is: top-jet denominators; simultaneous companion ideals;
Deuring residues; a reflected hypergeometric equation; the Wronskian
constant; and the competing adjacent jet. Initial exact p^2 probes at
p11,p19 and other small primes found that every lift killing the first
numerator failed to kill its companion. The argument below explains this
for all primes without extrapolating the finite bank.

## 2. Exact numerator normalization

Define integer polynomials

```
P(a)=sum_(j=0)^k binom(k+j,j) binom(2k-j,k-j) a^j,
R(a)=(-1)^k P(1-a),
T(a)=a^k P((a-1)/a).                                  (3)
```

The expression T is a polynomial: its jth summand is an integer multiple
of `a^(k-j)(a-1)^j`. If the two normalized differences at a node are u,v,
the top reciprocal coefficient is

```
(-1)^k sum_j binom(k+j,j)binom(2k-j,k-j) u^(k-j)v^j
                                      / (uv)^(m+k).   (4)
```

At the three nodes 0,1,a, its numerators are respectively P, R, and
`(-1)^k T`; all denominators are units. Thus their simultaneous minimum
valuation is exactly the minimum for P,R,T. The common scale p^e supplies
the denominator exponent `(2m+k)e=(3m-1)e` at top order.

The elementary eighth-round binomial identities give

```
P(a)=(-1)^k H_p(a) modulo p,
P(a)=(-1)^k P(1-a) modulo p,
T(a)=(-1)^k P(a) modulo p.                            (5)
```

For example, the two binomial factors in (3) reduce respectively to
`(-1)^j binom(k,j)` and `(-1)^(k-j) binom(k,j)`. Reflection follows from
translation of the coefficient `[X^(p-1)](X(X-1)(X-a))^k`; its degree is
below2p-1. Reciprocity of P is exact, since its coefficients are symmetric,
and then gives the third identity in (5).

Consequently all three top numerators are units together in the ordinary
case, or are all divisible by p in the Deuring-zero case. The question is
whether all can gain a second digit. It is enough to examine P and R.

## 3. The integral divided companion and its Wronskian

By (5),

```
G(a)=(P(a)-R(a))/p                                   (6)
```

lies in Z[a]. P and R have the same leading coefficient
`binom(2k,k)`, so `deg G<=k-1`. This degree drop is essential.

Let

```
L0=a(1-a)D^2+(1-2a)D+k(k+1),       D=d/da.
```

There are exact characteristic-zero identities

```
L0 P=p P',
L0 R=-p R',
L0 G=P'+R'.                                          (7)
```

Here is a direct proof, requiring no external differential-equation
theorem. If c_j is the jth coefficient of P, then for 0<=j<k,

```
(j+1)(2k-j)c_(j+1)=(k-j)(k+j+1)c_j.
```

This is precisely the coefficient recurrence for
`a(1-a)P''+(-2k-2a)P'+k(k+1)P=0`. Since 2k=p-1, this is the first
identity in (7). Replacing a by1-a gives the second. Subtraction and
division by the integer p gives the third.

Reduce modulo p. Then `R=P`, `L0 P=0`, and `L0 G=2P'`. For
`W=P G'-P'G`, the differential equation gives

```
(a(1-a)W)'=P L0G-G L0P=2P P'=(P^2)' modulo p.        (8)
```

The polynomial `a(1-a)W-P^2` has degree at most2k=p-1. A polynomial
of degree below p with zero derivative is constant in F_p[a]. This step
would be invalid without the degree bound. Evaluating at a=0 gives the
constant `-P(0)^2=-1` modulo p, because
`P(0)=binom(p-1,k)=(-1)^k` modulo p. Therefore

```
a(1-a)(P G'-P'G)=P^2-1 modulo p.                    (9)
```

At every Deuring root a distinct from zero and one, (9) becomes

```
a(1-a)P'(a)G(a)=1 modulo p.                          (10)
```

In particular both P' and G are units there. This simultaneously proves
simplicity of the residue root and the required nonvanishing companion.
For every p-adic lift above it, the difference P(a)-R(a)=pG(a) has valuation
exactly one. Since both P and R are already divisible by p, their minimum
valuation must be exactly one. T is also divisible by p by (5). Hence

```
min(v_p(P(a)),v_p(R(a)),v_p(T(a)))=sigma.              (11)
```

This is an all-lift statement. Hensel cancellation of P alone cannot
raise the common ideal beyond one digit. No claim about valuations of
individual numerators has been substituted for (11).

## 4. Exact precision and the adjacent-order check

All normalized reciprocal coefficients are p-adically integral, since
the three normalized differences are units. At e>=1 every lower observed
order l<k contributes at most `(2m+k-1)e`. The actual largest top-order
contribution is `(3m-1)e-sigma` by (11), and sigma is at most one.
Because e>=1, it dominates or ties every lower order. THM-4443's attained
inverse-denominator equality proves (1). At e=0 the confluent determinant
is a unit, giving the separate zero list. Formula (2) follows from the
determinant and the exact largest factor, without an inference about the
remaining partition.

There is also a positive adjacent-jet fact, useful for checking the
competition rather than merely bounding it. Put

```
A(a)=[X^k]((X-1)(X-a))^k=(-1)^k H_p(a),
B(a)=[X^(k-1)]((X-1)(X-a))^k.
```

The exact polynomial identity is

```
(k+1)B=-a(kA+(1-a)A').                               (12)
```

To derive it, regard A as a homogeneous degree-k polynomial in the two
roots c,a, and differentiate both homogeneity and simultaneous translation.
At c=1 these give `A_c=kA-aA'` and
`A_c+A'=-(k+1)[X^(k+1)]((X-1)(X-a))^k`.
The latter coefficient is B/a, by the reciprocal identity interchanging
the roots1,a under X->a/X. This proves (12).

At a Deuring root, k+1=-k modulo p and A'=P' is a unit by (10).
Equation (12) implies that B is a unit. Translating the cubic raised to k,
its coefficient of degree p-2 becomes `B-A x_i`; no other degree can
contribute modulo p in the degree range3k. When A=0, this is B at every
endpoint. The same Frobenius reciprocal argument as in the eighth report
then shows that **all three normalized adjacent coefficients of order
k-1 are units**. Their actual contribution is exactly `(3m-2)e`; it ties
the top contribution precisely at e=1 on the Deuring-zero branch.

This includes p=3,k=1, where the adjacent order is simply order zero.
In that case P=2+2a, R=2a-4, G=2, so there is no small-prime exception
hidden in the divided-companion argument.

## 5. Scope, controls, and the repaired question

The [eighth report](overnight8_20260906_jets_residue.md)
already supplies the primary-source identification of Deuring zeros with
supersingular Legendre parameters; that interpretation is inherited here
and is not a proof input to (1). The new map keeps the integral top-numerator
pair, passes to a divided reflection, and uses its Wronskian constant to
recover the simultaneous valuation lost by reducing to H_p. The result
upgrades a ceiling test to exact precision because the missing cap is one.
It does not reconstruct every unit-sensitive minor of a different observer.

The quantifier `m=(p+1)/2` is mandatory. This is not a theorem for every
prime at one fixed arbitrary multiplicity, for non-equilateral triples,
or for the general n-node coefficient packet. The residue-field extension
and higher-genus versions are not needed or claimed. Full intermediate
Smith factors outside the already proved p7 family remain a separate
question.

The standalone source checks exact integer coefficient recurrences,
all three differential identities, the strict degree bounds, the entire
Wronskian congruence, and (12). It tests sixteen named odd primes from3
through103, all118 Deuring root classes in that bank, and the unique first
Hensel digit that kills P at each root. The companion R survives at that
digit in every case. At primes at most31 it additionally checks every
lift above every root: exactly600 lifts, namely
`3+21+33+57+207+279`. The all-prime theorem rests on (7)-(11), not this
finite bank.

A separate rational product engine reconstructs the actual reciprocal
jets, including the adjacent coefficient at each node. It checks76
literal precision instances across four depths, signed denominator
positions, ordinary classes, zero classes, and higher unit lifts. Eighteen
full integer Hasse Smith matrices at p3,5,7,11 independently verify the
largest exponent, determinant, and penultimate determinantal ideal. The
source imports no repository or earlier producer code.

The strongest repaired question is now about the **intermediate** ideals:
can a similar divided companion control their joint unit cancellations?
For this precise largest-loss observer, the higher-digit question is
closed by the one-digit cap. A proposed higher-genus or arbitrary-node
extension must recheck both its differential identities and its degree
bound; zero derivative alone no longer forces a constant above degree p.

## 6. Reproduction and freeze

```
python -B 04-computation/overnight9_20260906_jets_deuring.py
python -B -O 04-computation/overnight9_20260906_jets_deuring.py
```

The source passes **1,292 optimization-live gates**. Normal and optimized
transcripts are byte-identical and LF-only:

```
source e7e4846f4bdd9cb5b2174571365192436c1578bd110dd71f6dde41a53faa3ab6
output 13b6475b7110129aee0e7c92bafb6aa995465fc3e153e965e1cfe176eb576f08
```

The independent referee derives the coefficient recurrence, both reflected
equations, the divided-companion integrality and degree, and the Wronskian
constant separately. It checks13 primes, all2,255 admissible p^2 lifts
through31 (both zero and ordinary classes),104 actual reciprocal observers,
and20 full local Smith matrices, using a separate standard-library engine.
The full prose and all three numerator normalizations are accepted,
including p3, the e0 mask, and the adjacent coefficient identity. Its
normal and optimized runs pass7,269 live gates with identical output:

```
audit source ce5b80b93a3fa39650e4cb01b6eaf92223476ea4eff012b147df501ef774c29a
audit output 0c8b33e3a0e421fee84c8c809fa72f594b2654a317015f82cc8cdfbc8322a1c3
```

All files remain outside the repository; no shared edits or Git mutation.

**Filing:** root integrated these audited artifacts in the ninth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
