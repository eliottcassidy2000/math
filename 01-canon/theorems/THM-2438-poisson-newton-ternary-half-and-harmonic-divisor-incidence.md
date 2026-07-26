---
id: THM-2438
title: "Poisson--Newton ternary half and harmonic divisor incidence"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT
  HOSTILE AUDIT. The exponential generating function of the forward
  differences of any sequence is its original exponential generating
  function multiplied by exp(-x). Consequently the Newton coefficients
  of the central binomial sequence are the central trinomial numbers,
  and the Newton coefficients of the two central Pascal halves count
  nonnegative and positive {-1,0,1}-walks. For any integer support H,
  its reciprocal sum is the extended limiting mean of its divisor
  incidence. It is finite exactly when the mean weak or strict
  multiplicative fibre loss caused by deleting H is finite, and in the
  convergent case both mean losses equal that reciprocal sum. This
  gives A014574 and Bertrand-series supports exact multiplicative-scar
  interpretations. It does not identify Cauchy and Dirichlet
  convolution or derive primes from additive startup defects.
source: mac-mini-2026-07-26-poisson-newton-harmonic-incidence
depends_on:
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - THM-2413-prime-index-affine-drift-and-twin-center-weld
  - THM-2433-operation-fibre-deletion-incidence-and-startup-scar
related:
  - THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry
script: 04-computation/poisson_newton_harmonic_incidence_thm2438.py
output: 05-knowledge/results/poisson_newton_harmonic_incidence_thm2438.out
script_sha256: 458b241b6580866b6f469a3d5cf52cd941f394b6f65cfb88092be26d3c1982a0
output_sha256: ade4f8795643d295169576f79ad04d00f3fe7a7b66ca6440d764d485cba94121
hash_basis: working-tree bytes (LF)
---

# THM-2438 -- Newton exposes ternary drift; harmonic mass measures the multiplicative scar

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; UNDER INDEPENDENT
HOSTILE AUDIT.**

Two analogies from the operation-graph program become exact after retaining
the right incidence coordinates:

```text
uniform additive step:
  sequence values --forward differences--> Newton coefficients;

nonuniform multiplicative step:
  support vertices --principal ideals--> divisor incidence.             (1)
```

The first transform sends the two central halves of a binary cube to the
two central halves of a ternary zero-drift walk.  The second sends the
reciprocal mass of an integer sequence to the mean multiplicative fibre
loss caused by deleting that support.  They are parallel incidence
statements, not an isomorphism between addition and multiplication.

## 1. The exact Poisson--Newton transform

Let `(a_n)_(n>=0)` take values in a characteristic-zero commutative
algebra and put

```text
d_k=Delta^k a_0
   =sum_(j=0)^k (-1)^(k-j) binom(k,j)a_j.                         (2)
```

As formal exponential generating functions,

```text
A(x)=sum_(n>=0) a_n x^n/n!,

D(x)=sum_(k>=0) d_k x^k/k!,

D(x)=exp(-x)A(x),                    A(x)=exp(x)D(x).             (3)
```

This is coefficientwise formal; no analytic convergence hypothesis is
needed.  Indeed,

```text
sum_(k>=0) sum_(j=0)^k
  (-1)^(k-j) binom(k,j)a_j x^k/k!

=sum_(j>=0) a_j x^j/j!
  sum_(m>=0)(-x)^m/m!

=exp(-x)A(x).                                                     (4)
```

Equation (3) is the sequence-level companion to THM-2412's polynomial
operator identity

```text
Delta=exp(D)-I.
```

It also makes the special role of two exact.  If `a_n=c^n`, then

```text
A(x)=exp(cx),                  D(x)=exp((c-1)x),

Delta^k a_0=(c-1)^k.                                             (5)
```

Thus `2^n` has the constant Newton spectrum `1`, whereas `4^n` has
the spectrum `3^k`.  The latter is not another unit eigenfunction.

## 2. The central binomial layer becomes central trinomial

Put

```text
C_n=binom(2n,n),

T_k=[u^0](u^(-1)+1+u)^k
   =sum_(j=0)^(floor(k/2)) k!/(j!^2(k-2j)!).                     (6)
```

The exponential generating functions are

```text
sum_(n>=0) C_n x^n/n! = exp(2x) I_0(2x),

I_0(2x)=sum_(j>=0) x^(2j)/(j!)^2.                                (7)
```

Multiplying (7) by `exp(-x)` and applying (3) gives

```text
sum_(k>=0) Delta^k C_0 x^k/k!
 =exp(x)I_0(2x)
 =sum_(k>=0) T_k x^k/k!.
```

Therefore

```text
Delta^k binom(2n,n)|_(n=0)=T_k.                                  (8)
```

The first values are

```text
T_k=1,1,3,7,19,51,141,393,1107,3139,... .                       (9)
```

This is an exact binomial-transform identity.  It is not a numerical
resemblance between two short prefixes.

## 3. The supplied Pascal halves are Newton shadows of ternary walks

Retain THM-2412's sequences

```text
A_n=sum_(j=0)^n binom(2n,j)
   =(4^n+C_n)/2
   =1,3,11,42,163,638,...,

M_n=(4^n-C_n)/2
   =0,1,5,22,93,386,...,

B_n=M_(n+1)
   =1,5,22,93,386,... .                                          (10)
```

Equations (5) and (8) imply

```text
Delta^k A_0=(3^k+T_k)/2,

Delta^k M_0=(3^k-T_k)/2,                                        (11)

Delta^k B_0
 =Delta^k M_0+Delta^(k+1)M_0
 =(4*3^k-T_k-T_(k+1))/2.                                       (12)
```

These coefficients have a literal walk interpretation.  Among the `3^k`
words in `{-1,0,1}^k`, exactly `T_k` have coordinate sum zero.  Negation
pairs every positive-sum word with a negative-sum word.  Hence

```text
(3^k+T_k)/2
 =#{epsilon in {-1,0,1}^k:sum epsilon_i>=0},

(3^k-T_k)/2
 =#{epsilon in {-1,0,1}^k:sum epsilon_i>0}.                      (13)
```

Thus the Gregory--Newton coefficients of the **binary** weak and strict
Pascal halves are the weak and strict halves of the **ternary** drift
alphabet.  Forward differencing does not preserve the original binary
alphabet: changing `n` adds two binary slots, and the new-versus-old
comparison has the three local states `-1,0,1`.

This sharpens the tournament analogy.  A tournament still has two
orientations per arc and therefore `2^binom(v,2)` labelled states.
Equation (13) describes the Newton transform of a two-slot binary layer;
it does not add a third orientation to a tournament.  If ties are admitted,
they are a genuine third state and must be retained rather than silently
identified with either orientation.

## 4. A support sequence as a multiplicative scar

Let `H` be any subset of the positive integers.  On the full positive
multiplicative carrier, let

```text
F_U^<=(n)=#{{a,b}:a<=b, ab=n},
F_U^<(n) =#{{a,b}:a< b, ab=n},

S=N minus H,

D_H^<=(n)=F_U^<=(n)-F_S^<=(n),
D_H^<(n) =F_U^<(n) -F_S^<(n).                                  (14)
```

Define three locally finite incidence counts:

```text
I_H(n)=#{h in H:h|n},

K_H(n)=#{{a,b}:a<b, a,b in H, ab=n},

Q_H(n)=1 if n=h^2 for some h in H, and 0 otherwise.              (15)
```

Then, pointwise for every `n`,

```text
D_H^<=(n)=I_H(n)-K_H(n),

D_H^<(n) =I_H(n)-K_H(n)-Q_H(n).                                 (16)
```

Indeed, a lost unordered factor pair with exactly one endpoint in `H`
contributes once to both `D` and `I`; a distinct pair with both endpoints
in `H` contributes once to `D` but twice to `I`, producing `K`; and a
diagonal hole pair contributes to the weak fibre and to `I`, but is absent
from the strict fibre, producing `Q`.

Equation (16) is the infinite-support form of THM-2433's Burnside deletion
formula.  Every target fibre is finite, so no global finiteness assumption
on `H` is required.

If `H` is finite with `m=|H|`, then for

```text
N>=max(H)^2
```

the cumulative losses are exactly

```text
sum_(n<=N) D_H^<=(n)
 =sum_(h in H) floor(N/h)-m(m-1)/2,

sum_(n<=N) D_H^<(n)
 =sum_(h in H) floor(N/h)-m(m+1)/2.                             (17)
```

For the proper-factor convention `a,b>=2` used in THM-2433, remove the
unit-cofactor occurrence for each hole.  The corresponding constants in
(17) become `m(m+1)/2` and `m(m+3)/2`; the limiting mean below is unchanged.

## 5. Harmonic mass is exactly the limiting mean scar

Summing the first count in (15) in the other order gives

```text
1/N sum_(n<=N) I_H(n)
 =sum_(h in H,h<=N) floor(N/h)/N.                                (18)
```

For each fixed `h`, the summand tends to `1/h` and is bounded by `1/h`.
Finite truncation followed by monotone/dominated convergence therefore
gives the extended limit

```text
lim_(N->infinity) 1/N sum_(n<=N) I_H(n)
 =sum_(h in H) 1/h,                                             (19)
```

where both sides may be `+infinity`.

The same phase transition, with the same finite limit, holds for the
actual weak and strict deletion losses:

```text
sum_(h in H)1/h < infinity

iff sup_N 1/N sum_(n<=N)D_H^<=(n)<infinity

iff sup_N 1/N sum_(n<=N)D_H^<(n)<infinity,                       (20)
```

and in the convergent case

```text
lim_(N->infinity) 1/N sum_(n<=N)D_H^<=(n)

=lim_(N->infinity) 1/N sum_(n<=N)D_H^<(n)

=sum_(h in H)1/h.                                               (21)
```

To prove (20)--(21), first note from (16) that

```text
D_H^<= <= I_H <= 2D_H^<=,

D_H^< <= I_H-Q_H <= 2D_H^<.                                    (22)
```

Also

```text
sum_(n<=N)Q_H(n)<=sqrt(N).                                      (23)
```

Thus a divergent reciprocal sum forces both mean losses to diverge.
Suppose instead that

```text
R_H=sum_(h in H)1/h<infinity.
```

Let `A_H(x)=#(H intersection [1,x])`.  Then `A_H(x)/x->0`.
The number of ordered pairs `(a,b) in H^2` with `ab<=N` is

```text
P_H(N)=sum_(a in H,a<=N) A_H(N/a).                              (24)
```

Given `epsilon>0`, choose `X` so that
`A_H(y)<=epsilon y` for `y>=X`.  Split (24) at `a=N/X`:

```text
P_H(N)/N
 <=epsilon sum_(a in H)1/a
   +sum_(a in H,a>N/X)1/a
 ->0.                                                           (25)
```

Hence the cumulative collision tax `K_H` is `o(N)`, while (23) handles
the diagonal.  Equations (16), (19), and (25) prove (21).

This is the rigorous meaning of “the reciprocal of an integer sequence is
a subset of the harmonic numbers” in the multiplicand graph: `hN` has
natural density `1/h`, so the reciprocal sum is the mean number of deleted
factor incidences.  Multiplication's base is nonuniform precisely through
these ideal densities.

For comparison, if `H` is finite, deleting it from the positive additive
fibre gives

```text
D_(+,H)^<=(z)=D_(+,H)^<(z)=|H|,          z>2max(H).               (26)
```

Every additive translate has eventual target density one, whereas the
multiplicative translate `hN` has density `1/h`.  This is an exact
structural asymmetry, not a claim that addition causes prime occurrence.

## 6. Abel--Dini and Bertrand forms of the criterion

Let

```text
A_H(N)=#(H intersection [1,N]).
```

Finite partial summation gives

```text
sum_(h in H,h<=N)1/h

=A_H(N)/N
 +sum_(n=1)^(N-1) A_H(n)/(n(n+1)).                              (27)
```

Equivalently, for `0<r<1`,

```text
sum_(h in H) r^h/h
 =integral_0^r (sum_(h in H)t^h) dt/t.                           (28)
```

All terms are nonnegative, so the limit as `r` increases to one is the
reciprocal sum, finite or infinite.  Equations (27)--(28) are the exact
Abel--Dini support criteria behind (19); no cancellation theorem is being
smuggled in.

If the increasing enumeration `(h_k)` satisfies

```text
h_k asymp k log(k)(log log(k))^alpha,                            (29)
```

then limit comparison and the substitution `u=log log x` give

```text
sum_k 1/h_k converges
iff alpha>1.                                                     (30)
```

Thus the same Bertrand-series threshold is exactly the threshold between a
finite and an infinite mean multiplicative scar.

The triangular support gives a sharp elementary control:

```text
T_k=k(k+1)/2,

sum_(k>=1)1/T_k
 =sum_(k>=1)2(1/k-1/(k+1))
 =2.                                                            (31)
```

Accordingly, a uniform integer has asymptotically two triangular divisors
on average, counted with support multiplicity.  By contrast
`1+1/2+1/3+1/4+1/5>2` records five dense startup ideals rather than one
quadratically sparse support.

## 7. A014574 acquires an exact mean-scar meaning

Let `C_twin` be A014574, the centers `c` for which `c-1,c+1` are prime.
THM-2413 proves the termwise identity

```text
1/c
 =1/2(1/(c-1)+1/(c+1))
  -1/(c(c^2-1)),                                                (32)
```

and, using Brun's theorem, the cited convergence formula

```text
R_twin:=sum_(c in C_twin)1/c
 =B_2/2-sum_(c in C_twin)1/(c(c^2-1))
 <infinity.                                                      (33)
```

Applying (19)--(21) gives the new operational interpretation

```text
lim_(N->infinity)
  1/N sum_(n<=N)#{c in C_twin:c|n}
 =R_twin,

lim_(N->infinity)
  1/N sum_(n<=N)D_(C_twin)^<=(n)
 =lim_(N->infinity)
  1/N sum_(n<=N)D_(C_twin)^<(n)
 =R_twin.                                                       (34)
```

Thus the reciprocal A014574 support is literally the mean multiplicative
scar made by deleting all twin centers.  This does not prove that there
are infinitely many twin primes: both a finite and an infinite support can
have finite reciprocal mass.  It also does not turn A014574 into the
A373813 line-cover sequence; THM-2413's affine-drift fibres remain a
different object.

## 8. Equality and failure boundaries

1. Equation (3) uses exponential generating functions.  Multiplying an
   ordinary generating function by `exp(-x)` is a different transform.
2. The ternary words in (13) are Newton coefficients of binary layers, not
   a bijection between the original sets and words at the same index.
3. `4^n` has difference eigenvalue three.  Only `2^n` has eigenvalue one.
4. The divisor mean counts support multiplicity.  Replacing a sequence by
   its set of values first removes collisions; THM-2352 explains why that
   distinction can change convergence thresholds.
5. For divergent reciprocal mass, the deletion-loss means diverge, but
   their difference from the marked incidence need not be `o(1)`: dense
   hole--hole collisions can carry a positive proportion.
6. Cauchy convolution on additive fibres and Dirichlet convolution on
   multiplicative fibres are not conjugate by the Stirling transform.
7. The min--max switching energy of THM-2412 retains cycle-product
   obstructions not present in the one-dimensional ternary walk count.

## 9. Exact companion

Run

```bash
python3 04-computation/poisson_newton_harmonic_incidence_thm2438.py
python3 -O 04-computation/poisson_newton_harmonic_incidence_thm2438.py
```

The dependency-free companion:

- verifies `300` Newton identities through order `59`;
- enumerates the ternary weak, strict, and zero layers through length `9`;
- checks the formal Poisson--Newton transform on `100` coefficients;
- independently checks the two Bessel coefficient expansions;
- audits (16) on all `128` hole subsets of `{1,...,7}` through target
  `120`;
- verifies the finite cumulative formula and additive tail for every
  nonempty such subset;
- checks `894` exact Abel partial-summation and divisor-incidence identities;
- verifies the triangular telescoping law through term `500`; and
- sieves all `8,169` A014574 centers through `1,000,000`, checking (32)
  and a direct marked-incidence control.

Every truth-bearing check raises explicitly under optimized Python.
Normal and optimized transcripts must byte-match the stored output after
LF normalization.

## 10. Audit boundary

The formal transform, central-trinomial coefficient extraction, ternary
reflection count, pointwise collision tax, finite summatory constants,
harmonic mean limit, and `o(N)` collision argument are proved above and
verified independently inside the exact companion.  Promotion beyond
candidate status still requires a separate hostile reader to reconstruct
the infinite-support limit and check that the full/proper factor
conventions are never mixed.  Brun convergence in (33) remains cited
through THM-2413.
