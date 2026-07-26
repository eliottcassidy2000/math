---
id: THM-2448
title: "Finite-hole relative atoms, prime rays, and cover bounds"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every
  finite nonempty hole set H in the multiplicative carrier {2,3,...},
  artificial weak and strict atoms have an exact valuation-box
  criterion. If H is prime-free they form finite sets bounded sharply
  by max_(h in H) h*spf(h)<=M^(3/2). If H contains a prime, an exact
  finite cage C_H generates pairwise disjoint prime rays contained in
  both atom sets; after removing those rays, the weak and strict
  remainders are sharply bounded by M^3 and M^4. Their Dirichlet series
  therefore split into a finite cage polynomial times a prime tail
  plus a finite polynomial. Additive artificial atoms instead stop
  sharply at 2L+1 and 2L+2. Composite labels in the induced
  divisibility Hasse diagram obey a sharp valuation-free M^2 bound,
  and prime-free artificial twin centres obey c/2,c/3 in H and c<=2M.
source: codex-2026-07-26-finite-hole-relative-atoms
depends_on:
  - THM-2433-operation-fibre-deletion-incidence-and-startup-scar
  - THM-2000-support-harmonic-abel-dini-figurate-surface
related:
  - THM-362-natural-operation-graph-shadows
  - THM-2413-prime-index-affine-drift-and-twin-center-weld
  - THM-2438-poisson-newton-ternary-half-and-harmonic-divisor-incidence
script: 04-computation/finite_hole_relative_atoms_prime_rays_thm2448.py
output: 05-knowledge/results/finite_hole_relative_atoms_prime_rays_thm2448.out
script_sha256: 7ab8a72fc292e66340d7b5cffbfe15ccbb34ff401986ef670ddb2c109a72222d
output_sha256: 41e89b4ddf6bcae57ded4e3f9020231b8500770eed3539cd580e0b1e8ad5bd93
hash_basis: working-tree bytes (LF)
---

# THM-2448 -- finite holes create prime rays, not arbitrary multiplicative noise

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2433 computes the complete startup example `H={4,6}`.  The uniform
picture is sharper:

```text
prime-free holes:
  every artificial multiplicative atom is a bounded finite accident;

a prime among the holes:
  a finite divisor cage emits exact, pairwise disjoint rays m*q
  along all sufficiently large primes q;

additive holes:
  no rays occur -- every artificial atom remains in a linear startup
  interval.                                                        (1)
```

This separates the nonuniform divisor geometry of multiplication from the
uniform tail geometry of addition.  It does not derive primes from additive
startup defects.

## 1. Artificial atoms and their valuation boxes

Let

```text
H be a finite nonempty subset of {2,3,...},
S={2,3,...} minus H,                         M=max H.     (2)
```

For a retained target `z in S`, let its ambient weak proper-factor fibre
consist of unordered pairs

```text
{a,b},       2<=a<=b,       ab=z,                         (3)
```

and let the strict fibre require `a<b`.  Define:

```text
W_H={retained z:
       the ambient weak fibre is nonempty
       but its subfibre with a,b in S is empty},

T_H={retained z:
       the ambient strict fibre is nonempty
       but its subfibre with a,b in S is empty}.          (4)
```

Thus `W_H,T_H` contain only deletion-created atoms.  Ambient primes, and
ambient prime squares in the strict case, are deliberately excluded.

Write

```text
z=product_(i=1)^r p_i^(alpha_i),

B_alpha=product_i {0,1,...,alpha_i},

d_beta=product_i p_i^(beta_i).                           (5)
```

The involution `beta -> alpha-beta` swaps the two co-factors.  Then

```text
z in W_H
iff
  z notin H, z is composite, and
  for every beta in B_alpha minus {0,alpha},
     d_beta in H or d_(alpha-beta) in H;                 (6)

z in T_H
iff
  z notin H, z is neither prime nor a prime square, and
  (6) holds after all fixed points 2beta=alpha are omitted. (7)
```

**Proof.** The proper exponent box parametrizes the ordered proper-factor
fibre.  Quotienting by its involution gives (3); its fixed points are exactly
the diagonal square factorizations.  A retained pair disappears exactly when
at least one of its two vertices lies in `H`.  The ambient nonemptiness
conditions give the composite and prime-square clauses. QED.

The box is the full structural coordinate.  Looking only at the least-prime
pair is necessary but not sufficient: for `H={8,12}`, `z=24`, the pairs
`2*12` and `3*8` hit holes but `4*6` survives.

## 2. Prime-free holes have a sharp finite bound

Suppose `H` contains no prime.  Put

```text
B_H=max_(h in H) h*spf(h),                               (8)
```

where `spf` denotes the least prime factor.  Then every
`z in W_H union T_H` has the form

```text
z=p*h,        h in H,        p=spf(z)<=spf(h),           (9)
```

and consequently

```text
z<=B_H<=M^(3/2),                 Omega(z)>=3.             (10)
```

Both bounds are sharp: for `H={p^2}`, the target `z=p^3` belongs to both
`W_H` and `T_H` and satisfies equality in (8)--(10).

**Proof.** In the least-prime factor pair

```text
z=p*(z/p),
```

the prime `p` is not a hole.  Equations (6)--(7) therefore force
`h=z/p in H`.  Every prime factor of `h` is at least `p`, proving
`p<=spf(h)`.  A prime-free `h` is composite, so
`spf(h)<=sqrt(h)` and `Omega(h)>=2`.  This proves (10).  The displayed
singleton family proves sharpness. QED.

This already explains why the THM-2433 startup deletion has only finitely many
artificial atoms: `{4,6}` is prime-free.

## 3. The finite cage and its disjoint prime rays

For arbitrary `H`, define the divisor cage

```text
C_H={m>=2: every nonunit divisor of m belongs to H}.      (11)
```

In particular `C_H subset H`, so it is finite.  Moreover,

```text
C_H is nonempty iff H contains a prime.                   (12)
```

Indeed a prime in `H` belongs to `C_H`; conversely every `m in C_H` drags
each of its prime divisors into `H`.

Define the remote ray set

```text
R_H={m*q: m in C_H, q>M prime}.                           (13)
```

Then

```text
R_H subset W_H intersection T_H,                          (14)
```

and the products in (13) are pairwise distinct:

```text
m_1 q_1=m_2 q_2       implies       (m_1,q_1)=(m_2,q_2). (15)
```

**Proof.** Consider a proper factorization of `m*q`.  Since `q>M>=m`,
exactly one side contains `q`; the other side is a nonunit divisor of `m`
and hence lies in `H`.  Thus every weak or strict pair is killed.  The pair
`m<q` shows that the ambient strict fibre is nonempty.  For (15), `q_1`
cannot divide `m_2<=M`, so it must equal `q_2`, after which `m_1=m_2`.
QED.

The cage need not contain only primes.  For

```text
H={2,3,6},
```

one has `6 in C_H`, so every `6q` with prime `q>6` lies in both atom sets.
This is the minimal hostile control against a semiprime-only ray picture.

## 4. Everything off the rays is sharply bounded

The remote decomposition is exact:

```text
W_H minus R_H subset [2,M^3],

T_H minus R_H subset [2,M^4].                            (16)
```

Both exponents are sharp.  For `H={p}`:

```text
p^3 in W_H minus R_H,                   p^3=M^3,

p^4 in T_H minus R_H,                   p^4=M^4.         (17)
```

The target `p^4` is not weakly artificial because the diagonal
`p^2*p^2` survives.  This is exactly the weak/strict square wall.

**Proof of the weak bound.** Factor `z` into an increasing list of primes
with multiplicity,

```text
z=p_1 p_2 ... p_k,        p_1<=...<=p_k,                (18)
```

and let `a_j=p_1...p_j`.  Suppose `z in W_H` and `z>M^3`.  Let `j` be the
first index with `a_j>M`.

If `j<k` and `p_j>M`, a proper factor pair can be chosen with both sides
larger than `M`, contradicting (6).  Hence `p_j<=M`, so

```text
M<a_j<=M^2,             z/a_j>M,
```

again giving a surviving proper pair.  Therefore `j=k`.  Put
`m=a_(k-1)<=M` and `q=p_k`; then `q>M`.  For every nonunit divisor `d|m`,
the complementary factor `(m/d)q` exceeds `M`.  Equation (6) forces
`d in H`, so `m in C_H` and `z=mq in R_H`, a contradiction.  Thus the
first inclusion in (16) holds.

**Proof of the strict bound.** Suppose instead `z in T_H` and `z>M^4`.
Use the same first crossing.  If `p_j>M` and `j<k`, the factor
`p_j` has a distinct complementary factor larger than `M`; equality could
occur only for the excluded prime-square ambient fibre.  Otherwise

```text
M<a_j<=M^2<z/a_j,
```

so the two retained sides are distinct.  Hence again `j=k`, and the same
divisor argument gives `m in C_H` and `z in R_H`.  This proves the second
inclusion.  The factorizations in (17) prove sharpness. QED.

## 5. Exact Dirichlet phase transition

For `Re(s)>1`, put

```text
C_H(s)=sum_(m in C_H)m^(-s),

P_>M(s)=sum_(q>M prime)q^(-s).                           (19)
```

Pairwise ray disjointness and (16) give the exact decompositions

```text
sum_(z in W_H) z^(-s)
 =C_H(s)P_>M(s)+E_W(s),

sum_(z in T_H) z^(-s)
 =C_H(s)P_>M(s)+E_T(s),                                 (20)
```

where `E_W` and `E_T` are finite Dirichlet polynomials supported in
`[2,M^3]` and `[2,M^4]`, respectively.

Consequently:

- `H` is prime-free iff `W_H,T_H` are finite;
- in that finite phase both reciprocal sums converge (and the support series
  are Dirichlet polynomials);
- if `H` contains a prime, both reciprocal sums diverge and both support
  series have abscissa of absolute convergence `1`.

The last assertion follows from (12), any one fixed cage ray, and the
divergence of the reciprocal-prime sum.  In the repository's nonnegative
support-phase convention the finite/Dirichlet-polynomial phase is labelled
abscissa `0`; classically a finite Dirichlet polynomial has abscissa
`-infinity`.  This convention boundary is explicit rather than silently
identifying the two usages.

Equation (20) is the finite-hole counterpart of THM-2438's theorem that
harmonic support mass is limiting mean divisor incidence.  A prime hole
creates a positive finite cage coefficient multiplying the prime harmonic
boundary; a prime-free hole module cannot.

## 6. Addition has only a linear startup interval

Let now `H_+` be a finite nonempty subset of the positive integers,
`L=max H_+`, and define artificial weak and strict **summand** atoms exactly
as in (4), using positive pairs `a+b=z`.

Then

```text
weak artificial summand atom  -> z<=2L+1,

strict artificial summand atom -> z<=2L+2.              (21)
```

Both cutoffs are sharp for `H_+={1,2,...,L}`.

**Proof.** At `z>=2L+2`, the weak pair

```text
(L+1, z-L-1)
```

has both entries outside `H_+`.  It becomes strict for `z>=2L+3`.  At the
two preceding endpoints, every weak, respectively strict, pair has its
smaller entry at most `L`; the interval hole set kills all of them. QED.

Thus additive irregularities propagate into labelled fibre deficits, as in
THM-2433, but cannot emit atom rays.  Multiplicative rays arise from divisor
incidence, not from a causal conversion of addition into multiplication.

## 7. Composite induced covers obey a sharp `M^2` law

Return to the multiplicative `S` in (2), but now forget the co-factor and
take the divisibility order induced on `S`.  For retained `x<z` with `x|z`,
put `q=z/x`.  Then

```text
x prec_S z
iff
  x*d in H for every d|q with 1<d<q.                    (22)
```

In particular prime `q` gives the generic cover.  If `q` is composite and
(22) holds, then

```text
x*z<=M^2,                    z<=M^2/2.                   (23)
```

Both inequalities are sharp for

```text
H={6},             x=2,             z=18.               (24)
```

**Proof.** Every intermediate divisor vertex is uniquely `xd` with
`1<d<q`, proving (22).  If `p|q` is prime, both `p` and `q/p` are proper
divisors (possibly equal), so

```text
xp in H,             z/p=x(q/p) in H.
```

Their product is `xz`, proving the first inequality; `x>=2` proves the
second.  In (24) the sole intermediate vertex is `6`, and both bounds are
equalities. QED.

For `H={4,6}`, criterion (22) recovers exactly THM-2433's four exceptional
composite-labelled covers:

```text
(2,8;4),       (2,12;6),       (3,12;4),       (2,18;9). (25)
```

## 8. Twin centres meet the hole module in two coordinates

Let `C_all` denote the full A014574 carrier:

```text
C_all={c>=1: c-1 and c+1 are prime}.                    (26)
```

Suppose `H` is prime-free and

```text
c in C_all intersection (W_H union T_H).
```

Then

```text
c/2 in H,              c/3 in H,              c<=2M.   (27)
```

**Proof.** The exceptional centre `4` and the centre `6` cannot be
artificial under prime-free deletion: their proper factors are prime.  Every
later twin centre is divisible by six.  Applying (6) or (7) to the pairs
`2*(c/2)` and `3*(c/3)` forces the nonprime co-factors into `H`; the first
one gives `c<=2M`. QED.

For the startup module,

```text
H={4,6}:              W_H=T_H={8,12},                  (28)
```

and `12` is the unique artificial A014574 centre.  If one classifies all
targets with empty strict retained fibre rather than the artificial set
`T_H`, ambient primes and prime squares must be added separately.

## 9. Graph interpretation and hostile boundaries

For an internal operation shadow

```text
x R_odot z       iff       z=x odot s for some retained co-input s,
```

associativity sends a two-step path with witnesses `s,t` to the composite
witness `s odot t`.  The shortcut carrying that **labelled composite
witness** exists exactly when `s odot t` is retained.  For positive addition
and multiplication, cancellation makes this the corresponding iff for that
fixed input and output.  In a general associative operation, some different
retained witness can produce the same collapsed shortcut, so associativity
alone gives only the forward labelled-composition statement.  This is the
precise transitivity-like feature noticed in the summand and multiplicand
graphs.

It is not a tournament statement: the relation need not orient every pair,
and forgetting the labelled co-input loses the valuation-box cage.  The
following hostiles mark the essential losses:

```text
H={2,3,6}:    the composite cage element 6 emits 6q rays;

H={8}, z=16: the diagonal 4*4 survives weakly but is omitted strictly;

H={8,12}, z=24:
              killing 2*12 and 3*8 does not kill the surviving 4*6 pair.
                                                               (29)
```

## 10. Exact companion

Run

```bash
python3 04-computation/finite_hole_relative_atoms_prime_rays_thm2448.py
python3 -O 04-computation/finite_hole_relative_atoms_prime_rays_thm2448.py
```

The dependency-free verifier:

- performs `2,043,904` independent direct-fibre/valuation-box comparisons
  over all `2^11` hole sets in `{2,...,12}` and targets through `500`;
- performs `5,370,000` weak/strict prime-ray and finite-remainder checks on
  a `179`-set hole atlas through target `15,001`;
- checks `179` exact rational prime-ray prefix products;
- checks `40,290` additive hole/target cases and both sharp endpoints;
- checks `65,173` induced-cover cases and hostile controls;
- sieves all `8,169` twin centres through `1,000,000` and recovers `12`
  uniquely for `H={4,6}`; and
- verifies the sharp `p^3`, `p^4`, `M^(3/2)`, and `M^2` examples plus every
  hostile in (29).

Every assertion uses an optimization-safe `require`.  Normal and optimized
transcripts byte-match the stored output.
