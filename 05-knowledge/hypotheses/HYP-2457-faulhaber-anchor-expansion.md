# HYP-2457 - Faulhaber anchor expansion is an odd-moment triangular carrier

**Status:** OPEN synthesis; exact odd-moment identity and formal expansion
verified, uniform remainder/bracket proof still open.
**Source:** codex-2026-06-13.
**Companions:** HYP-2456, HYP-2455, HYP-2454, HYP-2453, HYP-2444,
HYP-2443, HYP-2239, HYP-2128.
**Computation:** `04-computation/triangular_faulhaber_anchor_expansion_codex.py`;
stored output `05-knowledge/results/triangular_faulhaber_anchor_expansion_codex.out`.

## Statement

Let `a_p(n)` be the real anchor solving

```text
sum_{j=0..n} (a+j)^p = sum_{j=1..n} (a+n+j)^p.
```

Set

```text
c = a+n,
u = n(n+1),
S_r(n)=sum_{i=1..n} i^r.
```

Then the balance equation is exactly

```text
D_p(c,n)
 = sum_{i=0..n}(c-i)^p - sum_{i=1..n}(c+i)^p
 = c^p - 2 * sum_{r odd} binom(p,r)c^(p-r)S_r(n).
```

Thus all even Faulhaber moments cancel.  The real root has the formal fixed-`p`
expansion

```text
c_p(n) = p*u + alpha_p + beta_p/u + gamma_p/u^2 + O(u^-3),
```

where

```text
alpha_p = (p-1)(p-2)/(12p),

beta_p  = -(p-1)(p-2)(2p^2-4p-1)/(180p^3),

gamma_p = (p-1)(p-2)(p+2)(40p^3-128p^2+53p+11)/(30240p^5).
```

Subtracting `n` gives the user's expansion:

```text
a_p(n)
= p*n^2 + (p-1)*n
  + (p-1)(p-2)/(12p)
  - (p-1)(p-2)(2p^2-4p-1)/(180p^3*n(n+1))
  + O(n^-4).
```

The factor `(p-1)(p-2)` in every displayed correction explains why `p=1` and
`p=2` are exact integer towers, while higher powers drift into a Bernoulli
fractional address.

## Exact p=1 And p=2 Recovery

For `p=1`,

```text
c_1(n)=u=n(n+1),
a_1(n)=c_1(n)-n=n^2.
```

This is the ordinary equal-sum tower from HYP-2454.

For `p=2`,

```text
c_2(n)=2u=2n(n+1),
a_2(n)=c_2(n)-n=2n^2+n.
```

This is the consecutive-square tower from HYP-2454.

The correction terms vanish for both `p=1` and `p=2`:

```text
alpha_p=beta_p=gamma_p=0.
```

This is now a structural fact, not a low-degree accident.

## Derivation Sketch

Pair the midpoint terms:

```text
(c-i)^p-(c+i)^p
 = sum_r binom(p,r)c^(p-r)((-i)^r-i^r).
```

Even `r` vanish.  Odd `r` contribute `-2i^r`, and the unpaired middle term
`c^p` remains.  This gives the odd-moment equation.

Now write `x=c/u` and `q=1/u`.  The first relevant Faulhaber moments are

```text
S_1 = u/2,
S_3 = u^2/4,
S_5 = u^3/6 - u^2/12,
S_7 = u^4/8 - u^3/6 + u^2/12.
```

Through `q^3`, the reduced equation is

```text
x^p - p*x^(p-1)
- binom(p,3)x^(p-3)q/2
- binom(p,5)x^(p-5)q^2/3
+ binom(p,5)x^(p-5)q^3/6
- binom(p,7)x^(p-7)q^3/4
=0.
```

Plugging

```text
x = p + alpha_p*q + beta_p*q^2 + gamma_p*q^3
```

makes the coefficients of `q^0,q^1,q^2,q^3` vanish.  The stored script checks
this exactly for `p=1..12`; the same algebra is a finite symbolic calculation
for each fixed `p`.

## Square-Pyramidal Cuboid Carrier

At `p=2`, Faulhaber's formula gives the square pyramidal numbers

```text
P_2(n)=1^2+...+n^2=n(n+1)(2n+1)/6.
```

Hence

```text
6P_2(n)=n(n+1)(2n+1).
```

Equivalently, six equal square pyramids, or one double-sized plus four regular
ones in the user's packing picture, fill the cuboid of dimensions
`n`, `n+1`, and `2n+1`.

HYP-2454's ordinary tower common sum is

```text
S_1(n)=n(n+1)(2n+1)/2 = 3P_2(n),
```

so the cuboid is also `2S_1(n)`.  This makes the p=2 tower a geometric
packing face of the same triangular carrier `u=n(n+1)`.

## Link Back To HYP-2454

HYP-2454 observed the finite bracket

```text
D_p(2pT_n,n)<0<D_p(2pT_n+1,n)
```

for `3<=p<=8` and `n<=40`, where `2pT_n=p*n(n+1)=p*u`.

HYP-2457 sharpens the route: at `c=p*u`, the leading odd-moment terms force a
negative defect, and the root moves by

```text
alpha_p = (p-1)(p-2)/(12p),
```

which lies in `(0,1)` for the checked higher powers.  The next proof target is
to turn the asymptotic drift into a uniform sign bracket for all `n` and
`p>=3`, or to isolate the first small exceptional case.

## LRC14 Transfer

This is another instance of the repo's hidden-address grammar:

```text
visible scalar balance -> attach the missing odd-moment/Faulhaber address
                        -> recover the exact obstruction.
```

For LRC14, a denominator being "blocked" is analogous to the scalar anchor
`c`.  The proof object should keep the odd-wall/resource ledger:

```text
blocked twists,
owner support,
divisor fiber,
shell-27 class,
carry residue,
endpoint atom,
moment/resource defect.
```

HYP-2456 showed that a visible triangular crossover word becomes exact only
after adding a Beatty/Pell carry coordinate.  HYP-2457 says the power-balance
anchor behaves similarly: the integer tower is the visible carrier, while
Bernoulli/Faulhaber odd moments are the hidden address deciding when integrality
survives.

## Tournament Analysis

The stored tournament deliberately does not use interval entries as vertices.
Alternatives considered: row entries, endpoints, centers, odd moments,
Bernoulli denominators, cuboid packings, simplex-packing analogues, LRC proof
obligations, and support lifts.

The chosen vertices are:

```text
odd_moment_projection,
u_anchor_expansion,
p1_p2_exact_towers,
square_pyramid_cuboid,
simplex_cuboid_packing,
bernoulli_denominator_address,
lrc14_moment_resource,
support_lift_transfer.
```

The pairwise observable is a majority comparison of

```text
(exactness, proof_readiness, geometry, lrc_transfer, computability, novelty),
```

with the listed order as the fixed Hamiltonian tie path.  The run is transitive:

```text
leader = odd_moment_projection
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles = 0
SCC sizes = [1,1,1,1,1,1,1,1]
Hamiltonian paths = 1
edge_flips_vs_exactness = 3
```

The useful challenged assumption is that the scalar anchor is not the proof
object.  The proof object is the retained odd-moment/support address explaining
when the anchor stays integral and when it drifts into Bernoulli fractions.

## Next Tasks

1. Prove the fixed-`p` asymptotic with a uniform `O(u^-3)` remainder after
   `gamma_p`, using the exact odd-moment equation and the implicit-function
   expansion around `x=p`.
2. Prove the HYP-2454 integer bracket for every `p>=3`, or split out the finite
   small-`n` exceptional check.
3. Compare the simplex-in-cuboid packing constants across dimensions with the
   Faulhaber odd-moment corrections.
4. Build the LRC14 odd-wall ledger: blocked unit twists plus owner/carry/fiber
   data, not scalar blocked denominators alone.
5. Look for a compact residue/tournament quotient that detects Bernoulli
   denominator drift without expanding the whole degree-`p` polynomial.
