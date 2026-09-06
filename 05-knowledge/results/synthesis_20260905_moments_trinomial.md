# Trinomial returns: exact tunability, hidden carries, and unbounded resonance

**Status:** PROVED and independently paper-audited by the root and wildcard
lanes for Sections 2--5; FINITE-EXACT for the two-rung coprimality census in
Section 6. Current-surface integration is separate. No general sharp
first-return bound or new Gaussian-moment nullcone theorem is claimed.

**Session:** synthesis-20260905, moments lane. Source tree at inheritance:
`3eb2b8a66e56`; concurrent later LRC/JC work does not change these inputs.

## 1. Inheritance and concept board

The closest proved mechanism is
[THM-2111, effective compound-root bound](../../01-canon/theorems/THM-2111-effective-compound-root-bound-for-one-variable-constant-terms.md):
the first nonzero constant term is exactly the order of a compound determinant,
and is at most `binom(M+N,min(M,N))`. Its proof was read, including the exact
simple-pole boundary. The sharper `M+N` statement remains open in the repo.
[THM-2260, positive finite-fibre capacity](../../01-canon/theorems/THM-2260-positive-finite-fibre-capacity-and-thin-predicate-boundaries.md)
already rules out reducing that binomial count by a generic image-size argument.

The canonical hostile is
[THM-2070, horizontal Wick embedding and dihedral cancellation](../../01-canon/theorems/THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md):
support return is not coefficient noncancellation. The corrected near miss is
[THM-1710, cyclotomic refutation](../../01-canon/theorems/THM-1710-tnc-cyclotomic-refuted-resultant-replaces.md):
tuned cancellation parameters need not be roots of unity. The recovered sidecar
is [THM-1680, trinomial gcd decision](../../01-canon/theorems/THM-1680-tnc-trinomial-gcd-decision.md),
especially its still useful two-rung question. Another proved neighbour is
[THM-2639, free equal-mass two-rung certificate](../../01-canon/theorems/THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate.md),
whose full-semigroup freeness hypothesis cannot be inferred from one slice.

The live board was: compound contact order; numerical-semigroup gaps;
coefficient-weighted return ideals; residue carries between return levels;
horizontal Wick transport. The inherited assumption to challenge was the
small-relation interpretation of tunability in
[THM-1755, angular tunable dichotomy](../../01-canon/theorems/THM-1755-tnc-angular-uniform-tunable-dichotomy.md).
The present proof closes that file's request for a precise three-charge
tunability criterion, with a different criterion and an explicit refutation
of every coefficient-independent bound on the shortest signed relation.

Relevant protocol cards: search the statement before the method; recover the
missing coordinate; and re-evaluate a certificate after a fibre-changing
operation. Exact statement, theorem-slug, and MISTAKES searches found no prior
correction of THM-1755's single `c=1` family or bounded-height assertion.

## 2. Exact classification of the first support return

Let

```text
f(u)=alpha*u^(-a)+beta*u^b+gamma*u^c,
alpha*beta*gamma != 0,
1<=b<c, a>=1, gcd(a,b,c)=1.
```

Reflection covers the case of two negative charges. Divide out common charge
content before applying the formulas; common charge scaling preserves every
constant term. A neutral monomial is detected at the first moment and requires
none of the following analysis.

Set

```text
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g,
S={A*y+B*z:y,z>=0 integers}.
```

Then `1<=A<B`, `gcd(A,B)=1`, and `gcd(a,g)=1`.
Let `B_m` denote the set of multiplicity vectors `(x,y,z)>=0` satisfying

```text
x+y+z=m,  -a*x+b*y+c*z=0.
```

These are **unordered multiplicity channels**, not ordered words. Their
multinomial weights are retained separately in the scalar moment.

**Theorem 1 (classification).** Define

```text
k0=min{k>=1:a*k in S}.
```

Then the first support return is `m0=g*k0`. Moreover,

```text
|B_m0|>=2  iff  a-A*B in S.                       (1)
```

Whenever (1) holds,

```text
k0=1, m0=g, g>B, 2*m0<=a+c.                     (2)
```

Thus the first return can be tuned to zero by nonzero complex coefficients
exactly on the semigroup locus (1). This is an all-exponent classification.

### Proof

The two balance equations are equivalent to

```text
(a+b)*y+(a+c)*z=a*m.
```

Since `gcd(a,g)=1`, every return has `m=g*k`. The equation then becomes
`A*y+B*z=a*k`. Conversely any such nonnegative pair gives the integer
`x=g*k-y-z`, and the original identity implies `a*x=b*y+c*z>0`.
Hence every such pair is a valid channel. Existence is finite, for example
from either positive-negative binomial return.

Here is the needed semigroup fact, proved directly. For every integer `n`,
write uniquely

```text
n=A*y+B*z, 0<=z<A, y integer.
```

Then `n in S` iff `y>=0`, and all representations are

```text
(y-B*j,z+A*j), 0<=j<=floor(y/B).
```

In particular `n` has at least two representations iff `n-AB in S`.
Put `F=AB-A-B`. The canonical representation of `F-n` is

```text
F-n=A*(-y-1)+B*(A-1-z).
```

Consequently exactly one of `n` and `F-n` belongs to `S`. This proves the
usual two-generator gap symmetry without an external input, including `A=1`.
In particular `F` is not in `S`.

Suppose the first returned number `a*k0` has at least two representations,
but `k0>=2`. Then `a*k0-AB in S`, while minimality makes both `a` and
`a*(k0-1)` gaps. Gap symmetry puts both `F-a` and `F-a*(k0-1)` in `S`.
Adding these three elements would give

```text
2*F-AB=F-A-B in S.
```

Adding `A+B` would then represent `F`, a contradiction. Therefore `k0=1`,
and the multiple-representation test is precisely `a-AB in S`. The reverse
implication is immediate from the representation formula.

Finally, collision implies `a>=AB`. But `g*A=a+b>a`, so `g>B`. Since
`B>=2`, the width is `a+c=g*B>=2*g`, proving (2).

To justify the word *tunable*, use nonzero scalar and variable dilations to
normalize the extreme coefficients to one. The first constant term becomes
a monomial in the remaining nonzero parameter times a polynomial with a
nonzero constant term. It is nonconstant exactly when there are at least two
channels, and hence has a nonzero complex root. This asserts existence of a
coefficient choice, not cancellation for every choice.

### Boundaries

- If the first channel is singleton, its nonzero multinomial times its torus
  coefficient monomial already detects `f`.
- The first support return is at most the smaller positive-negative pair
  return, and hence at most `a+c`.
- In (2), equality `2*m0=a+c` means `B=2`, hence `A=1`.
- Gaps and multiplicities do not decide whether a selected coefficient point
  cancels. That remains a polynomial-ideal question.

## 3. The all-level carry profile and the exact freeness boundary

Assume the first return is collided. Let `(y0,z0)` be the canonical
representation of `a`, so `0<=z0<A`, and write

```text
y0=B*h+r, 0<=r<B, h>=1.
```

**Theorem 2 (channel profile).** At every `k>=1`,

```text
|B_(k*g)|=k*h+floor(k*r/B)+floor(k*z0/A)+1.        (3)
```

The distinct sums of `k` channels from `B_g` number exactly `k*h+1`.
Therefore the number of channels that cannot be assembled from that first
slice is exactly

```text
floor(k*r/B)+floor(k*z0/A).                       (4)
```

All channels at every level are generated by the first slice iff

```text
r=z0=0  iff  a is a positive multiple of A*B.     (5)
```

If `|B_g|=2`, the all-level free-two-ray case is exactly `a=A*B`.

**Proof.** The canonical representation of `k*a` is

```text
z_k=k*z0-A*floor(k*z0/A),
y_k=k*y0+B*floor(k*z0/A).
```

Its number of representations is `floor(y_k/B)+1`, giving (3).
The first-level channels form an integer segment indexed by `0,...,h`.
Its `k`-fold sum is the segment indexed by `0,...,k*h`, proving the count
and (4). A nonzero residue produces a positive floor at some finite `k`;
zero residues produce none. This proves (5). In the two-channel case `h=1`,
zero residues say `a=AB`; the two endpoint rays have no additive relations
and generate the full semigroup. Conversely freeness on two first-level
rays forces (4) to vanish at every `k`.

At the second level, the lost-channel count is especially small but real:

```text
|B_(2*g)|=2*|B_g|-1+epsilon_y+epsilon_z,
epsilon_y=floor(2*r/B), epsilon_z=floor(2*z0/A),
epsilon_y,epsilon_z in {0,1}.                    (6)
```

The smallest witness found by the ordered height scan with two first-level
channels and both carries is `(-13,1,8)`. Here `(g,A,B)=(7,2,3)` and
`(y0,z0)=(5,1)`. The first channels are

```text
(1,5,1), (2,2,3).
```

Their pairwise sums account for three channels at mass fourteen, but two
additional channels are

```text
(1,13,0), (5,1,8).
```

Thus two first-level channels do not justify a three-term second-level
formula. This identifies the missing sidecar before applying THM-2639.

## 4. The old single-family and bounded-resonance claims are false

The `c=1` notation in THM-1755 corresponds to the middle charge `b=1` here.
Its statement that tunable triples are exactly `(-a,1,a+2)` is refuted by

```text
charges=(-3,1,9), m0=4,
f=u^(-3)+T*u+u^9,
CT(f^4)=4*(T^3+1),
CT(f^8)=28*(T^6+10*T^3+1).
```

On `T^3=-1`, the eighth constant term is `-224`. The first failed
implication was extrapolating the small height census to a single family.
The correct survivor is the tunable/singleton dichotomy, now characterized
by (1), not a bounded list of small resonance patterns.

**Theorem 3 (no bounded signed-relation height).** Choose any coprime
`1<=A<B`, and any `g>2B` coprime to `AB`. Put

```text
a=A*B, b=A*(g-B), c=B*(g-A).
```

These are primitive positive data with `b<c`. The support `(-a,b,c)` is
tunable, its first return is `g` with exactly two channels, and the shortest
nonzero signed integer charge relation has exact `l1` norm `2B`.

**Proof.** The primitive assertion follows by checking a prime dividing
`AB`: coprimality with `g` makes it miss at least one of `b,c`. Also
`a+b=gA`, `a+c=gB`, and `a=AB`, so Theorems 1--2 apply. Now let an arbitrary
integer vector `(R,S,T)`, permitting zero entries and either sign, satisfy

```text
-AB*R+A*(g-B)*S+B*(g-A)*T=0.
```

Rearranging gives

```text
g*(A*S+B*T)=AB*(R+S+T).
```

Since `gcd(g,AB)=1`, `g` divides `R+S+T`. If the vector had `l1<2B<g`,
its sum would have to be zero. Then `A*S+B*T=0`, and coprimality yields

```text
(R,S,T)=j*(B-A,-B,A),
```

whose norm is `2B*|j|`. No nonzero vector can therefore have norm below
`2B`, and the displayed vector with `j=1` attains it.

For an explicit unbounded family, take `A=n`, `B=n+1`, and `g=AB+1`,
with `n>=2`. The first violation of the proposed bound six is

```text
charges=(-12,27,40), first mass=13,
shortest relation=(1,-4,3), l1=8.
```

These examples already belong to THM-2639's proved free equal-mass class.
Their role is to refute bounded resonance as a necessary explanation of
tunability, not to challenge noncancellation or NC2. The general support
classification needs the arithmetic coordinates `(A,B,a,g)`, not a universal
small relation vector.

## 5. Typed transfers and the sharpened residual

**Numerical-semigroup transfer.** The source is the complete channel set of a
primitive trinomial; the target is representations of `a*k` by `A,B`; the map
is `(x,y,z)->(y,z)` at `m=g*k`, with inverse `x=g*k-y-z`. This is a bijection,
so it preserves support existence and multiplicity exactly. If one retains
only semigroup membership, it destroys channel multiplicity and all complex
coefficient phases. Restore the canonical residue pair and the weighted
polynomial `sum m!/(x!y!z!)*alpha^x*beta^y*gamma^z` before claiming
noncancellation. The decisive test is the two-extra-channel hostile above.

**THM-2639 transfer.** The source is the first-level integer segment; the
target is its generated graded subsemigroup. Addition preserves balance and
mass. The quotient that declares this subsemigroup complete loses exactly
the two floor terms in (4). It is valid at all levels iff `a` is a multiple
of `AB`; the two-generator free case requires `a=AB`. Thus (3)--(5) give
an exact entry test for the older two-rung mechanism and exhibit why its
assumptions cannot be relaxed to two channels in one row.

**Horizontal Wick transfer.** Choose `H>=c` and put

```text
P(Z,W)=alpha*Z^H*W^(H+a)
       +beta*Z^H*W^(H-b)
       +gamma*Z^H*W^(H-c).
```

Under the normalized Wick functional `L(Z^i W^j)=delta_(i,j)*i!`,

```text
L(P^m)=(H*m)!*CT(f^m).
```

This exact map preserves every zero/nonzero moment and first detection
order, not merely support. It retains the common height as a sidecar.
It does not control general multi-face Gaussian detection or move a
good-prime coefficient dependence out of THM-2022. THM-2070 supplies this
embedding generally; here it transports the new channel classification and
the exact carry obstruction.

**Effective residual.** Every singleton first return already obeys the
sharp width bound. Every collided first return has `m0=g` and `2g<=a+c`.
Thus the following stronger, precise assertion would prove the sharp width
bound for every trinomial:

```text
CT(f^g) and CT(f^(2g)) have no common torus zero
whenever a-AB belongs to <A,B>.
```

This is a sufficient route, not an equivalence: for `B>2`, a later return
could still lie within width `gB`. The next proof should retain the two
residue carries in (6), rather than assume the second row is a square of
the first. The present work does not prove this all-exponent coprimality.

## 6. Reproducible exact evidence

Run:

```text
python 04-computation/synthesis_20260905_moments_trinomial.py
python -O 04-computation/synthesis_20260905_moments_trinomial.py
```

Source:
[synthesis_20260905_moments_trinomial.py](../../04-computation/synthesis_20260905_moments_trinomial.py).
Output: [matching exact transcript](synthesis_20260905_moments_trinomial.out).
The [optimized transcript](synthesis_20260905_moments_trinomial_optimized.out)
is byte-identical. All nineteen verification sites use explicit `require`
calls which remain live under `-O`; there are no removable Python assertions.

SHA-256 of LF-normalized source and output bytes:

```text
source: 9b979a778acc1b725a34b4d1417245aab9850c45660a0e6d77db636f7b05653e
output: e21ec48936d5e52b935a68049b1fbcdc13d99f26101e0f6cc6ab52df272d78dd
```

The universe is all primitive triples `1<=a<=60`, `1<=b<c<=60`, with no
phase or coefficient specialization in the gcd calculation. There are
`88,830` supports: `86,555` singleton and `2,275` collided first returns.
The script independently computes channels by scanning the negative count
and solving the original charge equation, verifies every earlier mass is
empty, and compares against the semigroup classifier. All `11,375`
checks of (3) at `k=1,...,5` agree. The second-level excess counts are
`911` with zero carries, `1,153` with one, and `211` with two.

For each collided support, normalize extreme coefficients to one, remove
the common middle-coefficient monomial, and use its `B`th power as a single
variable. The first two return polynomials have gcd one over `Q` in all
`2,275` cases. This is FINITE-EXACT over arbitrary complex torus coefficients
on those supports; it is not a random coefficient census and not a proof
for unbounded exponents. The maximum first-row multiplicity in this box is
thirty, at `(-58,1,60)` and mass fifty-nine.

Controls include the old-family counterexample, the two-carry hostile,
four independent bounded shortest-relation searches, and scaling checks.
The paper proofs are independent of the census. No literature-priority
claim follows from a repository search. The original effective problem was
checked against the primary Erman--Smith--Varilly-Alvarado source,
[*Laurent polynomials and Eulerian numbers*](https://arxiv.org/abs/0908.2609),
Conjecture 1; its conditional Eulerian degree theorem does not supply the
missing all-trinomial two-rung proof.
