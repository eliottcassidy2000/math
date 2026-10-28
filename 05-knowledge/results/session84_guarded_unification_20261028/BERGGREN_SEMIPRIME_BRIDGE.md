# Berggren addresses for square-sigma semiprimes, and the exact 2–3–11 root

Date: 2026-10-28. **PROVED** elementary guarded bijections and operation laws;
**FINITE-EXACT** independent enumeration through source 1,000,000. Classical
Euclid/Berggren enumeration is **CITED**. No literature-priority claim and no
claim of infinitude for any prime-producing quadratic is made.

## 1. What is inherited, and what is being connected

Write `s(N)=sigma(N)−N`. For distinct primes `p<q`,

```
N=pq,  t=s(N)=p+q+1,  sigma(N)=(p+1)(q+1).
```

The source is a square-sum aliquot edge exactly when the last product is a
square. The repeated-prime diagonal is forbidden: `s(p²)=p+1`, not `2p+1`.
These are Session83's class-II sources and centered-prime-pair identity.
The source shape is not preserved by its aliquot successor.

The precise prior repository mechanisms, read from clean source commit
`ae7c073e864d05be74fed85ba0a2bc164e026a6d` using `git show`, are:

* `01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md`:
  ordered Euclid/Gaussian lift, with content retained.
* `01-canon/theorems/THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md`:
  current U/A/D convention and classical unique tree addresses. Its
  non-backtracking and branch-character results are not required below.
* `05-knowledge/results/twin_pythagorean_124_20261026/tree/RESULT.md`:
  the D-ray, distinction between prime-to-prime tree edges and prime-filtered
  nodes, and the warning against calling a subset bijection a graph conjugacy.
* Session83 `FIBER_SHAPE_CENSUS.md`, §§3–4, and
  `CARRIER_AND_144_BRIDGE.md`: the centered Pythagorean equation and coprime3
  transport, respectively. The present result makes the former's tree
  address, scale, leg choice, and child-primality guards explicit.

Classical sources: Alperin's [*The Modular Tree of Pythagorus*](https://arxiv.org/abs/math/0010281)
gives the ternary-tree structure through the modular group; Romik's
[*The dynamics of Pythagorean triples*](https://arxiv.org/abs/math/0406512)
relates finite ternary addresses to a modified Euclidean algorithm. These
sources supply the tree background, not the new guarded aliquot formulas.

The live concept board is **distinct semiprimes / square sigma / primitive
tree address / nonprimitive scale / prime gap / coprime multiplication**.
The canonical hostile is a repeated prime or a composite candidate after a
tree step. The least-used sidecar is **which Pythagorean leg is the half-gap**.

## 2. Complete, injective map for odd sources

Use a positive primitive triple in parity orientation

```
(a,b,c)=Phi(m,n)=(m²−n²,2mn,m²+n²),
m>n>=1, gcd(m,n)=1, m+n odd.
```

There are exactly two charts, each with both displayed endpoints required
to be prime:

| Half-gap leg | Allowed scale | Prime endpoints | Square root of sigma |
|---|---|---|---|
| odd leg a | d>=1 | `p=d(c−a)−1=2dn²−1`, `q=d(c+a)−1=2dm²−1` | `db` |
| even leg b | even d>=2 | `p=d(c−b)−1=d(m−n)²−1`, `q=d(c+b)−1=d(m+n)²−1` | `da` |

In both charts, the aliquot edge is

```
pq -> 2dc−1,                  q−p=2d·(chosen half-gap leg).
```

**Theorem.** The guarded decorated data `(primitive triple,d,half-gap leg)`
above are in bijection with **all odd distinct-semiprime square-sigma
sources**. Their Berggren word can replace the primitive triple without
changing the bijection.

**Proof.** Forward, if `r` denotes the selected primitive leg and `z0` the
other, then

```
(p+1)(q+1)=d²(c²−r²)=d²z0².
```

The prime guards imply that `p,q` are genuine distinct primes, because
`q−p=2dr>0`. In the odd-leg chart `c±a` are even, so both primes are odd.
In the even-leg chart `c±b` are odd; even d makes both primes odd. Odd d
would make both endpoints even, so could never give distinct primes.

Conversely take odd primes `p<q` with square product and set

```
r=(q−p)/2, z=sqrt((p+1)(q+1)), c1=(p+q+2)/2.
```

Then `r²+z²=c1²`, with both legs positive. Let `d=gcd(r,z,c1)` and divide
by d. Exactly one reduced leg is odd; choose the chart according to whether
that leg is `r/d` or `z/d`. This gives a uniquely oriented primitive triple.
The original endpoints are `d(c−r0)−1,d(c+r0)−1`, and the parity argument
above supplies the even-d guard when needed. Reduction and leg parity are
unique; unique prime factorization recovers `p,q` from `N=pq`. Thus the map
is both surjective and injective. QED.

Primitivity alone is not enough: the root `(3,4,5)` at d1 in the odd-gap
chart gives `(p,q)=(1,7)`, not primes. Scaling is genuinely necessary:
`115=5·23` comes from `3(3,4,5)`, half-gap9 and sigma root12. Conversely
`119=7·17` comes from the primitive `(5,12,13)` without scaling.

### A factor-coordinate alternative and two exact obstructions

Equivalently let `g=gcd(p+1,q+1)`. The reduced factors are coprime squares:

```
p+1=g u², q+1=g v², gcd(u,v)=1, 1<=u<v.
```

The integer g is even; it need not be squarefree. Opposite parity of `(u,v)`
is the odd-gap chart, with `g=2d,(u,v)=(n,m)`. Both odd is the even-gap
chart, with `g=d,(u,v)=(m−n,m+n)`. These are exact alternative coordinates.

No such source contains the prime3. Indeed if p3, then
`4(q+1)=z²` forces `q+1=h²`; a prime `q=h²−1=(h−1)(h+1)` forces h2 and q3,
contradicting distinctness. Also g cannot be a square: then both p and q
would be one less than a square, and the only prime of that form is3.
Consequently the odd-gap chart cannot have `2d` square, and the even-gap
chart cannot have d square. These are necessary filters, not replacements
for primality.

### Fixed gaps are finite, so this is not the older twin-prime image

For a fixed odd-prime gap `q−p=2r`, the centered identity is equivalent to

```
(c1−z)(c1+z)=r².
```

Thus all candidates are obtained from the finitely many positive factor
pairs `A<B`, `AB=r²`, `A+B` even, by

```
c1=(A+B)/2, z=(B−A)/2,
p=c1−r−1, q=c1+r−1,
```

retaining only odd distinct primes. This is a complete finite algorithm for
each fixed gap, with `p<=(r²+1)/2−r−1`. Consequently bounded gaps give only
finitely many square-sigma semiprimes. In particular no twin-prime pair
works: at r1 there is no strict factor pair. The smallest possible odd-prime
gap is10, attained by7,17: r2 has no positive triangle, r3 gives only
`(p,q)=(1,7)`, and r4 only `(0,8)`. This result concerns the prime factors
of the semiprime, not whether its aliquot endpoints form a small-gap edge.
The earlier twin-prime/Berggren dictionary concerns prime inradius/exradius
or adjacent prime-leg nodes and therefore has a different prime predicate.

## 3. Every even source is three times one distinguished ray

The even case has `p=2`. Square sigma means

```
3(q+1)=z²  iff  q=12k²−1, z=6k,
```

with `k>=1` and q prime. Divisibility first gives z3u and q3u²−1;
q odd then makes u even, giving u2k. Thus the complete even family is

```
N=24k²−2 -> t=12k²+2, sigma(N)=36k².
```

The unhalved centered triangle is exactly

```
(q−2, 2z, q+4)
 = (12k²−3,12k,12k²+3)
 = 3 Phi(2k,1).
```

Its primitive content is exactly3, because `Phi(2k,1)` is primitive.
With the branch convention

```
U(m,n)=(2m−n,m), A(m,n)=(2m+n,m), D(m,n)=(m+2n,n),
```

we have `Phi(2k,1)=D^(k−1)(3,4,5)`. Every even semiprime square-sigma source
therefore lies on this precise scale3 D-ray image, and every prime-guarded
point of the image is such a source. The first four q are11,47,107,191;
the next is299=13·23, so even this ray does not preserve primality.

### The 2–3–11 anchor and the lifted target exclusion

At k1 the primitive root gives `q=11`, hence `22->14`. Since `3` does not
divide `22`, multiplication by3 gives `66->78` and total144. The number3
thus appears in two proved but different roles: primitive content of the
unhalved centered triangle, and a coprime multiplier with `sigma(3)=4`.
Neither role says the three integers2,3,11 themselves are a Pythagorean triple.

Generally q is `−1 mod3`, so the valid coprime lift is

```
6q -> 6q+12 = 72k²+6,        sigma(6q)=144k².
```

**Every one of these lifted targets has no class-II predecessor at all.**
For an even target t a distinct-semiprime predecessor must be `2(t−3)`,
with `t−3` an odd prime. Here `t−3=3(24k²+1)>3` is composite. This is an
all-fiber exclusion, stronger than saying the displayed source is classIII.
It does not imply a singleton fiber: at k2, t294 receives both `282=6·47`
and `293²=85849`, the latter a class-I prime square with aliquot sum294.
The k1 singleton `F78={66}` is the separately certified Session83 fact.

The other sigma144 source70=2·5·7 has a different signature cancellation:
`sigma(70)=3·6·8=144`; it is not a semiprime source and is not generated
by the map `2q -> 6q`. Sharing144 does not identify the mechanisms.

## 4. Exact child formulas: one prime survives, one must be retested

For any positive p<q with `(p+1)(q+1)=z²`, define candidate children

```
F−(p,q)=(q, p+4q+4−4z),    z−=2(q+1)−z;
F+(p,q)=(q, p+4q+4+4z),    z+=2(q+1)+z;
F0(p,q)=(p, 4p+q+4+4z),    z0=2(p+1)+z.
```

Each new pair `(P,Q)` satisfies `P<Q` and `(P+1)(Q+1)=Z²` with the indicated
positive root. For example, putting `x=p+1,y=q+1`, the first identity is
`y(x+4y−4sqrt(xy))=(2y−sqrt(xy))²`. Its endpoint order follows from
`x+3y−4sqrt(xy)=(sqrt(y)−sqrt(x))(3sqrt(y)−sqrt(x))>0`.
The other two cases follow by the plus signs.

For odd-source data, these are the actual Berggren children:

| Chart | F− | F+ | F0 |
|---|---|---|---|
| odd half-gap a | U | A | D |
| even half-gap b | D | A | U |

This follows either by substituting the Euclid branch formulas or by the
triple matrices. The scale and chosen-leg parity are preserved. The U/D
swap is the coordinate effect of interchanging which leg is the half-gap.
Thus one prime endpoint survives each step, but **only the algebraic
square-product identity, not the arithmetic divisor-sum interpretation,
survives without testing the new endpoint**.

All three branches can fail, in both charts:

| Chart | Branch | Parameters `(m,n;d)` | Prime pair -> child candidate |
|---|---|---|---|
| a | U | `(4,3;1)` | `(17,31)->(31,49)` |
| a | A | `(2,1;6)` | `(11,47)->(47,299)` |
| a | D | `(2,1;3)` | `(5,23)->(5,95)` |
| b | U | `(2,1;12)` | `(11,107)->(11,299)` |
| b | A | `(2,1;8)` | `(7,71)->(71,391)` |
| b | D | `(2,1;12)` | `(11,107)->(107,299)` |

Here `49=7²`, `299=13·23`, `95=5·19`, and `391=17·23`.
The filtered subset is not ancestor-closed either: the root in chart a,
d1 gives the invalid pair `(1,7)`, while its U-child gives `(7,17)` and
the valid edge `119->25`, total144. Therefore **a failed prime guard cannot
be used to prune the entire descendant subtree**. Keeping one endpoint is
not keeping both, and a later step can discard the previously bad endpoint.
For p2, F− and F+ give an even new endpoint greater than2, so only F0
can remain a prime pair. Its new q is `12(k+1)²−1`: the same D-ray advance.
These are candidate-pair tree operations, **not aliquot iterates**.

## 5. Nonprimitive scaling is an affine transport, not multiplication of N

If `(p,q)` is one odd admissible pair, multiplying its centered triangle by
a positive integer h gives candidate endpoints

```
p_h=h(p+1)−1, q_h=h(q+1)−1.
```

If and only if both candidate endpoints are prime, the original semiprime
interpretation survives. For those admissible scales,

```
t_h=h(t+1)−1,
N_h=h² sigma(N)−h(t+1)+1,
sigma(N_h)=h² sigma(N).
```

The candidate-pair maps compose by multiplication of h, but the prime
subset is not closed under them. Example: `115->29` at h2 becomes
`517=11·47 ->59` with square total576. At h4 the candidates are23 and95;
their product2185 has actual sigma2880, not the formal square2304. The
first lost condition is primality of95, not the Pythagorean identity.

Scale and leg choice cannot be forgotten. At the same primitive root and
same scale d6, chart a gives `517=11·47 ->59`, total576, whereas chart b
gives `265=5·53 ->59`, total324. The centered hypotenuse30, target59, word,
and scale agree, but the prime gaps36 and48 and square roots24 and18 differ.
The missing sidecar is the selected half-gap leg. Conversely a semiprime
source retains all these data through its unique prime factorization.

## 6. Verification, scope, and stopping boundary

Run the dependency-free checker from this directory:

```
python3 check_berggren_bridge.py > results.json
python3 -O check_berggren_bridge.py > /tmp/berggren_optimized.json
cmp results.json /tmp/berggren_optimized.json
```

The normal and optimized executions byte-match. The explicit source universe
is every distinct-prime product `pq<=1,000,000`: **209,867** prime pairs.
Exactly **901** have square sigma, partitioned as **540** odd-gap-chart,
**285** even-gap-chart, and **76** even-source D-ray images. Independent
direct divisor sums and trial primality check every accepted source.

An independent Berggren forward traversal and Euclid-pair enumeration agree
on all **15,919** primitive triples with `c<=100,003`. This bound suffices
for the odd-source comparison: p>=5, so
`dc=(p+q+2)/2 <= (floor(1,000,000/5)+7)/2`. For p<=sqrt(B), the function
`p+B/p` decreases, justifying the upper bound at p5. The checker compares
the entire generated source set with the prime-pair census, checks inverse
addresses and injectivity, all six child failures, candidate-child formulas,
the first D-ray prime failure, scale loss, and the lifted-fiber exclusion.
An independent factor-pair versus bounded-p scan verifies all fixed-gap
candidate sets for `1<=r<=100`, including the least admissible odd gap10.

The exact proofs are global; these counts and minimum witnesses are
finite-exact. No density, prime infinitude, full aliquot-orbit conjugacy,
global graph equivalence, or new tournament structure is asserted.
