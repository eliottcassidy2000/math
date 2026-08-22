# Keller level four: a slice kills multiplicity without expanding the norm

**Status: promoted as THM-3504 (PROVED + VERIFIED-EXACT).**

## Inheritance pass

- **Closest proved mechanism.** THM-3495 turns an irreducible source divisor
  under a finite cubic map into `uP^e`, then uses a specialization to force
  `e=1`.  THM-3498 supplies the next polynomial `G=L^43N(J)` and removes the
  old `L` denominator.
- **Canonical hostile.** A modular specialization can be squarefree while
  saying nothing about a genuinely arbitrary multivariate factorization.
  Here it is lawful only after finite-image geometry reduces the global
  possibilities to one prime power.
- **Corrected near miss.** The degree-81 eliminant gate proves generic block
  separability, not squarefreeness or irreducibility of `G`.
- **Least-used sidecar.** A valuation at the *other* boundary, `A=infinity`
  on `b=c=1`, gives a sharp anti-aliasing degree bound for the norm slice.

The session used the META-PATTERNS cards “Refine and saturate before
transporting a factor or shadow” and “Use redundant paths as detectors.”

## Concept board

| Object | Question | Preserved invariant | Lost coordinate | Cheapest test |
|---|---|---|---|---|
| `G=L^43N(J)` | Is it a new prime? | finite-norm zero support | global coefficient expansion | squarefree lawful slice |
| `V(J)->A^3` | What is image degree? | irreducibility under finite image | multiplicity `e` | proper-power hostile |
| slice `b=c=1` | Can values determine a polynomial? | pole order at infinity | other target directions | prove degree `<p` first |
| old factors `L,H,J` | Is the image new? | specialization of divisibility | factors that collapse on slice | three modular gcds |
| fourth eliminant | Is the discriminant recursion lawful? | degree/squarefreeness | image equation | inherit THM-3498 |

## The structural reframe

The apparent task was “factor the enormous global numerator of `N(J)`.”  The
finite-map geometry makes that unnecessary.  Over `U=A^3\V(L)`, `F` is finite
etale of degree three.  Since `V(J)` is irreducible and meets `F^-1(U)`, its
finite image has one irreducible hypersurface closure `V(P)`.  Norm support
and THM-3498's `gcd(G,L)=1` force

```text
G=uP^e,        1<=e<=3.
```

Thus a specialization need not prove multivariate irreducibility.  It only
has to distinguish exponent one from exponents two and three.  Squarefreeness
does exactly that.

The nonempty-open witness is explicit.  THM-3495 has `q=(3,-1,0)` with
`H(q)=0`; writing `p=F(q)` gives

```text
p=(10,-46,33),                 L(p)=-504,
F(p)=(-1854753363,121225664,-19180),
L(F(p))=-69753247104.
```

The first nonzero `L` value makes `J(p)=0` through `N(H)=J/(2^35L^7)`;
the second puts `p` in the finite locus needed for the next image.

## Why 543 values are a proof rather than a fit

On `b=c=1`, put the free target coordinate equal to `A`.  Then

```text
L=A(27A-2),       S=27A-1,       E=Lw^3+w-2.
```

At infinity all three roots have `t=1/A` valuation `2/3`.  The reduced
inverse coordinates satisfy

```text
v(q_x)=2/3,       v(q_y)>=-2/3,       v(q_z)>=-2.
```

The frozen support of `J` has

```text
max(-i+j+3k)=228.
```

So `J(q)` has pole order at most `152` per sheet, its norm at most `456`,
and `L^43` adds `86`.  Therefore `deg G(A,1,1)<=542`.  The prime `1009` is
larger than the bound.  This is the indispensable sidecar: without it,
values on `F_1009` would determine only a polynomial function modulo
`A^1009-A`.

## Exact outcome

Using `543` regular interpolation points and `12` unused regular points gives

```text
deg G_1009=542,
gcd(G_1009,G_1009')=1,
factor degrees=1,2,2,4,12,21,500,
gcd degrees against L,H,J = 0,0,0.
```

The coefficient-ledger hash is

```text
47fba77866ee50d00fcae28b834e8a0b0c18a4cf52c2cd1b9c05155410c91d00.
```

The reducible factorization is not a defect: reduction of an absolutely
irreducible characteristic-zero polynomial may split.  What matters is that
every exponent is one.  Therefore `e=1`, which returns to the geometric route
and proves `G` absolutely irreducible with generic image degree one.

## Redundant paths and hostiles

The exact evaluator has two aggregation plans for all 66,146 terms of `J`:
one sums along `y` in 4,160 fixed `(x,z)` groups, and the other sums along `z`
in 5,657 fixed `(x,y)` groups.  Held-out values agree between them.  The norm
is independently computed by a closed cubic formula and by the determinant
of the literal `3x3` multiplication matrix.

Three hostile controls are load-bearing:

1. using only `542` nodes fails at the omitted point, detecting exact degree
   `542`;
2. squaring the slice is detected by a derivative gcd of degree `542`;
3. multiplying by old `H` is detected by an old-factor gcd of degree `14`.

The first attempted run also served as a convention hostile: it exposed that
the standard closed norm formula used `w^3+pw+q=0`, while the evaluator stored
`w^3=pw+q`.  The direct determinant caught the sign translation before any
interpolation result was accepted.

## Consequence and boundary

THM-3504 now gives

```text
closure(F(V(J)))=V(G),
S_(F^4)=V(LHJ G),
[Delta_4]=[2G],
```

with four distinct prime nonproperness components.  This extends the fixed
component-count pattern through exponent four.

The computation does not print `G`, determine its primitive integral scalar,
or expose its global multidegree and singularities.  Depth five requires a
new norm of `G`, new boundary valuations, a degree-243 separability witness,
and another multiplicity/distinctness audit.  No all-level law or new
Jacobian-conjecture conclusion is inferred.

