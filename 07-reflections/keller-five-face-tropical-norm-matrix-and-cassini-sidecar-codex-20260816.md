# The fixed Keller norm tower needs five Newton faces, not one exponent

**Current-status update:** THM-3522 proves renewal propagation at every
polynomial rung, and THM-3527 now realizes the fixed packet through `R8`.
Polynomiality of the next rung and the all-level law remain open.  The text
below is a historical pre-3522 analysis; use the current theorems and the
Pell-57 sidecar for present status.

## Inheritance pass

- Closest proved mechanism: THM-3495's exact inverse chart and
  `L^7N(H)=J/2^35`, followed by THM-3498's old-`L` valuation and localization.
- Canonical hostile: MISTAKE-413's warning that constant units and nonmonic
  resultant powers cannot be discarded.
- Corrected near miss: the visible top face and pole exponent do not determine
  the next top face.
- Least-used sidecar: the opposite toric end of the inverse cubic, where
  `(Bw)^3-3Bw-2=(Bw-2)(Bw+1)^2` transports complete Newton faces.
- Incoming strengthening: THM-3504 proves `G=L^43N(J)` is the fourth image
  prime and `S_(F^4)=V(LHJ G)`.  THM-3506 starts strictly after that result.

## Concept board

1. old-boundary pole `e`;
2. exceptional-line power `m` in `x^e(3xz-2y)^m`;
3. the quadratic-edge weight `beta=i-j-2k`;
4. the Chebyshev-end weight `gamma=i-j-5k`;
5. finite-sheet units;
6. composition degree `3^n`;
7. the Cassini determinant of the projectivized face state.

The decisive change was to stop treating item 1 as a closed state.  Even the
pair `(e,m)` is propagated only after three hidden faces are supplied.

## What the inverse chart actually says

For `(a,b,c)=(A/t,B,Ct)` as `t->0`, the inverse cubic has one linear root
and two quadratic roots.  The linear root reads `min(i-k)`; the conjugate
pair reads `min(i-j-2k)`.  Thus

```text
e_next=e-lambda_min-beta_min,
m_next=(the y exponent on the opposite lambda face).
```

For the fixed packet these hidden faces are

```text
lambda_min=-e,
beta_min=-5e+2m,
in_min-lambda(P) ~ y^(3e-2m)z^e,
```

so the transform is not `e->6e+1`, `m->2e+1`, but

```text
(e,m) |-> (7e-2m,3e-2m).
```

The exceptional constants are geometric.  On the quadratic edge,

```text
y^2+27z=-9*79,
y^2+108z=-9*313,
```

which explains the `79` and `313` already visible in the top coefficient of
`J`.

At the opposite toric end, the limiting Chebyshev labels `2,-1,-1` read the
`z`-top and `gamma` faces.  The latter has the closed form

```text
z^e(27x^2z+y^3)^(e-2m/3).
```

Substitution at the simple and double labels produces respectively
`y^2+108z` and two copies of `y^2+27z`.  This is why the beta face has the
asymmetric exponent ratio `1:2`; it is branch multiplicity, not curve
fitting.

## Exact packet and the first hostile point

The five exact packets are

```text
L: (e,m)=(1,0),
H: (e,m)=(7,3),
J: (e,m)=(43,15).
```

Their beta faces factor as

```text
L: y^3z,
H: y^15z^4(y^2+27z)^2(y^2+108z),
J: y^99z^28(y^2+27z)^10(y^2+108z)^5.
```

The exact transform therefore gives

```text
G=L^43N(J): (e,m)=(271,99),
```

refuting `(259,87)`.  It also derives, without a global expansion of `G`,

```text
in_min-lambda(G) ~ y^615z^271,
in_min-beta(G) ~ y^615z^172(y^2+27z)^66(y^2+108z)^33.
```

Three good-prime evaluations prove the finite old-`L` sheet is a unit.
Hence `v_L(N(G))=-271` and `R_5=L^271N(G)` is polynomial and coprime to `L`.
The two transported faces of `G` already imply the exact exposed face

```text
in(R_5) ~ x^1699(3xz-2y)^615.
```

No finite-sheet test has yet promoted `1699` to the exact pole exponent of
`N(R_5)`.

## Conditional induction and the precise missing closure

The transform proves three output faces from five input faces.  To iterate,
the output must renew two further faces:

```text
z-top: x^(2e-4m/3)z^(2e-2m/3),
gamma-bottom: z^e(27x^2z+y^3)^(e-2m/3).
```

These are exact for `L,H,J`.  THM-3513 subsequently proves them for `G` by
the hybrid weights `gamma-k` and `gamma-3k`, making the next full matrix step
lawful for this fixed input.  The corresponding faces of `R_5` and later
rungs are still unproved; a uniform preservation proof would be needed for
the all-level induction.

The cheapest next closure tests are therefore not a global construction of
`G`:

1. obtain `max z` and `min gamma` on enough exact specializations to rule out
   lower competitors;
2. evaluate the corresponding initial resultants at the labels `2,-1`;
3. separately test the finite old-`L` value of `R_5` before calling `1699` a
   pole exponent.

## Cassini, reduced fractions, and primitive triples

Conditionally let `v_(n+1)=Mv_n` with

```text
M=[[7,-2],[3,-2]],       v_0=(1,0).
```

Then

```text
det(v_n,v_(n+1))=3(-8)^n,
e_n=1 mod 6,
m_n=3 mod 6 for n>=1,
gcd(e_n,m_n)=1.
```

The gcd proof is sharp: it divides the preceding Cassini determinant
`3*8^(n-1)`, while the mod-6 classes exclude both primes `2` and `3`.
Therefore `m_n/e_n` is reduced, but consecutive cross-determinant has
magnitude `3*8^n`, so these fractions are not Farey neighbors and do not form
a Stern--Brocot edge path.

The projective map

```text
r |-> (3-2r)/(7-2r)
```

reverses orientation and contracts `[0,1]` to `[1/5,3/7]`; its ratios
alternate around `(9-sqrt(57))/4`.  Since the parameters are coprime odd,
they give primitive Pythagorean triples after parity division:

```text
(20,21,29),
(812,645,1037),
(31820,26829,41621), ...
```

This is a genuine Fibonacci-style Cassini/projective sidecar.  It does not
turn the norm tower into the Fibonacci sequence, a Farey path, or an
unconditional ternary-tree classification.

## Monoid comparison

The fibre degree is multiplicative under composition: the fixed orbit has
degree `3^n`.  The Newton pair is not that grading.  It is a signed
projective boundary state, updated by a norm and vulnerable to face and
finite-sheet cancellation.  The useful common grammar is only this:

```text
composition transports a multiplicative degree;
the function-field norm transports a multi-face Newton packet.
```

Conflating them would turn a fixed-map calculation into an unsupported
classification claim.  THM-3506 intentionally does not do that.
