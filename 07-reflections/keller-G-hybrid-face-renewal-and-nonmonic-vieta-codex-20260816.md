# The two missing faces of `G` meet at one hybrid Newton vertex

**Status:** PROVED + VERIFIED-EXACT + INDEPENDENTLY REPLAYED as THM-3513.
This note concerns only the fixed polynomial
`G=L^43 N(J)`.  It does not supply the next finite-sheet unit, a fifth image
prime, an all-level norm induction, or a general Jacobian-conjecture claim.

## Inheritance pass

- Closest proved mechanism: THM-3506's inverse-chart five-face transform and
  exact `(271,99)` packet for the three already transported faces of `G`.
- Canonical hostile: MISTAKE-415, where an exposed pair was treated as a
  closed state before every inverse-branch face had been checked.
- Corrected near miss: the visible `z`-top monomial of `J` is not enough by
  itself; it must also be the unique extremum for the actual hybrid weights
  read by the two target degenerations.
- Least-used sidecar: the nonmonic Vieta factor.  The leading coefficient
  `D=27A^2C+B^3` of the residual cubic is exactly what reconstructs the full
  gamma-face binomial rather than one sampled monomial.

The live concept board was: `J`'s `z`-top, `J`'s gamma face, the two hybrid
weights, the two residual cubics, nonmonic norm normalization, and the next
finite sheet.  The first six close the renewal gate; the seventh remains a
separate obligation.

## 1. The hidden common endpoint

Write a monomial of `J(x,y,z)` as `x^i y^j z^k`, and retain

```text
gamma   = i-j-5k,
delta_6 = i-j-6k = gamma-k,
delta_8 = i-j-8k = gamma-3k.
```

The frozen `66,146`-term ledger has

```text
max k       = 76,
min gamma   = -314,
```

and the unique `k=76` term is

```text
c_J x^66 z^76,       c_J=2^15 3^171.
```

The exact gamma face has `34` terms and meets the `z`-top face only at this
monomial.  Therefore

```text
delta_6 >= -314-76   = -390,
delta_8 >= -314-3*76 = -542,
```

with equality in either inequality only at `x^66z^76`.  A direct scan of all
`66,146` terms independently gives the same singleton faces.  The next
weights are `-389` and `-539`, respectively, so there is no equal-weight
competitor hidden behind the endpoint argument.

This is the economical connection:

```text
source:       z-top face and gamma face of J
map:          intersection through delta_6=gamma-k and delta_8=gamma-3k
preserved:    the unique coefficient c_J
lost:         all higher hybrid layers
sidecar:      strict next-weight gaps
decisive test: direct full-ledger hybrid scan.
```

## 2. The `c`-top degeneration

Put

```text
(a,b,c)=(A,B,C s^3),       w=q/s,       s -> infinity.
```

The inverse cubic has leading equation

```text
27A^2C q^3-2=0.                                      (1)
```

Exact weighted initial forms of the inverse numerators give, on every one
of its three roots,

```text
(x,y,z) ~ (q/s, -s/q, -C s^6/q^3).                   (2)
```

Thus this degeneration reads `delta_6`.  The singleton face contributes

```text
J(x,y,z)
 ~ c_J C^76 q^-162 s^390
 = c_J (27/2)^54 A^108 C^130 s^390,                  (3)
```

where (1) makes the leading value root-independent.  Meanwhile

```text
L^43 ~ (27A^2C^2)^43 s^258.                           (4)
```

Taking the product of the three values in (3) and multiplying by (4) gives

```text
G(A,B,Cs^3)
 ~ (3^1128/2^117) A^410 C^476 s^1428.                (5)
```

Since `G` is polynomial by THM-3498, powers of `s` are exactly three times
the `c`-degrees of its monomials.  Equation (5), together with uniqueness of
the `delta_6` face, proves the complete face

```text
max deg_c(G)=476,
in_max-c(G)=(3^1128/2^117) a^410 c^476.               (6)
```

There is no `b`-dependent equal-degree term: (5) is the complete leading
coefficient in the generic variables `A,B,C`, not one specialization.

## 3. The target-gamma degeneration

Now put

```text
(a,b,c)=(At,B/t,C/t^5),       w=qt,       t -> 0,
D=27A^2C+B^3.
```

The three residual roots satisfy the nonmonic cubic

```text
Dq^3-3Bq-2=0.                                        (7)
```

Again the exact inverse numerators reduce on (7) to

```text
(x,y,z) ~ (qt, -1/(qt), -C/(q^3t^8)).                 (8)
```

This degeneration reads `delta_8`.  Every branch therefore contributes

```text
c_J C^76 q^-162 t^-542.                               (9)
```

The nonmonic leading coefficient in (7) is load-bearing.  Vieta gives

```text
product(q)=2/D,
product(q^-162)=(D/2)^162.                            (10)
```

Also `L^43~(CD)^43t^-344`.  Multiplying (9)--(10) yields

```text
G(At,B/t,C/t^5)
 ~ (3^513/2^117) C^271 D^205 t^-1970.                 (11)
```

The exponent of `t` is exactly the target weight `i-j-5k`, so (11) proves
the complete face

```text
min gamma(G)=-1970,
in_min-gamma(G)
 =(3^513/2^117)c^271(27a^2c+b^3)^205.                 (12)
```

Discarding `D` as a harmless unit would destroy the face support and repeat
the normalization pattern warned about by MISTAKE-413.  Here it is retained
before any specialization.

## 4. What closes, and what does not

Equations (6) and (12) are exactly the two missing renewal faces for the
THM-3506 state `(e,m)=(271,99)`:

```text
2e-4m/3 = 410,
2e-2m/3 = 476,
e-2m/3  = 205,
-8e+2m  = -1970.
```

Consequently this fixed `G` now has the full five-face packet.  THM-3506's
already proved transform to the exposed pair `(1699,615)` no longer depends
on an unproved renewal assertion for `G`.

The next norm still has a distinct gate.  Its two divergent old-`L` sheets
are controlled by the exposed face, but an exact unit on the finite inverse
sheet is still required before one may identify the pole exponent `1699`,
clear the denominator exactly, or call the cleared norm a fifth image
equation.  Nothing here proves that gate, irreducibility, degree `243`
separability, an all-level recurrence, a statement about arbitrary Keller
maps, or any general Jacobian-conjecture conclusion.

## 5. Exact companion

The companion

```text
04-computation/keller_G_renewal_faces_independent_probe_20260816.py
```

reconstructs `J` directly from the canonical global norm route rather than
importing THM-3506's face script.  It pins the full `J` ledger, checks the two
hybrid faces by direct scan and endpoint intersection, computes every
relevant invariant initial form, reduces the inverse leading coordinates
modulo both residual cubics, retains the Vieta leading coefficient, and pins
both rational face scalars.  Ordinary and optimized replay are required to
match the stored transcript exactly.

The LF SHA-256 hashes of the script and output are respectively
`f9e82f502026dfe499ebba9290295f98056d1b7dba7c893184d9871a032be01f` and
`becaa80c075bd46e4193b216406c2152f3d5d8565f6116ba6db9b712409badaa`;
the semantic ledger hash is
`2887be137141414185ab305e7e0416754f73a009625924a5b6c6c1a268101dbd`.
An independent parent replay reproduced that semantic ledger exactly in
about four minutes ten seconds.  This is an exact replay audit, not an
all-level or fifth-image audit.
