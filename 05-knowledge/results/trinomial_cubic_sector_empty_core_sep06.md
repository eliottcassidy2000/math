# Exact overlap with the real-rooted-core duplication sector

**Status: PROVED ANALYTICALLY + FINITE-EXACT controls.** This note determines
which first-cancellation roots of two fixed trinomial families satisfy the
real-rooted-core hypothesis of
[THM-4440, signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md).
It does not extend that theorem to a nonreal-rooted core.

The source is the real-rooted polynomial coefficient implication
`[s^k]H=0 => [s^(2k)]H²<0`, with `ord_0 H<k<deg H`.
The target is the actual anchored doubled sign for the
[quadratic family](trinomial_width15_empty_core_returns_sep06.md) and
[cubic family](trinomial_cubic_empty_core_sep06.md). The closest mechanism is
THM-4440; the canonical missed coordinate is the factor `tau^-1` between
canonical second-row normalization and the first anchor. The least-used
sidecar is the location of the first-row roots relative to the core's
real-rooted sector. The live concepts are carry gauge, core discriminant,
first-row root placement, weighted trace inertia, and translated positivity.

## 1. Exact gauge and sign transfer

For `a=15` or `21`, consider

```text
f(u)=alpha*u^-a+beta*u^(2g-a)+gamma*u^(3g-a),
alpha*beta*gamma!=0,       tau=alpha*gamma²/beta³.
```

Choose any nonzero `kappa` with `kappa^g=beta/gamma` and put
`mu=beta^-1*kappa^(a-2g)`. Then the exact scalar/variable gauge is

```text
f_tilde(u)=mu*f(kappa*u)=u^-a*(tau+u^(2g)+u^(3g)).
```

It preserves `tau`. With `R(s)=tau+s²+s³`, one has
`CT(f_tilde^g)=[s^a]R(s)^g` and
`CT(f_tilde^(2g))=[s^(2a)]R(s)^(2g)`.
The admissible families have `0<a<3g`, and `R(0)!=0`.
Consequently THM-4440 applies to `H=R^g`, `k=a`, whenever `R` is
real-rooted and the first moment vanishes.

Let `X` be the actual first channel's coefficient monomial. Since this
channel has total mass `g` and charge zero, its transform is
`X_tilde=mu^g X`. In the displayed gauge it equals a nonzero integer
power of the real scalar `tau`. Thus

```text
CT(f_tilde^(2g))/X_tilde² = CT(f^(2g))/X²,
X_tilde²>0  when tau<0.
```

The strict real coefficient sign from THM-4440 therefore transfers to
the actual first-anchor quotient. In both families this quotient is
`tau^-1 Q_g(tau)`. The carry factor is retained; the canonical `Q_g`
alone has the opposite sign at a negative root.

The cubic discriminant is

```text
disc_s(tau+s²+s³) = -tau*(4+27tau).
```

Hence, for real `tau<0`, the exact real-rooted sector is
`-4/27<=tau<0`. At its left endpoint
`R(s)=(s+2/3)²(s-1/3)`, which is allowed: THM-4440 does not require
simple roots. Reversing the core gives `1+s+tau*s³` and the same sector.
For a real core on exponents `{0,2,3}`, every scalar/variable gauge
normalizes by real operations to this same core, so such a gauge cannot
enlarge the sector. This is a statement about the compressed cubic
real-core consumer; a different theorem on a different object is not
excluded.

The map preserves the rootwise doubled sign and constant-term zero
questions. It does not preserve arbitrary raw coefficient signs. The
needed sidecar is the core discriminant, followed by root placement.
The cheapest decisive test is exact evaluation of the first polynomial
at `-4/27`, together with monotonicity throughout the sector.

## 2. The quadratic first-row family

For `g>=8`, `gcd(g,15)=1`, put `u=g-5>=3`. The normalized first row is

```text
P2(tau)=6tau²+20u*tau+u(u-1).
```

Its two roots are distinct and negative by the quadratic-family proof.
Its derivative is positive on `[-4/27,0]`, while `P2(0)>0`. At the left
endpoint,

```text
F2(u)=P2(-4/27)=u²-(107/27)u+32/243,
F2(3)=-670/243,       F2(4)=68/243.
```

The polynomial `F2` is strictly increasing for `u>=4`. Thus exactly one
first root lies in the sector for `g=8`; none does for any other
admissible integer `g` (the next is `11`). There is no boundary equality.
In particular, the all-height negative trace/positive norm certificate
adds both root signs beyond the real-core consumer for every admissible
`g>=11`.

## 3. The cubic first-row family

For `g>=11`, `gcd(g,21)=1`, put `u=g-7>=4`. The normalized first row is

```text
P3(tau)=72tau³+504u*tau²+84u(u-1)tau+u(u-1)(u-2).
```

The separate cubic proof establishes three distinct negative roots.
On `[-4/27,0]` one has

```text
P3''(tau)=432tau+1008u>0,
27P3'(-4/27)=2268(u-4)²+11844(u-4)+11216>0.
```

Thus `P3` is strictly increasing there. Since `P3(0)>0`, its sector-root
count is one or zero according as the following endpoint value is
negative or positive:

```text
F3_boundary(u)=u³-(139/9)u²+(2066/81)u-512/2187.
```

This symbol denotes the endpoint cubic, not the discriminant factor
called `F3` in the main proof. The two endpoint values are

```text
F3_boundary(4)=-177848/2187,
F3_boundary(13)=-178820/2187.
```

Its derivative is an upward quadratic, positive at `0`, negative at
`4`, and positive at `13`. It therefore has one zero in `(0,4)` and
one in `(4,13)`; only a local minimum occurs on `[4,13]`. The maximum
on that interval is at an endpoint, and both endpoint values are
negative. For `u=14+s`,

```text
F3_boundary(14+s)=s³+(239/9)s²+(14666/81)s+161272/2187>0
```

for every `s>=0`. For integer parameters there is no equality or
unclassified gap. Exactly one first root lies in the real-core sector
for

```text
g=11,13,16,17,19,20;
```

none does for any admissible `g>=22`. At least two cancellation roots
of every support, and all three for `g>=22`, are outside the real-core
hypothesis. Their negative doubled signs require the separate
three-minor certificate. This is an exact scope separation on actual
first-row roots, not a counterexample to THM-4440.

## 4. Exact sidecar and boundary of the conclusion

The [source](../../04-computation/trinomial_cubic_sector_empty_core_sep06.py)
and [saved output](trinomial_cubic_sector_empty_core_sep06.out) check all
endpoint identities, derivative signs, the positive shift, and named
exact Sturm counts. The finite list is `g=8,11,13` in the quadratic
family and `g=11,13,16,17,19,20,22,23` in the cubic family. The interval
and all-height statements are proved above; the finite controls do not
stand in for those proofs.

```sh
python3 -B 04-computation/trinomial_cubic_sector_empty_core_sep06.py
python3 -B -O 04-computation/trinomial_cubic_sector_empty_core_sep06.py
```

General higher-channel negativity outside the real-core sector remains
**OPEN**. The surviving connection is a complementary sector theorem
plus a complete quotient trace certificate. Neither real-rootedness of
the first row nor coprimality alone implies the missing rootwise sign.

The sidecar has 55 explicit gates; normal and optimized outputs are
byte-identical. Semantic digest
`f6d97101bd54257e2cee92743863aa409c08cc035082e719c3abd64004d68908`.
Raw SHA-256 source/output:

```text
0a4565ea2b852e8342502ac40b0aa75e81d4586bce2be14c0b73945ff973e177
45d4bad25282b5d8a6985a51e95a96f28ddf474940f24ebe3dd1798903774b82
```
