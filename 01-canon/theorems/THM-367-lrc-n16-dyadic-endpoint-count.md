---
id: THM-367
name: lrc-n16-dyadic-endpoint-count
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S391
depends_on:
  - THM-357
  - THM-360
  - THM-366
results:
  - 05-knowledge/results/lrc_n16_dyadic_endpoint_formula_s391.out
---

# THM-367: n=16 Dyadic Endpoint-Protection Count

## Statement

Work at Lonely Runner denominator `n=16`.  Let the owner speed be a pure
dyadic speed

```text
u = 2^k.
```

Its strict forbidden endpoints have labels

```text
e(m,eps) = (16m + eps)/(16u),
0 <= m < u, eps in {-1,+1}.
```

A protector speed `p` protects `e(m,eps)` exactly when

```text
dist_0(p*(16m + eps) mod 16u) < u,
```

where `dist_0` is circular residue distance to `0`.

Let `r = p mod 16u`.

1. If `r=0`, then `p` protects all `2u` endpoints.
2. If `r != 0`, write `r=2^j q` with `q` odd.  If `j >= k`, then `p`
   protects no endpoint.
3. If `j < k`, put `L=k-j`.  The number of protected endpoints is:

```text
L >= 3:                                      2^(k-2)
L = 2 and q mod 16 in {1,3,13,15}:          2^(k-1)
L = 1 and q mod 16 in {1,15}:               2^k
all other L=1,2 odd classes:                0.
```

Consequently, for maximal-owner lower protection `1 <= p < u`:

```text
u = 2,4,8  -> lower protectors do not cover all endpoints;
u = 16     -> lower protectors cover all endpoints, but need exactly 9.
```

For `u=16`, the unique necessary lower residue set is

```text
{1,3,5,7,8,9,11,13,15}.
```

Each of these nine residues has a private endpoint relative to the other
lower residues, and together they cover all `32` endpoints.

## Proof

The endpoint-protection criterion is THM-357 specialized to `n=16` and owner
`u`: a point `e(m,eps)` is strictly protected by `p` exactly when

```text
|p*(16m+eps) - a*16u| < u
```

for some integer `a`, which is the displayed circular residue condition.

If `r=0 mod 16u`, every residue is `0`, so all endpoints are protected.

Assume now `r != 0`, and write `r=2^j q` with `q` odd.  If `j >= k`, then
`r=2^k b` for some integer `b` not divisible by `16`, since `r` is not
`0 mod 16u`.  Modulo `16u=2^(k+4)`, the term `2^k b*16m` vanishes and the
endpoint residue is

```text
2^k b eps.
```

Its circular distance to `0` is at least `2^k=u`, because `b` is nonzero
modulo `16`.  The strict inequality `<u` is impossible.

It remains to count the case `j<k`.  Divide the strict inequality and modulus
by `2^j`.  With `L=k-j`, the reduced problem is modulo `2^(L+4)` with strict
radius `2^L`:

```text
dist_0(q*(16m+eps) mod 2^(L+4)) < 2^L.
```

Only `m mod 2^L` matters, and each reduced value of `m` has `2^j` lifts to the
original range `0 <= m < 2^k`.  For a fixed sign `eps`, as `m mod 2^L` varies,
the residues `q*(16m+eps)` are precisely the `2^L` residues congruent to
`q eps mod 16`.

For `L>=3`, in each odd congruence class modulo `16`, the strict radius
window around `0` contains `2^(L-3)` such residues.  The two signs give the
two classes `q` and `-q`, hence `2^(L-2)` reduced endpoints.  Lifting by
`2^j` gives

```text
2^j * 2^(L-2) = 2^(k-2).
```

For `L=2`, the reduced modulus is `64` and the strict radius window is

```text
{0, +/-1, +/-2, +/-3}.
```

Only the odd classes `+/-1` and `+/-3` occur.  Thus the two signs contribute
two reduced endpoints exactly when

```text
q mod 16 in {1,3,13,15},
```

and none otherwise.  Lifting gives `2^j * 2 = 2^(k-1)`.

For `L=1`, the reduced modulus is `32` and the strict radius window is

```text
{0, +/-1}.
```

Only the odd classes `+/-1` occur.  Thus the two signs contribute two reduced
endpoints exactly when

```text
q mod 16 in {1,15},
```

and none otherwise.  Lifting gives `2^j * 2 = 2^k`.

This proves the count formula.

For the lower-cover consequences, direct substitution into the formula gives
that owners `2,4,8` leave endpoints uncovered by all `p<u`; explicitly, the
endpoints with odd `m` remain uncovered.

For `u=16`, the lower residues `p=1,...,15` have the following private
endpoints relative to all other lower residues:

```text
p=1:   (1,-1), (15,+1)
p=3:   (5,+1), (11,-1)
p=5:   (3,+1), (13,-1)
p=7:   (7,-1), (9,+1)
p=8:   (2,+/-1), (6,+/-1), (10,+/-1), (14,+/-1)
p=9:   (7,+1), (9,-1)
p=11:  (3,-1), (13,+1)
p=13:  (5,-1), (11,+1)
p=15:  (1,+1), (15,-1).
```

Therefore any lower cover must contain all nine residues

```text
{1,3,5,7,8,9,11,13,15}.
```

A direct check from the same residue inequality shows those nine residues
cover all `32` endpoints.  Hence the exact minimum lower-cover size for
`u=16` is `9`.

## Interpretation

At denominator `16`, the no-gate branch is already closed by THM-366: without
a speed divisible by `16`, the odd unit points `a/16` survive as boundary
witnesses.  THM-367 describes the local debt created by entering the gated
branch.

The first dyadic owner whose lower protectors can close its endpoint layer is
`u=16`, and the closure costs nine lower residues.  This is the local form of
the S389/S390 Cayley-Dickson analogy:

```text
the 16-gate kills the old unit skeleton, but exports structured dyadic debt.
```

The computation also shows a warning for the full proof: higher pure dyadic
owners have self-similar nine-covers, so a proof cannot rely only on private
endpoints at every dyadic height.  It needs a global debt-flow invariant using
maximality, speed budget, and endpoint-layer divergence.

## Verification Record

`04-computation/lrc_n16_dyadic_endpoint_formula_s391.py` verifies the formula
by direct endpoint enumeration for `k=0,...,8`, checking all `8176` residues
`p mod 16u`.

The same script records:

```text
u=2,4,8: lower protectors do not cover all endpoints;
u=16:    exact lower-cover number 9;
u=32:    exact lower-cover number 9 in the bounded solver;
u=64:    exact lower-cover number 9 in the bounded solver.
```

It also prints a constructive nine-cover family for pure dyadic owners
`u=16,32,64,128,256,512`, and scans odd payloads `u=16w` for
`w in {1,3,5,7,9,15}`.

Stored output:

```text
05-knowledge/results/lrc_n16_dyadic_endpoint_formula_s391.out
```

## Related

- THM-357: Lonely Runner endpoint-protection trichotomy.
- THM-360: unit endpoint divisibility filter.
- THM-366: small-denominator divisibility sieve.
- HYP-1857: n=16 Cayley-Dickson dyadic debt law.
- HYP-1858: n=16 dyadic endpoint-debt proof route.
- `04-computation/lrc_n16_dyadic_endpoint_formula_s391.py`.
