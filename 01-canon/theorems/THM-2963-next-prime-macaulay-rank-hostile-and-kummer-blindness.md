---
id: THM-2963
title: "Next-prime Macaulay rank hostile and Kummer blindness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The proposed
  all-width next-prime rank-at-least-35 criterion is
  false: at depth zero, support (0,5,6,10), and p=41, the full degree-seven
  Macaulay map has rank exactly 34, with two explicit projective common roots
  and a nonzero 34-minor. Two further rank-34 cells occur in the frozen
  8,879-case universe. All factorial arguments are below the selected prime,
  so every Kummer carry vanishes and the initial form is the whole map.
  The next prime restores rank 36 in all three cells, but this is finite-exact
  evidence only. The valid theorem is the one-way THM-2947 certificate: any
  denominator-safe prime with rank at least 35 proves characteristic-zero
  nonvanishing on the physical no-real-support locus. No uniform one-prime,
  two-prime, or rank-at-least-34 theorem is claimed.
source: codex-gmc-next-prime-rank-hostile-2026-07-29
audit: >
  Independent reconstruction, without importing the companion, reproduces
  the three rank-34 cells, the explicit p=41 roots and determinant-12
  minor, the 8,879-case census and histogram, and the repaired one-way
  THM-2947 implication. Normal, optimized, and stored transcripts, declared
  LF hashes, stderr, and documentation checks all pass.
depends_on:
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
related:
  - THM-2949-fixed-fifth-compound-cofactor-and-width-three-twelve-closure
  - THM-2955-width-twenty-fixed-cofactor-modular-depth-gate
  - THM-2957-first-gap-widths-fifteen-twenty-modular-depth-ladder
script: 04-computation/gmc_all_width_modular_rank_hostile_thm2963.py
output: 05-knowledge/results/gmc_all_width_modular_rank_hostile_thm2963.out
script_sha256: e31072e4c2f60e0f74abf7f2db9f1e00d541c9e771dde158d3df2afc82e1e74c
output_sha256: 81318604fafcf2482ae59ba475f3f3adbdeb5e0eb6e1e37447ecb8f0f4eaa5ac
hash_basis: working-tree bytes (LF)
---

# THM-2963 -- next-prime Macaulay rank hostile and Kummer blindness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This note audits the proposed all-width statement

> for support `(0,a,b,M)` at depth `n>=0`, if `p` is the first prime
> strictly larger than `4(n+M)`, then the full degree-seven Macaulay map
> has rank at least `35` over `F_p`.

That statement is false already at width ten.

## 1. Inheritance and the lost coordinate

The closest proved mechanism is
`THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate`.
Over the reals, absence of real support pairs all local factors by complex
conjugation, so the degree-seven rank lies in

```text
{30,32,34,36}.
```

Thus one real rank-`35` witness forces rank `36`. Reduction modulo a prime
preserves a nonzero minor, but it destroys this conjugation pairing: a
finite-field closed point may have degree one, so modular rank `35` is
possible, and two degree-one points or one degree-two point give rank `34`.

The least-used tempting sidecar was Kummer carry control. Here it points in
the opposite direction. The condition `p>4(n+M)` puts every factorial
argument strictly below `p`, so all relevant `p`-adic factorial valuations
are zero. There is no Kummer channel suppression: the prime initial form is
the whole reduced moment map.

## 2. Independent constructor

The companion

```text
04-computation/gmc_all_width_modular_rank_hostile_thm2963.py
```

does not import any fixed Macaulay chart or cofactor. For

```text
support=(n,n+a,n+b,n+M)
```

it:

1. constructs the four-variable factorial moment forms directly;
2. substitutes the fourth coefficient as `-(x+y+z)`;
3. constructs every `Q*S_5`, `C*S_4`, and `F*S_3` row of the full
   `46 x 36` degree-seven map; and
4. performs exact finite-field row reduction.

The monomial orders are generated recursively by ordinary weak
compositions, independently of the fixed-cofactor scripts.

## 3. Minimal exact hostile

Take

```text
n=0,       (a,b,M)=(5,6,10),       p=41.
```

This is exactly the prescribed prime because `41` is the first prime
larger than `4(n+M)=40`.

In recursive weak-composition monomial order, the direct reduced forms are

```text
Q:  6,37,37,35,0,9

C:  15,14,19,10,30,21,9,9,12,10

F:  3,34,7,6,28,0,35,34,3,35,20,21,30,19,3.
```

Their frozen payload digest is

```text
dc29c15792f64a9759c3a450626ab2a413ad197a4cb9217f1de09f13b811b586.
```

As a non-truth-bearing compatibility check, evaluating the independent
THM-2943 denominator-cleared constructor at the same cell gives the same
three forms up to the nonzero scalars

```text
Q:30,                  C:21,                  F:27 mod 41.
```

Its full matrix also has rank `34`. Thus the hostile is not caused by a
normalization or ordinary/divided-power convention mismatch.

There are two distinct common projective roots:

```text
P_1=(1:22:0),              P_2=(1:39:38)       in P^2(F_41).
```

Evaluation at `P_1` and `P_2` gives two independent annihilating covectors
on `S_7`. Hence the full Macaulay row rank is at most `34`.

The reverse inequality has a literal minor certificate. With global row
order

```text
Q rows 0..20, C rows 21..35, F rows 36..45,
```

the rows

```text
0,1,2,3,4,5,21,22,6,7,8,9,10,23,24,11,12,13,14,
25,27,15,16,17,26,29,18,28,19,36,35,37,38,39
```

and target columns `0..33` have determinant

```text
12 mod 41.
```

Therefore

```text
rank_F41(Phi_7)=34                                      (1)
```

exactly. The `Q,C` submap has rank `30`, so this is a genuine length-two
intersection with the quartic, not degeneration of the quadratic/cubic
complete intersection.

Every factorial argument used above is at most `40`; hence every
`41`-adic factorial valuation is zero. This explicitly refutes a
carry-based proof of modular corank at most one under the stated prime
choice.

## 4. Finite hostile census

The normal and optimized executions agree byte-for-byte on the frozen
universe

```text
all n=0 supports,                3<=M<=30;
all depths 1,M-1,M,2M,          3<=M<=15;
3000 seeded random cases,       3<=M<=120, 0<=n<=2000.
```

After duplicate removal this is `8,879` cases, with

```text
rank 34:       3
rank 35:      89
rank 36:   8,787
rank <=33:     0.
```

The three rank-`34` cases are exactly

```text
(n,a,b,M;p)=(0,5,6,10;41),
             (0,4,7,15;61),
             (0,8,16,19;79).                         (2)
```

The first and third consist of two rational common points:

```text
p=41: (1:22:0), (1:39:38);
p=79: (1:29:51), (1:49:20).
```

The `p=61` case has no rational common point. Since its `Q,C` submap has
rank `30` and the full corank is two, it is one degree-two closed point.

This supports, but does **not prove**, the weaker experimental bound

```text
rank_Fp(Phi_7)>=34.                                   (3)
```

No theorem is claimed from `(3)`.

An additional exhaustive scout at `n=0`, `31<=M<=50`, checked `15,540`
supports and found rank histogram

```text
rank 35: 82,             rank 36: 15,458,
```

with no further rank-`34` case. This extra range is scout evidence rather
than part of the frozen transcript.

## 5. Second-prime repair: exact finite evidence, not a theorem

No alternate chart can repair `(1)`: it is the rank of the **full** map,
so every `35`-minor vanishes modulo `41`.

The next prime does repair every case in `(2)`:

```text
p=41 ->43: rank 36,
p=61 ->67: rank 36,
p=79 ->83: rank 36.                                  (4)
```

Thus the two-successive-prime bank covers the frozen `8,879`-case
universe. There is currently no proof that it works for every width and
depth. Two different primes can both divide a nonzero integer resultant,
so `(4)` cannot be promoted by a coprimality slogan alone.

The theorem-grade repaired criterion is only:

> If **some** prime avoiding the factorial denominators gives
> `rank(Phi_7)>=35`, then some integer `35`-minor is nonzero. Over the
> physical real no-support locus, THM-2947 upgrades rational rank at least
> `35` to rank `36`, hence the genuine resultant is nonzero.

This is a sufficient certificate bank. It is not an equivalence and not a
uniform one-prime or two-prime theorem.

## 6. Strongest survivor and next decisive test

The failed connection is:

```text
source:     carry-free factorial reduction mod p
target:     degree-seven Macaulay initial form
preserved:  exact coefficients and any surviving minor
lost:       real conjugate-pair parity
failure:    two finite-field common points can survive simultaneously
sidecar:    a second prime, or a prime-independent resultant divisor bound
test:       search simultaneous rank<=34 at the first two primes.
```

The best structural next target is therefore not “Kummer forces corank one.”
It is one of:

1. prove a finite prime-bank separation theorem for the integer
   `35`-minor ideal/resultant;
2. bound the common prime divisors of two calibrated minors; or
3. find a simultaneous two-prime hostile.

The present exact hostile should be kept as the mandatory regression test
for any such theorem.

## 7. Reproduction

```text
python 04-computation/gmc_all_width_modular_rank_hostile_thm2963.py
python -O 04-computation/gmc_all_width_modular_rank_hostile_thm2963.py
```

Both reproduce
`05-knowledge/results/gmc_all_width_modular_rank_hostile_thm2963.out`.

LF-normalized SHA-256:

```text
script  e31072e4c2f60e0f74abf7f2db9f1e00d541c9e771dde158d3df2afc82e1e74c
output  81318604fafcf2482ae59ba475f3f3adbdeb5e0eb6e1e37447ecb8f0f4eaa5ac
```

The independent audit reconstructed the three hostile cells directly from
the factorial moment forms, including the degree-two closed point over
`F_61`, and independently recovered the complete rank histogram.  It also
checked that denominator clearing, rather than a characteristic-zero
conjugation argument modulo `p`, is exactly what makes the rank-`>=35`
criterion one-way.
