---
id: THM-3645
title: "Berggren two-cube local prime support and mod-24 cones"
status: >
  PROVED + LOCAL-ALGEBRA + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Every mod-prime obstruction to the slope conic is supported on nVT.
  A screen through n forces exactly two mod-24 candidate cones.  The first
  slope escaping the fixed p<=997 screen is (512,1019), not the previously
  reported (1012,1039).  No p-adic or Pell-orbit sufficiency is claimed.
source: kps-s189 / THM-3640 local-screen structure, 2026-08-21
audit: >
  PASS after two wording repairs -- an independent reconstruction proved the
  full support-prime iff theorem, exceptional primes, mod-24 cones, and
  screen-through-n completeness; replayed all finite counts and digests;
  checked all 1375 candidate-prime pairs through (n,p)=(31,97) and all 65769
  nonzero residue pairs with zero disagreements; and verified the first
  fixed-screen escape and infinite CRT family.  Normal and optimized streams
  match the stored transcript after LF normalization.
depends_on:
  - THM-3640-berggren-positive-cube-slope-atlas-through-401
related:
  - THM-3375-berggren-positive-two-cube-pell-ray
  - THM-3594-berggren-positive-cube-slope-atlas-through-201
script: 04-computation/berggren_two_cube_local_prime_support_mod24_thm3645.py
output: 05-knowledge/results/berggren_two_cube_local_prime_support_mod24_thm3645.out
script_sha256: 0193e4dfe3a6124985ac8108b9daf2473ace85b45f569cf6f3e486f243c4fc68
output_sha256: ec67a7c70dbbbe7fb6ebc43acccde79f386da115d82e01c445b88ae0f83a0b7a
semantic_sha256: 3844db20b97ae505d83bd2131daabb9eec87bbfcdb725258f188299624e0b4b4
hash_basis: raw LF bytes for files; canonical JSON for semantic ledger
---

# THM-3645 -- two-cube local prime support and mod-`24` cones

**PROVED + LOCAL-ALGEBRA + FINITE-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  Consider a primitive parity-correct slope

```text
n odd,  m even,  n/2<m<n,  gcd(m,n)=1.                 (1)
```

Put

```text
U=n-m,       V=2m-n,       T=2m+n.                     (2)
```

Then `U,V,T` are positive odd integers,

```text
gcd(U,V)=1,       4m^2-n^2=VT,                          (3)
```

and `n,V,T` are pairwise coprime.  The mod-prime screen used in THM-3640 is
the solubility of

```text
3W^2 = n^2 V T z^2 + B        (mod p),
B=4(2m^2+2mn-n^2).                                     (4)
```

## 1. Complete support theorem

For every odd prime `p!=3`, `(4)` is soluble unless `p|nVT`.  At a support
prime its exact criterion is

```text
p|n or p|V:       (6/p)=1,
p|T:              (-2/p)=1.                            (5)
```

At the two exceptional primes,

```text
p=2:  always soluble,
p=3:  soluble iff m=n!=0 (mod 3).                       (6)
```

This is an if-and-only-if classification of **mod-prime** solubility, not a
claim about solubility over every `Q_p`.

To prove it, first note from `(1)--(3)` that an odd prime divides at most one
of `n,V,T`.  If `p` divides none of them, the left side after moving the
`z^2` term is a nondegenerate binary quadratic form over `F_p`.  Every such
two-dimensional form represents every scalar: in the anisotropic case this
is surjectivity of the norm `F_{p^2}^* -> F_p^*`, and in the isotropic case it
is immediate after splitting the form.  Thus `(4)` is soluble.

If `p|n` or `p|V`, substitution in `(4)` gives a square equation equivalent
to `(6/p)=1`.  If `p|T`, it gives a square equation equivalent to
`(-2/p)=1`.  Direct reduction of the eight primitive residue pairs modulo `3`
gives `(6)`; modulo `2`, `(W,z)=(0,0)` works.

## 2. Why a screen through `n` is complete

Suppose `(4)` survives every prime `p<=P` and `P>=n`.  The `p=3` condition
forces

```text
3|U,       3 does not divide nV,       3|T.              (7)
```

Also `0<V<n`, while

```text
T/3=(2m+n)/3<n.                                         (8)
```

Consequently every prime divisor of `n`, `V`, or `T` other than the already
handled prime `3` is at most `n`.  By Section 1 there can be no unseen larger
mod-prime obstruction.  Thus the finite screen is complete whenever `P>=n`.

## 3. The two residue cones

The first condition in `(5)` puts every prime factor of `nV` in

```text
p=1,5,19,23 (mod 24),                                   (9)
```

so `n,V` lie in the same four-element subgroup modulo `24`.  The second puts
every prime factor of `T` in `1` or `3 mod 8`, hence

```text
T=1 or 3 (mod 8).                                       (10)
```

Equation `(7)` says `n=V mod 3`, and `m=(n+V)/2` even says that `n,V` have
opposite residues modulo `4`.  The four remaining pairs are

```text
(n,V)=(1,19),(19,1),(5,23),(23,5) (mod 24).             (11)
```

The first two make `T=2n+V` equal to `5` or `7 mod 8`, contradicting `(10)`.
Therefore every completely screened slope lies in exactly one of

```text
(n,V)=(5,23) or (23,5) (mod 24).                        (12)
```

Equivalently, the primitive candidate universe is reduced to the two affine
cones

```text
m=12a+24b+14,  n=24a+24b+23,
m=12a+24b+26,  n=24a+24b+29,                            (13)
```

subject to positivity and `gcd(U,V)=1`.  These are necessary cones only.

## 4. Exact censuses and the fixed-screen boundary

The structural screen independently reproduces the complete THM-3640 local
ledger:

```text
n<=401: 8195 candidates, 104 survivors,
survivor digest 2b1e24ba...448ff8e.                     (14)
```

Extending the complete screen through `n<=997` gives

```text
50507 candidates, 481 survivors,
V=5 mod 24: 238,       V=23 mod 24: 243.                (15)
```

The first new surviving denominator after `401` is `431`, with the eight
local survivors

```text
(218,431),(230,431),(266,431),(290,431),
(314,431),(374,431),(398,431),(410,431).                (16)
```

None of `(14)--(16)` includes the generalized-Pell or compiler-orbit gate.

The condition `P>=n` is sharp.  For the fixed screen `p<=997`, the first
escape from `(12)` is

```text
(m,n)=(512,1019),       (U,V,T)=(507,5,2043).           (17)
```

Direct exhaustive residue evaluation shows that `(17)` passes all `168`
primes through `997`.  But `1019|n` and

```text
(6/1019)=-1,                                             (18)
```

so it fails at the first omitted denominator prime.  This corrects the
earlier exploratory report that called `(1012,1039)` minimal.

There are infinitely many escapes from the **fixed** screen, which should not
be confused with infinitely many admissible Pell slopes.  Let

```text
M=product of all primes p<=997,
U_s=507+2Ms,  V_s=5,
m_s=512+2Ms,  n_s=1019+4Ms       (s>=0).                (19)
```

Every member is primitive and parity-correct, has the same residues as `(17)`
at each screened prime, and still has `(n_s,V_s)=(11,5) mod 24`.  Hence all
members pass the fixed screen while violating `(12)`.

## 5. Reproduction and boundary

Run

```text
python3 04-computation/berggren_two_cube_local_prime_support_mod24_thm3645.py
python3 -O 04-computation/berggren_two_cube_local_prime_support_mod24_thm3645.py
```

For every candidate through `n=31`, the companion compares the least
structural obstruction through `97` with the least brute-force obstruction
through `97`; it also directly brute-checks `(17)` at all `168` screened
primes and at `1019`.  The hostile audit strengthens the small-universe check
to every candidate-prime and nonzero residue-prime pair.  Normal and optimized
streams agree byte for byte and reproduce the stored transcript after LF
normalization.

The theorem classifies the mod-prime gate only.  It proves neither `Q_p`
solubility nor generalized-Pell class existence, compiler-orbit admissibility,
an infinite family of admissible slopes, or a density statement.
