---
id: THM-4055
title: "Sixty-phase/dyadic response fibre law and fixed Rule 30 firewall"
status: >
  PROVED universal phase-fibre law + FINITE-EXACT fixed-bank application +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. Pairing a C_60 address
  with a least dyadic response phase C_(2^e) gives the fibre product over
  C_(gcd(60,2^e)); every fixed sixty-phase fibre contains exactly
  2^max(e-2,0) response-phase states. A sharp torsor sidecar reconstructs
  those states, while the response itself factors through C_60 iff e<=2.
  Applied to THM-4047's fixed Rule 30 bank, exactly 82 of 100001 certified
  column tails factor through C_60; the other 99919 require one, two, or
  three bits for lossless phase reconstruction. In particular ell_29 has
  ell_29(t+60)=ell_29(t)+1 in F_2 for every t>=90. This does not control the
  moving center, temporal balance, or any Rule 30 prize.
source: long-precise-frontiers / 2026-08-24
audit: >
  PASS. The primary path enumerates the complete combined clocks for all six
  dyadic periods in the certified bank, checks the fibre-product image,
  sidecar bijection, representative carry law, frozen THM-4047 hash and
  histogram, and reconstructs ell_29 from the closed packed front strip. An
  independent path enumerates compatible CRT pairs, inverts the sidecar
  modulo each fibre, reads THM-4047's independent output, and evolves ordinary
  spatial Rule 30 rows. Both paths certify a closed-strip repeat from time 90,
  and include a least-period-eight hostile that
  separates phase states from scalar response values. Normal and optimized
  streams byte-match the frozen outputs and the semantic digests agree.
depends_on:
  - THM-4047-rule30-left-front-affine-monodromy-clock
related:
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
  - THM-4035-sixty-clock-separation-and-finite-kakeya-spine
  - THM-4038-ap-deficit-holonomic-sixty-phase-law
  - THM-4042-prime-sector-ap-cover-exact-clock-and-holonomic-law
  - THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall
  - THM-4048-rule30-periodicity-balance-and-model-firewalls
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
script: 04-computation/sixty_dyadic_response_fibre_thm4055.py
output: 05-knowledge/results/sixty_dyadic_response_fibre_thm4055.out
script_sha256: c6de48e956c425b21d20f4dc364e10d41e3fd7d7401c486c0aee08e416c1f535
output_sha256: e6ed3dfa72f85e33d4b702f1cfcaa48965624f44e8a2f458601b5cd44d58728f
independent_audit_script: 04-computation/sixty_dyadic_response_fibre_thm4055_independent_audit.py
independent_audit_output: 05-knowledge/results/sixty_dyadic_response_fibre_thm4055_independent_audit.out
independent_audit_script_sha256: e7ed7c344aec61d98d6cf39f6102c5bd20f218fe9b765f6b5dfe9fc7157d5646
independent_audit_output_sha256: ee5f6b95286b47cde49308a5a128044dcaa0fe75ea2a18d4a23573a9d6bba0f5
semantic_sha256: 4de3162dfeaaa77c75256bdd86592d46110745992ff575f2c8c8953d378a55f0
hash_basis: raw LF bytes
---

# THM-4055 -- the exact sixty/dyadic fibre law

**PROVED universal law; FINITE-EXACT fixed-bank application; VERIFIED-EXACT
and INDEPENDENTLY AUDITED.** A sixty-phase address and a dyadic response do not
share a clock once the dyadic period exceeds four. The exact common object is
a fibre product, and its lost coordinate is a finite torsor whose size is
sharp for phase reconstruction.

## 1. Universal phase-fibre theorem

Let `w:Z -> A` be a response with least period

```text
p=2^e,       e>=0.
```

Put

```text
g=gcd(60,p)=2^min(e,2),
m=p/g=2^max(e-2,0),
L=lcm(60,p)=60m.                                      (1)
```

The simultaneous phase map is an isomorphism onto the fibre product

```text
C_L -> C_60 x_(C_g) C_p,
[t]_L |-> ([t]_60,[t]_p),                             (2)
```

where the right side consists exactly of pairs `(r,a)` satisfying
`r=a mod g`. Consequently every fixed `r in C_60` has exactly `m` compatible
`p`-phases.

Choose the canonical integer lift `r_0 in {0,...,59}`. The sharp lossless
sidecar on that fibre is

```text
q=(t-r_0)/60 mod m.                                   (3)
```

It reconstructs both combined and response phase:

```text
t=r_0+60q mod L,       a=r_0+60q mod p.               (4)
```

Conversely, for a compatible pair `(r,a)`, choose an integer lift `a_0` of
`a`. When `m>1`,

```text
q=((a_0-r_0)/g)*(60/g)^(-1) mod m.                    (5)
```

For `e>=3`, equation `(5)` is

```text
q=15^(-1)*(a_0-r_0)/4 mod 2^(e-2).                    (6)
```

For `e<=2`, `m=1` and the sidecar is simply `q=0`; no inverse modulo one is
invoked.

The origin in `(3)` depends on the lift: replacing `r_0` by `r_0+60h`
translates `q` by `-h`. Thus the fixed-phase fibre is canonically a
`C_m`-torsor, not a canonically pointed copy of `C_m`. In canonical
representatives its addition law contains the carry

```text
(r,q)+(s,u)
 = ((r+s) mod 60, q+u+floor((r+s)/60) mod m).         (7)
```

For `e>=3`, this is not a direct-product group splitting: `gcd(60,m)>1`, so
`C_60 x C_m` is noncyclic, whereas `C_L` is cyclic.

The response has a phase-only evaluator `w(t)=F(t mod 60)` if and only if

```text
p divides 60  <=>  e<=2  <=>  m=1.                   (8)
```

When `m>1`, at least one fixed-sixty-phase restriction has least sampled
period `m`. Therefore `(3)` is phase-minimal: any sidecar which reconstructs
the `p`-phase on every fixed `r` needs at least `m` states, and `q` attains
that bound.

### Response-value guardrail

The `m` states in the theorem are **phase states**, not necessarily `m`
distinct scalar values. A least-period-eight binary word

```text
10000000
```

is constant on the `r=1` fibre but nonconstant on the `r=0` fibre. Thus a
particular evaluator may compress some fibres. Least global period guarantees
only that at least one fibre retains the full sampled period. The phase
sidecar is universally lossless; it is not asserted to be scalar-minimal on
every fibre.

## 2. Exact fixed-column Rule 30 consequence

THM-4047 proves, for the single-seed fixed left-front columns
`ell_r(t)=a_t(-t+r)` with `0<=r<=100000`, the all-future least-period census

| least period `p` | 1 | 2 | 4 | 8 | 16 | 32 |
|---:|---:|---:|---:|---:|---:|---:|
| columns | 16 | 10 | 56 | 668 | 87118 | 12133 |
| phase states over fixed `C_60` phase | 1 | 1 | 1 | 2 | 4 | 8 |
| lossless phase-sidecar bits | 0 | 0 | 0 | 1 | 2 | 3 |

By `(8)`, exactly

```text
16+10+56=82                                             (9)
```

certified fixed-column tails factor through `C_60`. The remaining

```text
668+87118+12133=99919                                  (10)
```

do not. Their exact response phase is recovered by adjoining respectively
one, two, or three bits to the sixty-phase address. These bit counts are
minimal for phase reconstruction; binary response values can of course use
fewer states on individual fibres.

The first period-eight doubling column gives a fully response-visible
hostile. THM-4047 certifies `ell_29` from time `90`, and its least-period word
indexed by absolute residues `0,...,7` is

```text
(1,1,0,1,0,0,1,0).                                   (11)
```

Thus `ell_29(t)=W_(t mod 8)` for `t>=90`, where `W` is the displayed word.

The second half is the bitwise complement of the first. Since
`60=4 mod 8`,

```text
ell_29(t+60)=ell_29(t)+1       in F_2,  for every t>=90. (12)
```

Here both phase states over **every** fixed sixty-phase address give distinct
response values. Equation `(12)` is therefore a sharp all-future witness that
the sixty address loses one dynamically relevant bit.

## 3. Proof

The generalized Chinese remainder theorem identifies the image of
`C_L -> C_60 x C_p` with the pairs agreeing modulo `g`; both sides have
`L=60p/g` elements, proving `(2)`. Fixing `r` leaves `p/g=m` choices. If
`t=r_0+60k`, then two integers `k,k'` give the same `p`-phase exactly when

```text
p divides 60(k-k')  <=>  m divides k-k'.             (13)
```

This proves `(3)--(5)` and phase minimality. Equation `(7)` follows by
reducing `r+s` to its canonical representative.

If `p` divides `60`, the sixty phase determines the `p`-phase, so `w`
factors. Conversely, if `w` factors through `C_60`, then `60` is a period of
`w`. Every period is divisible by its least period `p`, proving `(8)`.

For completeness, suppose `m>1` and every fixed-phase restriction had period
properly dividing `m`. Since `m` is a power of two, every such period divides
`m/2`; hence `60m/2` would be a period of `w`. For `e>=3`,

```text
60m/2 = 15p/2 = p/2 mod p,                            (14)
```

contradicting the least period `p`. At least one restriction therefore has
least sampled period `m`.

Equations `(9)--(12)` now follow from THM-4047's independently audited exact
period census and its all-future certificate for column `29`. No finite tail
fit is used.

## 4. Validity boundary and cross-clock firewall

This theorem compares address maps only. It does not identify their causes or
their consumers.

- THM-4038's consecutive seven-sector AP law has a `C_60` **template** clock,
  but its evaluated deficit retains the unbounded height `n`. Its persistent
  owners are rational mechanical/Christoffel words; “Sturmian 60-periodic
  sequence” is not the precise object. THM-4042 also blocks extrapolating the
  denominator-lcm explanation to arbitrary prime-sector laws.
- THM-4035's pointed Fibonacci and triangular states can relabel `C_60`, but
  their scalar shadows lose phase and neither supplies Rule 30's dyadic
  sidecar.
- THM-4035's broad twisted-cubic spine and planar hostile share the same
  phase. Their different four-fold determinants come from representation
  weights, not from the clock.
- Under THM-4044's field hypotheses, its polynomial observer has kernel
  `((P^60-1)^k)` and loses a boundary Hasse jet. That ideal-theoretic loss is
  not the torsor loss in `(3)`.
- THM-4027 proves that THM-4026's Sun counterexample target is represented
  modulo every modulus. A sixty-clock residue atlas cannot restore its missing
  integer lift height.

Finally, THM-4047's center is the moving diagonal

```text
c_t=ell_t(t).                                         (15)
```

The fixed-column certificate gives no uniform onset bound reaching that
diagonal. THM-4050's marked stopping cylinder and THM-4048's observer/carry
identities likewise do not turn `(12)` into temporal periodicity, balance,
bounded radius, or a query lower bound. All three Rule 30 prizes remain
**OPEN**.

## 5. Replay

From the repository root:

```text
python3 -B 04-computation/sixty_dyadic_response_fibre_thm4055.py
python3 -B -O 04-computation/sixty_dyadic_response_fibre_thm4055.py
python3 -B 04-computation/sixty_dyadic_response_fibre_thm4055_independent_audit.py
python3 -B -O 04-computation/sixty_dyadic_response_fibre_thm4055_independent_audit.py
```

Each normal/optimized pair is byte-identical, and the independent semantic
digest agrees with the primary digest. **QED.**
