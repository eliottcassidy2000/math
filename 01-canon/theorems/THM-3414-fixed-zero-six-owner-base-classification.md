---
id: THM-3414
title: "Fixed-zero six-owner base classification"
status: >
  PROVED structural theorem + COMPUTER-ASSISTED PROVED global finite-profile
  atlas + INDEPENDENTLY AUDITED.  For every q>=2, a fixed-source-centre
  zero cover of all cyclic sheets by at most six transverse owners exists if
  and only if q is divisible by one of 15,16,18,20,24.  The reverse direction
  is all-q: it uses THM-3408's exhaustive quotient-order cutoff, a globally
  quantified arbitrary-companion lemma, and exact rational certificates over
  finitely many exceptional order profiles, not a bounded modulus scan.  No
  arbitrary common-centre, mobile-centre, physical-time, or LRC(14) conclusion
  is claimed.
source: root-2608-crouzeix-puzzle-2026-08-15
audit: exact standard-library Fraction simplex with direct post-verification; 368 admissible exceptional pairs, 64 pair-pruning lcm groups, 172 anchor lcm groups, 9 exceptional-family lcm groups, 14940 residual six-exceptional-owner multisets, positive dilation controls, and q21/q22/q102 hostiles; independent all-m/stratum/multiset/padding/dilation/scope audit and normal/-O replay clean
depends_on:
  - THM-3408-fixed-zero-additive-order-duality-and-six-core-corridor
  - THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight
related:
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3402-atomized-sheet-covers-and-constructive-cochain-locus
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
script: 04-computation/lrc_fixed_zero_six_owner_base_classification_thm3414.py
output: 05-knowledge/results/lrc_fixed_zero_six_owner_base_classification_thm3414.out
script_sha256: 57147225890f0d99fd9c0f13d70992dca2ef6113e34554bdff0b6df68f9e63e4
output_sha256: 48d2c34af30912f0d09cf34aad0a2e1e9d67fa6f781dc4383f06b2cab8369d91
semantic_sha256: 681e1b1f8066b54970ed791f14867df708e9808c8a66dce034002f0404e66b3f
certificate_sha256: bc080056351c09c2a68eb27ac7574db3d273fe8de85e7509c1bcc3497e34fe25
hash_basis: LF-normalized bytes
---

# THM-3414 -- fixed-zero six-owner base classification

**PROVED structural theorem + COMPUTER-ASSISTED PROVED global finite-profile
atlas + INDEPENDENTLY AUDITED.**

## 1. Statement

For `q>=2` and a transverse speed residue `u mod q`, put

```text
D_(q,u)(0)={ell in Z/qZ: ||u ell/q||<1/14}.                         (1)
```

Then

```text
there exist k<=6 transverse owners u_1,...,u_k with
union_i D_(q,u_i)(0)=Z/qZ

iff

15|q or 16|q or 18|q or 20|q or 24|q.                              (2)
```

Equivalently, the fixed-zero transverse sheet-cover rank satisfies

```text
r_0(q)<=6 iff q is divisible by a member of A={15,16,18,20,24}.     (3)
```

This is a theorem at the literal source centre `0`.  A common nonzero centre
is a different gauge problem by THM-3405 and MISTAKE-384.  Statement `(3)`
does not classify the full zero-cochain locus, mobile centres, physical-time
owners, or the LRC(14) decrement.

## 2. Inheritance and connection contract

THM-3401 gives covers of ranks `6,5,5,6,6` at the five bases in `A`.
THM-3408 supplies the exact additive-order density

```text
alpha(e),
e_q(m,n)=m/gcd(m,q/n),                                              (4)
```

for an owner of quotient order `m` on the sheet stratum of additive order
`n`, together with lcm descent and the exhaustive all-integer cutoff

```text
B={21,22,26,28,33,35,39,42,44,46,50,52,56,57,63,65,70,
   74,76,77,78,84,91,99,102,104,114,117,130,143,156,186,190},       (5)
H={22,33,44,46,50,102},                                             (6)
```

such that, for every base-free order `m>=2`,

```text
m notin B implies alpha(m)<=7/44,
max_(m in B) alpha(m)=1/5,                                         (7)
```

and every base-free cover by at most six owners contains an `H` owner.

| field | exact connection |
|---|---|
| source | a putative fixed-zero cover with quotient orders `m_i|q` |
| target | a finite family of rational profiles on strata `S_(q/d)(q)` |
| map | choose `L` from exceptional owners and send `m` to `m/gcd(m,d)` for proper `d|L` |
| preserved | exact owner mass on every selected stratum and the cover obstruction |
| destroyed | mode orientation and pointwise alignment inside a stratum |
| required sidecar | THM-3408's global outside-`B` bound and literal blocks for positive constructions |
| cheapest controls | `q=21` shared-prime loss, `q=22` corridor-without-cover, `q=102` sharp fractional equality |

The carrier is a rectangular owner-order/stratum-order payoff matrix.  There
is no intrinsic binary orientation, so tournament language would add no
content.

## 3. The globally quantified finite-companion lemma

Let `L` be base-free and let

```text
Delta(L)={d:d|L, d<L}.                                              (8)
```

For rational weights `lambda_d>=0` summing to one, define

```text
beta_(L,lambda)(m)
 =sum_(d in Delta(L)) lambda_d alpha(m/gcd(m,d)).                    (9)
```

If `L|q`, the term indexed by `d` in `(9)` is exactly the THM-3408 mass of
the owner on `S_(q/d)(q)`, because `q/(q/d)=d`.  The divisor is proper, so
`q/d>1`; every selected stratum is legitimate.

For any set `Delta subseteq Delta(L)`, put

```text
C(L,Delta)={m>=2:
  lcm(L,m) is base-free and, for some d in Delta and e in {1} union B,
  m|de and m/gcd(m,d)=e}.                                          (10)
```

This is finite.  More importantly, it controls **every** possible companion
order, with no bound on `q` or `m`.

Indeed, set `g=gcd(m,d)` and `e=m/g`.  If `e in {1} union B`, then

```text
m=eg,                 g|d,                  hence m|ed.             (11)
```

Thus an actual base-free companion with exceptional projected order belongs
to `(10)`.  If it does not belong to `(10)`, then for every support divisor
`d` its projected order is base-free and lies outside `{1} union B`; `(7)`
gives

```text
alpha(m/gcd(m,d))<=7/44.                                           (12)
```

Consequently a certificate need check only the finite rows `(10)`, plus the
single generic row `7/44`.  This is the load-bearing reason the computation
below proves an all-`q` statement rather than another finite census.

## 4. Unit forcing and pair pruning

Suppose `q` is base-free and a cover has at most six owners.  Call an owner
**exceptional** when its quotient order lies in `B`, and let `b` be the number
of exceptional owners.  On the unit sheet stratum, a family with `b<=1` has
total mass at most

```text
1/5+5(7/44)=219/220<1.                                              (13)
```

Hence `b>=2`.  THM-3408 also forces at least one of those owners to lie in
`H`.

There are exactly `368` unordered base-free-lcm pairs from `B`, with
repetition.  For a fixed pair `F=(x,y)`, let `L=lcm(x,y)`.  The exact atlas
groups `155` pairs into `64` common-`L` classes and supplies weights in `(9)`
and a rational `t` satisfying

```text
beta(m)<=t for every m in C(L,Delta(L)),
7/44<=t,
beta(x)+beta(y)+4t<1.                                               (14)
```

By Section 3, `(14)` bounds four arbitrary further owners, including orders
not in `B` and orders much larger than `L`.  Thus no cover can contain any of
those `155` pairs.  All inequalities are verified exactly; the largest left
side in this stage is

```text
7903/7920<1.                                                        (15)
```

The remaining `213` pairs form a graph with the following nontrivial maximal
cliques:

```text
(21,26,39,42,46,57,74,78,91,102,114,186)
(21,33,39,57,63,77,91,99,117,143)
(26,35,46,50,65,70,74,91,130,190)
(26,28,46,52,56,74,76,91,104)
(26,46,57,74,76,91,102,114,186)
(21,33,39,46,57,91,102)
(22,46,102)
(44,46,102).                                                       (16)
```

The orders `84` and `156` are isolated.  Repetition is allowed exactly at
the `29` looped vertices listed by the companion; in particular `22` and
`44` cannot repeat.  Every exceptional-order multiset in a putative cover
must be a multiset clique of this residual graph.

## 5. The anchor atlas

Fix the exact exceptional-owner count `b`, where `2<=b<=5`, and choose a
residual pair `F=(x,y)` containing an `H` owner.  Every other exceptional
order in the family is a common neighbour of `x,y`.  All remaining owners
are ordinary, meaning outside `B`.

For each common-`L` group of anchor pairs, the atlas supplies one probability
vector `lambda`, an ordinary bound `t_0`, and a neighbour bound `t_F` for each
pair, such that

```text
beta(m)<=t_0 for every finite exceptional projection of an ordinary m,
7/44<=t_0,
beta(z)<=t_F for every compatible residual common neighbour z,

beta(x)+beta(y)+(b-2)t_F+(6-b)t_0<1.                               (17)
```

Using `6-b` ordinary slots is safe even when the original family has fewer
than six owners, since every bound is nonnegative.  The exact inventory is

| `b` | certified anchor pairs | common-`L` groups | worst verified left side | uncertified anchor pairs |
|---:|---:|---:|---:|---|
| 2 | 60 | 45 | `1179/1232` | none |
| 3 | 60 | 45 | `6147/6160` | none |
| 4 | 57 | 42 | `121741/122760` | `(46,46),(46,102),(102,102)` |
| 5 | 55 | 40 | `173249/174240` | `(46,46),(46,74),(46,91),(46,102),(102,102)` |

The smallest atlas margin is therefore

```text
1-6147/6160=13/6160.                                               (18)
```

Exact residual-graph enumeration leaves only five `b=4` multisets and twenty
`b=5` multisets with no certified anchor pair.  The five are all size-four
multisets on `{46,102}`.  The twenty are the union of

```text
all size-five multisets on {46,102},
all size-five multisets on {46,74,91} containing at least one 46.   (19)
```

Their overlap is the all-`46` multiset.  Nine further common-`L` certificates
exclude `(19)`: three groups for `b=4`, with worst value `591/616`, and six
groups for `b=5`, with worst value `47/57`.

Finally, if `b=6`, there are no arbitrary companions.  The exact multiset
enumeration has `14,940` residual graph cliques containing `H`.  For each
family the companion finds a single proper divisor `d|L` such that

```text
sum_(m in family) alpha(m/gcd(m,d))<1.                              (20)
```

The largest minimum in `(20)` is only `1/2`, attained by

```text
family=(21,33,39,77,91,143),       L=3003,       d=3.               (21)
```

Sections 4--5 exhaust `b=2,...,6`, contradicting the fractional-cover
necessary inequality in every base-free case.  This proves the reverse
implication of `(2)`.

## 6. Positive direction

THM-3401 gives the following fixed-zero covers:

```text
15: (1,2,3,4,5,7)
16: (1,3,5,7,8)
18: (1,5,6,7,9)
20: (1,3,4,7,9,10)
24: (1,5,7,8,11,12).                                               (22)
```

If `a in A` divides `q` and `v` is a speed in the `a`-cover, lift it to

```text
u=(q/a)v.                                                          (23)
```

Then `u` remains transverse and

```text
||u ell/q||=||v ell/a||.                                           (24)
```

Thus every base block pulls back under `Z/qZ -> Z/aZ`; the lifted family
covers all of `Z/qZ` with at most six owners.  This proves the forward
implication and completes `(2)`.

## 7. Equality and failure boundaries

- `q=15` is the first positive rank-six edge and supplies the canonical
  trimode witness in `(22)`.
- `q=21` has exact fixed-zero rank eight.  Its unit owners lose the shared
  prime-order sheet strata; high unit density alone is not a cover.
- `q=22` lies in the forced `H` corridor but has exact rank seven.  Therefore
  membership in `H` is necessary in the base-free case, not sufficient.
- At `q=102`, THM-3408's exact fractional game has value `4/25`, with tight
  owner orders `(2,3,17,102)` and tight sheet orders `(6,34,51,102)`.  This is
  the sharp known fractional hostile and remains strictly below `1/6`.

The atlas proves only the rank-at-most-six threshold.  It does not assert the
conjectural all-`q` formula for the fractional value `tau(q)`, nor does it
reconstruct literal covers from density profiles.

## 8. Exact computation contract

The companion is standard-library only.  It:

1. rebuilds the `368` admissible pairs and verifies the residual graph and
   maximal cliques;
2. solves `64+172+9` rational profile programs with a deterministic exact
   two-phase simplex;
3. independently substitutes every returned `Fraction` vector into every
   finite candidate, generic, neighbour, and fixed-family inequality;
4. exhausts all `14,940` residual six-owner multisets and their proper
   divisor witnesses;
5. checks fifty positive dilation lifts, exact `q=21,22` literal ranks, and
   both sides of the `q=102` rational saddle certificate; and
6. freezes every rational support and bound in `certificate_sha256` and the
   full control surface in `semantic_sha256`.

Reproduce with

```bash
python 04-computation/lrc_fixed_zero_six_owner_base_classification_thm3414.py
python -O 04-computation/lrc_fixed_zero_six_owner_base_classification_thm3414.py
```

The two outputs must be byte-identical.  The computation has no modulus
cutoff: its only finite universes are the theorem-proved exceptional set `B`,
the exact finite-companion sets `(10)`, and their residual multisets.
