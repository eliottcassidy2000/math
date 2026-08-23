---
id: THM-3706
title: "LRC Dedekind checksum and unbounded Bockstein depth"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  THM-3679's scalar-address seam
  has no row-independent finite 13-adic death depth, even among positive
  primitive nine-coordinate rows with the fixed strict blocker valuation
  profile (1,2,3).  Separately, the genuine Cover14 deep-well ray
  {1,...,12,182m}, already proved lonely by canonical THM-732, has unbounded
  total-speed 13-adic valuation.  Along the same ray both the fixed observer
  s(14,182m+1) and the explicit lonely-phase observer
  s(14m,182m+1) have exactly one 13-adic digit for every 13-unit m.  This
  one-digit Dedekind versus arbitrarily deep checksum split extends
  arithmetically to n=p+1 for every prime p>=5.  These are no-go results for
  fixed-depth or these Dedekind-observer completion strategies, not an LRC
  counterexample or an LRC(14) proof.
source: codex-lrc14-20260822 / THM-3679 checksum crossed with the 182 far ray
audit: >
  PASS -- an independent reconstruction checked both typed families, the
  unique-lift count, Cover14 typing, the collision-safe THM-732 citation,
  both Dedekind reciprocity formulas, their rational valuations, the general
  prime analogue, raw hashes, and normal/optimized transcript identity.
depends_on:
  - THM-3679-lrc-p-adic-scalar-seam-and-total-speed-checksum
  - THM-732-disc-v-bernoulli-edge-pair-dedekind-form-exact-certificates-far-element-tail
  - THM-2057-scaled-zeta-core-one-tail-closure
related:
  - THM-566-lrc14-covering-sets-have-no-uniform-bounded-denominator-witness
script: 04-computation/lrc_dedekind_checksum_unbounded_depth_thm3706.py
output: 05-knowledge/results/lrc_dedekind_checksum_unbounded_depth_thm3706.out
script_sha256: 396d1df0195ac5067ae899314c6c376656f229e1abc89043e5319f2ccaa88f12
output_sha256: 28b417f88a774485307855adc4c32a543b2006deb437dde35853e5be41bf2d4b
hash_basis: raw LF bytes
---

# THM-3706 -- the Dedekind observer sees one digit of an unbounded seam

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**
THM-3679 kills its last scalar
address hostile at the row-dependent depth

```text
n_*=nu_13(W)+1,                 W=sum_i w_i.         (1)
```

This theorem proves that “row-dependent” is essential.  It then separates
three objects that must not be conflated: a typed nine-coordinate relation
row, a genuine thirteen-speed Cover14 row, and a classical Dedekind sum.

## 1. Unbounded depth in the exact THM-3679 row type

For `m>=1`, put

```text
w(m)=(1,2,3,4,5,11;13,26,182m).                   (2)
```

The first six coordinates are units modulo thirteen, the last three are
divisible by thirteen, the entries are positive and pairwise distinct, and
the row is primitive.  Its checksum is

```text
W(m)=13(5+14m).                                    (3)
```

For every `h>=1`, there is a positive `m` for which

```text
nu_13(W(m))=h.                                     (4)
```

Indeed, for `h=1` take `m=1`.  For `h>1`, the coefficient `14` is a unit
modulo `13^(h-1)`, so there is one residue class solving
`5+14m=0 mod 13^(h-1)`.  Of its thirteen lifts modulo `13^h`, exactly one
also solves the next congruence; any of the other twelve has exact valuation
`h-1`.  A positive representative proves (4).  THM-3679 then says the scalar
all-unit seam survives through depth `h` and dies exactly at `h+1`.

The conclusion persists after fixing a strict blocker-depth profile.  Put

```text
w_strict(m)
 =(1,2,3,4,5,11;
   11*13,12*13^2,13^2*182m),                       (5)
```

with `m` a 13-unit.  The three blocker valuations are exactly `(1,2,3)`,
and this is again a positive, pairwise-distinct, primitive row of exactly the
six-unit/three-blocker type.  Moreover,

```text
sum_i w_strict(m)=13^3(1+14m).                     (6)
```

For every `r>=0`, the same lifting argument gives a 13-unit `m` with
`nu_13(1+14m)=r` exactly.  Hence (6) has every valuation `3+r`, despite the
fixed blocker profile.  Small exact controls are

```text
ordinary family (nu_13(W),m):
(1,1),(2,8),(3,60),(4,1412),(5,10200),(6,238688),

strict (1,2,3) family:
(3,1),(4,25),(5,12),(6,4237),(7,2040),(8,716065).  (7)
```

Therefore no fixed finite bank of Bockstein layers can kill THM-3679's
scalar seam uniformly, even on this narrow positive primitive type.  This
does not weaken THM-3679: its finite death theorem is explicitly per row.

## 2. Cover14 itself does not bound the checksum depth

The preceding rows are finite-algebraic relation rows; no scalar-cover
semantics has been assigned to them.  A separate genuine Cover14 family
shows that the small-modulus covering predicate alone cannot repair the
uniform-depth proposal.  Let

```text
V_m={1,2,...,12,182m}.                              (8)
```

Every `V_m` is a primitive thirteen-speed Cover14 row: moduli `2,...,12`
have their own speed, and `182m` is divisible by both `13` and `14`.  Its
total speed is

```text
sum V_m=78+182m=13(6+14m).                         (9)
```

The lifting proof above again makes `nu_13(sum V_m)` arbitrarily large; the
constructed choices can be taken to be 13-units.  The first examples are

```text
(nu_13(sum V_m),m)
 =(1,1),(2,7),(3,72),(4,1255),(5,12240),(6,212167). (10)
```

Canonical THM-732 proves the entire deep-well ray (8) lonely by its exact
Bernoulli edge-pair certificate and far-element tail.  Thus these are not
counterexamples to LRC(14).  They are counterexamples to the claim that
Cover14, or even Cover14 plus already-known loneliness, forces a uniform
bound on the total-speed valuation.

The nine-coordinate families in Section 1 and the thirteen-speed family in
this section are separate examples.  Equation (9) does not transport a
THM-3679 current or identify the scalar address seam with THM-3701's scalar
mass line.

## 3. Both natural Dedekind observers see exactly one digit

Let `s(a,b)` denote the classical Dedekind sum with the sawtooth convention,
and put

```text
D_m=182m+1=14(13m)+1.                              (11)
```

The standard reciprocity law and `D_m=1 mod 14` give, directly,

```text
s(14,D_m)
 =91m(13m-14)/(6D_m).                              (12)
```

For completeness, reciprocity says

```text
s(14,D_m)+s(D_m,14)
 =-1/4+(1/12)(14/D_m+D_m/14+1/(14D_m)),
```

while `s(D_m,14)=s(1,14)=13/14`; simplification is (12).  If `m` is a
13-unit, then `D_m=1 mod 13` and `13m-14=-1 mod 13`.  Since `91=7*13`,

```text
nu_13(s(14,D_m))=1                                (13)
```

as a rational numerator-minus-denominator valuation.  This first sum is a
fixed numerator-14 observer, not the binding statistic at THM-2057's lonely
phase when `m>1`.

That distinction strengthens rather than weakens the no-go.  THM-2057's
explicit strict lonely phase on (8) is

```text
t_m=14m/D_m.
```

The corresponding phase-adapted Dedekind sum has the second reciprocity form

```text
s(14m,D_m)=91m(13-14m)/(6D_m).                    (13a)
```

For a 13-unit `m`, its displayed factor `13-14m` is also a 13-unit, so this
explicit lonely-phase observer again has valuation exactly one.  The lifted
choices realizing arbitrary depth can be taken to be 13-units, as are all
examples in (10).  Their total-speed checksums therefore have unbounded
valuation while both Dedekind valuations remain one.

At the first point on the ray,

```text
s(14,183)=-91/1098,
-12*183*s(14,183)=182.                             (14)
```

This rederives the older observed `182` Dedekind identity; the theorem does
not depend on a historical hypothesis or on an unproved modular-form
interpretation.

## 4. The split is uniform for n=p+1

The arithmetic mechanism is not special to thirteen.  Let `p>=5` be prime
and `m>=1`, and put

```text
n=p+1,
K=n(n-1)=np,
D_m=Km+1=n(pm)+1.                                  (15)
```

One reciprocity step gives, for positive `q`, the general identity

```text
s(n,nq+1)=nq(q-n)/(12(nq+1)),

s(n,D_m)=npm(pm-n)/(12D_m).                        (16)
```

For every positive p-unit `m`, all factors in the second expression except the
displayed `p` are p-units, so

```text
nu_p(s(n,D_m))=1.                                  (17)
```

The phase-adapted numerator `nm` satisfies `D_m=p(nm)+1`, and the same
identity with its two arguments exchanged in role gives

```text
s(nm,D_m)=pnm(p-nm)/(12D_m),                       (17a)
```

again of exact p-adic valuation one.

Meanwhile the natural Cover-`n` row with `n-1=p` speeds

```text
{1,2,...,p-1,Km}                                  (18)
```

has total

```text
p((p-1)/2+(p+1)m).                                 (19)
```

The row in (18) is positive, pairwise distinct, and primitive.  Each modulus
`2,...,p-1` has its own speed, while `K=np` covers the remaining moduli `p`
and `n`; hence it is genuinely Cover-`n`.

The coefficient `p+1` is a p-unit, and the root modulo `p` is nonzero.
The same unique-lift argument gives every positive p-adic valuation in (19)
while keeping `m` a p-unit.  Thus the one-digit Dedekind versus unbounded
checksum split is a prime-neighbour phenomenon, not a numerical accident at
`14`.

## 5. Consequence and boundary

This is the one-prime analogue of THM-566's bounded-denominator obstruction.
THM-566 loads every denominator in a proposed fixed rational-witness bank by
an LCM multiplier.  Here one lifted residue class loads arbitrarily many
digits of a proposed fixed p-adic bank.  Both mechanisms require the eventual
observer to adapt to the row or to use semantic structure stronger than the
bare Cover14 predicate.

The strongest survivors are:

```text
THM-3679: every fixed positive row dies at its own finite depth;
THM-3706: those death depths are unbounded, even at blocker profile (1,2,3);
THM-732: the Cover14 deep-well ray is nevertheless lonely;
Dedekind: both fixed and explicit lonely-phase statistics record only the first digit. (20)
```

Nothing here supplies the nonzero common current, cross-source owner
intertwiner, one-packet marginal, or cover-specific mass/endpoint inequality
still required by the live proof map.  In particular, it neither constructs
an LRC counterexample nor proves LRC(14).

The standard-library-only exact companion pins THM-3679 and canonical THM-732,
checks the displayed families, evaluates five Dedekind sums directly, tests
the lifting mechanism at sixty `(p,depth)` controls for twelve primes, and
matches under normal and optimized execution.

```bash
python -B 04-computation/lrc_dedekind_checksum_unbounded_depth_thm3706.py
python -B -O 04-computation/lrc_dedekind_checksum_unbounded_depth_thm3706.py
```

**QED.**
