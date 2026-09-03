# LRC(14) partial-attachment address nonfactorization

**Status: FINITE-EXACT REFLECTION / REPRODUCIBLE SIDECAR, NOT A
THEOREM.** **LRC(14) remains open.** The two rows below are previously known
safe controls. They audit one explicitly defined quotient; they do not exclude
a counterexample class, prove an arbitrary-height bound, or promote the collar
calculation to canon.

## Scope and inheritance

The exact semantics come from:

- THM-4330,
  `01-canon/theorems/THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve.md`,
  for the anchored `h=420` hostile controls;
- THM-4335,
  `01-canon/theorems/THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md`,
  for addressed open-tooth, farthest-reach renewal;
- THM-4345,
  `01-canon/theorems/THM-4345-lrc14-halfperiodic-anchor-strip-euclidean-remainder-and-current-envelope.md`,
  for the exact half-periodic current consumer; and
- THM-4348,
  `01-canon/theorems/THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow.md`,
  for the centered-residue successor convention and nesting warning.

The closest proved mechanism is addressed shortest-cover renewal. The
canonical hostile is the THM-4330 pair `P=1287,9009`. The corrected near miss
is that owner/status/count equality does not identify a partial physical
itinerary. The least-used sidecar is the located endpoint occurrence on a
missing component.

The result below is a **finite nonfactorization witness only for the declared
coarse quotient and the declared next-blocker/exit/current consumers**.

## Declared finite universe

Fix

```text
h=420,
D_i=11+1680i,                         0<=i<=6,
(C_0,C_1,C_2,C_3)=(525,945,1365,1575),
B={D_0,...,D_6,C_0,...,C_3},
W_P=B union {P},                       P in {1287,9009}.
```

The physical speed rows are `{840} union W_P`. These are the two THM-4330 /
THM-4335 safe controls in the live `420|h` anchored `2+12` setting. They
meet the inherited denominator filters and defeat the displayed
half-turn/unit/capacity gate package, but each already has a safe witness.

For each of the `2h=840` anchor-safe components `I_k`, let

```text
epsilon=floor(k/h),     j=k mod h.
```

The exact greedy trace uses open tail teeth and the inherited farthest-right
tie order.

## Source, quotient, and full attachment map

The source object is the family of exact component traces for a fixed row
`W_P`. A trace contains the selected tooth occurrences up to either a span,
renewal, or first uncovered exit.

The tested coarse quotient retains the multiset

```text
(epsilon, final status, completed role word)
```

and, separately, the final status at every physical component. On a missing
component its word is replaced by `()`: its component remains known through
the separate status map, but its partial word and endpoint itinerary are
erased.

Before quotienting, the row parameter `P` and its speed-to-role dictionary
are retained, and every selected outgoing endpoint has the record

```text
(h, epsilon, source owner role, side=R, occurrence ordinal,
 k, j, source tooth address n,
 next owner role, next tooth address m, final status).
```

Thus the unquotiented source preserves owner, side, occurrence, anchor clock,
phase sheet, physical component, and tooth address. Together with the row's
speed dictionary, it recovers the metric endpoint.

## Exact hostile pair

For `P=1287` and `P=9009`, the common coarse quotient is

```text
epsilon  status   completed role word       multiplicity
0        missing  ()                         363
0        span     (D0)                        55
0        renew    (C2,C0,D0)                   1
0        renew    (D0,C3)                      1
1        missing  ()                         363
1        span     (D0)                        55
1        renew    (C3,D0)                      1
1        renew    (D0,C0,C2)                   1
```

Its exact `repr` digest is

```text
33ac4293f77a205ca52441a1cd387d9c8877ba6d44fe9687e1c3d59319dc41d3
```

The componentwise status maps are identical:

```text
missing=726, span=110, renew=4.
```

All `110` completed spans and all four completed renewals agree, including
their addressed teeth. Nevertheless:

```text
changed component traces                    209
changed partial role words                  208
changed nonmissing component traces           0
attachment-record symmetric difference     2229
```

Every changed trace ultimately has status `missing`. This localizes the
quotient loss to erased *partial* attachment information. The larger raw
symmetric-difference count includes downstream occurrence-ordinal shifts; it
is not a count of 2,229 independent physical changes.

### Small exact consumer witness

For component `k=9`,

```text
I_9=[127/11760,139/11760],     epsilon=0, j=9.
```

Both rows reach the same source occurrence: the right side of
`C3=T(1575,17)`, occurrence ordinal zero, at

```text
x=239/22050.
```

Its selected successor changes:

```text
P=1287: C3(17) -> P(14) -> D5(92) -> C2(15) -> exit,
P=9009: C3(17) -> D4(73)                    -> exit.
```

Equivalently, the same addressed source

```text
(row speed map,420,epsilon=0,C3,R,occurrence=0,k=9,j=9,n=17)
```

maps to `(P,14)` in the first row and `(D4,73)` in the second. The exact
local exits are

```text
P=1287: 211/19110,  clearance 1/14, binding speed 1365;
P=9009: 1023/94234, clearance 1/14, binding speed 6731.
```

On the common wall-free intervals adjacent to `x`, with neighboring walls

```text
1021/94234 < x < 766/70637,
```

the lower depth, upper depth, and signed current are

```text
                 left of x       right of x
P=1287           (3,1,2)         (2,1,1)
P=9009           (2,1,1)         (1,1,0).
```

The global exact anchor-core energies also differ:

```text
E_1287 = 103565349065690759276041319 / 79794403580513459456309691,
E_9009 = 58621224727881646861397 / 45136410657303198493260,

E_1287-E_9009
 = -11655225554668397928097 / 13640068988121958881420460.
```

The exact `q`-cubic consumer differs by

```text
-394136575599552666161942 / 1196916053707701891844645365.
```

Therefore this coarse quotient does not factor the next-blocker, local-exit,
one-sided-current, global-current, or `q`-cubic consumers.

### Safe-control negative result

The globally first missing component in the ordered scan is unchanged:

```text
k=5, witness=939/141274, clearance=1/14, binding speed=10091.
```

All changed traces are missing, while every completed cover is unchanged.
Thus attachment changes do **not** automatically change every consumer. In
particular, this experiment does not justify storing every address
unconditionally, and it supplies no new safety claim: both inputs were safe
before the comparison.

## The `7u|h` collar calculation

**Status of this subsection: PROVED ELEMENTARY CALCULATION INSIDE THIS
REFLECTION, CHECKED EXACTLY ON A FINITE GRID AND THE LIVE INSTANCE, BUT NOT
PROMOTED TO A CANON THEOREM.** The general identities are proved by direct
substitution; the script's finite loop is a hostile check, not the proof of
their universal quantifier.

Let `h=7Lu`, with positive integers `L,u` and odd `u`. At a signed
`u`-wall

```text
t_0=(14n+sigma)/(14u),                 sigma in {-1,+1},
```

the anchor is zero. Put

```text
delta=1/(28h)=1/(196Lu),
t_in =t_0-sigma delta,
t_out=t_0+sigma delta.
```

Direct substitution yields

```text
||(2h)t_in||=||(2h)t_out||=1/14,
||u t_in||   =1/14-1/(196L),
||u t_out||  =1/14+1/(196L),
||13u t_in|| =1/14+13/(196L),
||13u t_out||=1/14-13/(196L).
```

Thus `u` kills the inward attachment and `13u` kills the outward
attachment. The anchor-safe component width is

```text
W_A=3/(49Lu).
```

The full-component containment slacks are

```text
u inward:   (28L-13)/(196Lu)>0                    for L>=1,
13u outward:(28L-169)/(2548Lu)>0                  iff L>=7.
```

Hence for `L>=7`, the two adjacent safe components are entirely covered,
inward by `u` and outward by `13u`. The outward slack is negative at
`L=6`, namely `-1/(15288u)`, so `L=7` is sharp for this particular
two-component containment calculation. This blocks only the naive
two-boundary desingularization; it does not prove the row unsafe elsewhere.

For the live values `h=420,u=3,L=20`,

```text
t_0=1/42,
t_in=93/3920,     u-distance=279/3920<1/14,
t_out=281/11760,  13u-distance=267/3920<1/14,
inward slack=547/11760,
outward slack=391/152880.
```

There are `4u=12` distinct adjacent collar components at this scale.

## Connection contract

```text
source:
  exact THM-4335/4348 greedy traces for the two h=420 safe controls

target:
  the declared completed-cover quotient and occurrence-labelled
  endpoint-to-component attachment map

map:
  retain phase/status/completed support and componentwise status;
  replace every missing partial word by (), erasing its itinerary

preserved by the quotient:
  abstract role support, h=420 clock, phase sheet, strict/open convention,
  componentwise final status, completed spans/renewals, multiplicities

destroyed by the quotient:
  missing-prefix owner sequence, endpoint occurrence and tooth address,
  successor, arrival itinerary, metric wall position, local exit,
  signed current and nonlinear q-cubic value

required sidecar for the demonstrated consumers:
  row/speed map plus
  (h,epsilon,role,side,occurrence,k,j,n,next_role,next_n,status)

hostile:
  P=1287 versus P=9009; equal quotient, unequal k=9 successor,
  local exit and current. The global first ordered exit remains equal.

cheapest decisive replay:
  evaluate the successor of b(1575,17) on component 9 and the
  one-sided current around 239/22050.
```

The full attachment map is not claimed sufficient for LRC(14). In fact,
identical completed renewal attachments coexist with unequal currents here,
so metric occupation remains a separate, consumer-dependent coordinate.

## Next decisive `828`-component test

Use the actual `h=420,u=3` row fragment with `13u=39`, and adjoin ten odd
tails satisfying the inherited denominator/capacity gates. Delete the `12`
collar components covered by the calculation above. On the remaining `828`
components, construct the exact bipartite multigraph

```text
endpoint occurrence (speed,n,side,epsilon)
       -> anchor component (k,j)
       -> selected successor / exit.
```

Record current and first exit before quotienting. The decisive outcomes are:

1. equal coarse states with unequal residual exits, proving another address
   coordinate necessary;
2. fibre-constant consumers after collar deletion, licensing a compression;
3. a uniformly missing residual component, giving a genuine safety
   certificate for that finite family.

The missing mathematical input is a counterexample-conditioned restriction,
or a finite realization, for the ten residual tails. The algorithmic
trace/current machinery already exists. Unconditional reuse bounds are not
available because the THM-4348 nesting families remain compatible with safe
rows.

## Reproduction

From the repository root:

```bash
python -B 04-computation/lrc14_partial_attachment_address_nonfactorization_probe_20260902.py
python -B -O 04-computation/lrc14_partial_attachment_address_nonfactorization_probe_20260902.py
python -B 04-computation/lrc14_partial_attachment_address_nonfactorization_independent_referee_20260902.py
python -B -O 04-computation/lrc14_partial_attachment_address_nonfactorization_independent_referee_20260902.py
```

A primary run with `PYTHONHASHSEED=913` also matches. The scripts explicitly
request LF stdout so normal, optimized, hash-seeded, and frozen output bytes
match on Windows as well as LF-native systems.

```text
04-computation/lrc14_partial_attachment_address_nonfactorization_probe_20260902.py
  ba8f5a8c88fd0b28bd7dcd223901e48b039851f8fef7660bce37de8b0faeaa77
04-computation/lrc14_partial_attachment_address_nonfactorization_independent_referee_20260902.py
  d4a12dcd58133e5234d5863975a07aa4f9701d6edb34c678a3e6e9c05ca6cc12
05-knowledge/results/lrc14_partial_attachment_address_nonfactorization_probe_20260902.out
  c63a7192ea37eb95c283d64d766daba6963e858fe065227f15611393c46024c4
05-knowledge/results/lrc14_partial_attachment_address_nonfactorization_independent_referee_20260902.out
  6077796f0ce2b6a613f295c6a7a7121bcddf7ce03ab98e78e274c64d21f419fb
```
