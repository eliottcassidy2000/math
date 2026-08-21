---
id: THM-3663
title: "LRC minimal carry-closed exceptional observer bank"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  Three consecutive low-digit columns of translates of THM-3660's signed
  exceptional detector form a 39-coordinate injective observer for all 169
  addresses under both split and cyclic subtraction.  Vertical carry acts by
  a coordinate permutation and reversal by a signed coordinate permutation.
  No one- or two-column carry-closed bank is injective; exactly the thirteen
  consecutive three-column banks work.  The 11- and 14-edge boundary banks
  fail even to distinguish (0,0) from (0,1).  This is an address-level
  operation congruence, not a physical-current observer or LRC(14).
source: kps-s192 / THM-3661 adaptive-observer continuation, 2026-08-21
depends_on:
  - THM-3660-lrc-exceptional-leakage-functional-and-fourteen-edge-boundary
  - THM-3661-lrc-exceptional-detector-simple-spectrum-convolution-rigidity
related:
  - THM-3662-lrc-eleven-cell-exceptional-flux-and-high-digit-variation-gate
script: 04-computation/lrc_minimal_carry_closed_observer_thm3663.py
output: 05-knowledge/results/lrc_minimal_carry_closed_observer_thm3663.out
script_sha256: 3cbcad785da8a0231a9482c8d263bdc2fcde7dc153c0dbc3ed4766a6ef8b5966
output_sha256: 96b3346662ddb574ce7e20e2f88a364557c973cfc590b08a955eaed244999e91
semantic_sha256: 953994bb4440423e1569e8c35165a2a60c8d7c3fa438b92508203ab64f3085a9
hash_basis: raw LF bytes
---

# THM-3663 -- a minimal carry-closed observer has three columns

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  The exceptional detector can be turned into a complete finite
observer while making the native vertical carry and reversal operations
visible on observer coordinates.

## 1. Translate observers for both address laws

Let `g=1_(X_+)-1_(X_-)` be THM-3660's detector on the common 169-element
address set.  Equip that set with either

```text
G_sp=F13^2
```

or the assembled cyclic law

```text
G_cy=C169,       a=r0+13r1.                           (1)
```

For an offset bank `A` define, using the chosen group subtraction,

```text
O_A(x)=(g(x-a))_(a in A).                             (2)
```

Take

```text
A={-1,0,1} x F13,                 |A|=39.             (3)
```

Here `-1` denotes low digit 12; the high digit ranges over all 13 values.
Exact enumeration gives

```text
|{O_A(x):x in G_sp}|=169,
|{O_A(x):x in G_cy}|=169.                            (4)
```

Thus (2) loses no address information under either law.  The two codebooks
are not equal—split and cyclic subtraction remain distinct—but each is
injective.

## 2. Carry and reversal are coordinate actions

Let `v=(0,1)`.  The bank (3) is closed under subtracting `v`.  Hence

```text
O_A(x+v)_a=O_A(x)_(a-v),                              (5)
```

for both group laws.  Vertical carry is therefore one fixed permutation of
the 39 coordinates.

Branch reversal is the same address involution in both presentations:

```text
J(r0,r1)=(12-r0,12-r1).                               (6)
```

The bank (3) is also closed under the group inverse `a->-a` for either law.
Since `g(Jx)=-g(x)`,

```text
O_A(Jx)_a=-O_A(x)_(-a).                               (7)
```

Thus reversal is a signed coordinate permutation.  Equations (4)--(7) make
the observer an exact congruence for these two address-level operations,
unlike a finite separator whose successor requires discarded state.

## 3. Minimality among carry-closed translate banks

Any offset bank closed under `a->a+v` is a union of complete vertical
13-cycles, hence has the form

```text
A_C=C x F13                                  (8)
```

for a set `C` of low-digit columns.  The detector occupies only the five
columns

```text
S={0,3,6,9,12}.                                      (9)
```

If `|C|<=2`, the translated supports cover at most ten low-digit columns.
At least one entire vertical fibre of 13 addresses therefore has the
all-zero observer code, so injectivity is impossible.

For `|C|=3` there are only

```text
binom(13,3)=286                                      (10)
```

banks.  Exact exhaustion under each group law gives precisely thirteen
successful banks:

```text
C={c-1,c,c+1},          c in F13.                    (11)
```

No other three-column bank is injective.  Hence 39 coordinates are necessary
and sufficient among all carry-closed translate banks.

## 4. Hostile boundary banks and the unrestricted interval

It is tempting to use only the 11 carry-visible edges of THM-3662 or all 14
edges of THM-3660 as translation offsets.  Both fail under both laws:

```text
O(0,0)=O(0,1)=the all-zero code.                      (12)
```

So the sparse boundary is an excellent linear flux detector but not a
complete nonlinear address observer.

For an arbitrary bank of `m` translates, the union of their detector supports
has size at most `8m`.  Injectivity permits at most one all-zero address, so

```text
8m>=168,       hence m>=21.                           (13)
```

Together with (3), the unrestricted translate minimum currently lies in

```text
21<=m_min<=39.                                       (14)
```

No claim of unrestricted optimality is made.

## 5. A permutation-rigidity corollary

THM-3661 gives a useful operation-compatibility boundary.  For either group
law, let `C_g` be convolution by `g`.  If a permutation operator `U` on the
169 addresses commutes with `C_g`, simple spectrum gives

```text
U=P(C_g).                                             (15)
```

Every polynomial in a convolution operator is itself convolution.  A
convolution matrix that is a permutation matrix has a delta-function kernel.
Therefore

```text
U is a group translation.                             (16)
```

So compatibility with this one sparse observable leaves no exotic
permutation dynamics.  This is a finite group-algebra rigidity statement;
it does not identify any physical LRC evolution with such a permutation.

## 6. Exact verification and scope

Reproduce with

```bash
python3 -B 04-computation/lrc_minimal_carry_closed_observer_thm3663.py
python3 -B -O 04-computation/lrc_minimal_carry_closed_observer_thm3663.py
```

The assertion-free companion source-pins THM-3661, constructs both
codebooks, verifies every coordinate action in (5)--(7), exhausts all
carry-closed banks of one through three columns, and checks both hostile
boundary banks.  Normal and optimized streams are byte-identical.  The
complete semantic digest is

```text
953994bb4440423e1569e8c35165a2a60c8d7c3fa438b92508203ab64f3085a9.
```

This theorem lives entirely on the static 169-address set.  It does not
construct a THM-2334 current-to-address map, retain source character, owner,
phase, chronology, or gauge, prove first-return compatibility, force high-
digit variation, normalize a translated detector to fixed `X`, prove
exceptional entry, or imply LRC(14).  **QED.**
