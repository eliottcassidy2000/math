---
id: HYP-2175
status: SUPPORTED by S613 quotient audit; theorem decomposition open
source: user-2026-06-03; codex-2026-06-03-S613
related:
  - HYP-2174
  - HYP-2172
  - HYP-2169
  - HYP-2168
  - HYP-2167
  - HYP-2166
  - HYP-2164
  - HYP-2156
  - HYP-2153
  - HYP-2154
  - HYP-2135
  - THM-401
  - THM-407
---

# HYP-2175: LRC dimension descent is residue-depth descent, not runner-count descent

## Claim

The useful way to reduce the n=14 LRC proof toward n=13, n=12, n=11, and
smaller parameters is not to imagine a literal projection of runner sets.  The
modulus changes from

```text
C = 2n - 1,
```

so the salient proof size is the arithmetic of `C`: unit-shell count, nonunit
shell strata, prime-power depth, and the number of D/N/pinch survivors after
the THM-401/THM-407 quotient pipeline.

The descent sequence near n=14 is therefore:

```text
n=14: C=27=3^3  -> depth-3 nonunit tower, gcd strata 1/3/9
n=13: C=25=5^2  -> depth-2 nonunit tower, gcd strata 1/5
n=12: C=23 prime -> clean unit-shell reset
n=11: C=21=3*7  -> squarefree composite re-entry, gcd strata 1/3/7
```

Thus n=12 is a better proof model than n=11 despite being larger.  Conversely,
n=11 is smaller but less clean because nonunit obligations reappear.

## S613 Evidence

`04-computation/lrc_dimension_descent_salience_s613.py` applies the same
unit-shell quotient, D/N gates, and pair-sum pinch gate for every `3 <= n <= 14`.
It recovers the known n=14 S608/S610 counts exactly:

```text
n=14, C=27: unit_rows=340928, D/N survivors=27733,
             strict=27730, floor=3, below=0
```

The adjacent descent ledger is:

```text
14->13: C 27=3^3 -> 25=5^2, survivors 27733 -> 4883, floor 3 -> 2
13->12: C 25=5^2 -> 23 prime, survivors 4883 -> 420
12->11: C 23 prime -> 21=3*7, survivors 420 -> 3415
11->10: C 21=3*7 -> 19 prime, survivors 3415 -> 146
 9-> 8: C 17 prime -> 15=3*5, survivors 80 -> 220, floor 2 -> 4
```

Every audited D/N survivor has a strict or floor pair-sum pinch; there are zero
below-floor quotient rows for `3 <= n <= 14`.

The floor atoms also line up with the collapse-family story:

```text
n=14: AP, V*=AP[12->24], 2*AP
n=13: AP, 2*AP
n=12: AP, 2*AP
n=11: AP, 2*AP
n= 8: AP, AP[6->12], (1,4,5,6,7,11,13), 2*AP
n= 6: AP, (1,3,4,5,9), 2*AP
n= 5: AP, (1,3,4,7), 2*AP
```

So the sporadic `p_0=0` additive chains are not noise.  They are the same
residue-depth phenomenon seen at small `C`, where nonunit shell structure and
all-order cancellation align just enough to leave extra floor atoms.

The incoming HYP-2169 doubling-sporadic law sharpens one slice of this:
`AP[(n-2)->2(n-2)]` is tight exactly when `3|(2n-1)`.  S613 is compatible with
that law: the doubled floor atoms at `n=8` and `n=14` occur precisely at
`C=15,27`, while the audit keeps the other sporadics visible as off-law
collapse atoms.

## Structural Reading

The n=14 quotient tower of HYP-2166 is:

```text
integer rows
  -> Res_27 shadow plus carry fiber
  -> unit/nonunit shell quotient under <2,-1>
  -> D/N gates
  -> pair-sum pinch quotient
  -> owner/certificate reattachment.
```

S613 says the reusable part of this tower is the obligation ledger, not the
literal residue set.  There is no natural proof-preserving projection

```text
Res_27 -> Res_25 -> Res_23 -> Res_21,
```

because the coimage object changes with `C`.  What can descend is the theorem
shape:

```text
prime C          : pure unit-shell theorem;
prime-power C    : one prime channel with depth e;
squarefree C     : several independent nonunit channels;
lift over C      : carry-fiber conservativity, as in HYP-2167.
```

In the coimage/Yoneda language, the residue quotient is the coimage and the
D/N/pinch/owner ledgers are probes.  Smaller n is useful only when the probes
factor through a simpler coimage.

## Proof Leverage

The next proof decomposition should be by `C=2n-1` type.

1. Prime reset lemma.  For prime `C`, all shells are unit shells.  The quotient
   is exactly one choice from each antipodal shell when the row size is `n-1`.
   In the audit, prime `C=23,19,13,11,7,5` leaves only AP and `2*AP` floor rows;
   `C=17` has two unit orbits but the same floor outcome.

2. One-prime tower lemma.  For `C=p^e`, the proof burden is the depth of the
   nonunit chain.  `C=25` has only AP/`2*AP`; `C=27` adds the V* floor atom and
   the carry-fiber issue.  HYP-2169 isolates the mod-3 doubled-sporadic slice
   of this tower.  The theorem to prove is not "n=14 is large" but "depth-3
   one-prime nonunit towers are conservative after owner reattachment."

3. Squarefree composite lemma.  `C=21` and `C=15` show nonunit obligations can
   reappear when n decreases.  They should be treated as independent prime
   channels rather than as a nested carry tower.

4. Collapse-family bridge.  The small sporadics `(1,3,4,7)`,
   `(1,3,4,5,9)`, `AP[6->12]`, and `(1,4,5,6,7,11,13)` are floor atoms of the
   same quotient audit.  Their relation-lattice all-order cancellation should
   be studied as the low-dimensional model of V*.

## Tournament Analysis

S613 uses n-level quotient summaries as tournament vertices.

```text
vertices: n=3..14
observable: (below, survivors, unit_rows, primitive_floor, floor,
             nonunit_shells, shell_orbits, n)
switch: lower proof burden first, with n as tie Hamiltonian path
score histogram: degrees 0..11 once each
SCCs: singleton/transitive
directed 3-cycles: 0
edge flips versus naive dimension descent: 3/66
Hamiltonian path easiest->hardest:
  n=3 -> n=4 -> n=5 -> n=6 -> n=7 -> n=9 -> n=10 ->
  n=8 -> n=12 -> n=11 -> n=13 -> n=14
```

The flips are exactly the salient warning: prime/composite changes can outrank
runner count.

## Honest Status

This is not a proof of LRC n=14.  It is a proof-routing theorem shape supported
by exact finite quotient audits through n=14.

What is proved by computation:

```text
1. the n=14 quotient counts match S608/S610 exactly;
2. all audited D/N survivors for 3<=n<=14 are strict or floor by pair-sum pinch;
3. the main nonmonotone events are precisely composite/prime changes in C=2n-1;
4. the small p_0 collapse sporadics occur as quotient floor atoms.
```

What remains open:

```text
turn the C-type decomposition into symbolic lemmas, especially the
prime-reset lemma, the one-prime-power tower lemma, and the n=14 carry-fiber
conservativity theorem of HYP-2167.
```

## See

`04-computation/lrc_dimension_descent_salience_s613.py`,
`05-knowledge/results/lrc_dimension_descent_salience_s613.out`,
`07-reflections/lrc-dimension-descent-salience-s613.md`,
`05-knowledge/hypotheses/INDEX.md` (HYP-2169),
`05-knowledge/hypotheses/INDEX.md` (HYP-2172),
`05-knowledge/hypotheses/INDEX.md` (HYP-2174),
`05-knowledge/hypotheses/HYP-2166-lrc-n14-res27-quotient-tower-conservativity.md`,
`05-knowledge/hypotheses/HYP-2167-lrc-n14-carry-fiber-conservativity.md`,
`01-canon/theorems/THM-401-pair-sum-sieve-modulus-is-2n-minus-1.md`,
`01-canon/theorems/THM-407-twisted-involution-shell-reduction-of-the-LRC-additive-residual.md`.
