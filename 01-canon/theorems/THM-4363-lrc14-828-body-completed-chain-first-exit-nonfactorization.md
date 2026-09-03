---
id: THM-4363
title: "LRC14 828-body completed-chain and first-exit nonfactorization"
status: >
  PROVED FINITE-EXACT FOUR-ROW NONFACTORIZATION + INDEPENDENTLY AUDITED;
  LRC(14) OPEN. For four explicit primitive h=420 minority-anchor rows,
  deleting twelve strictly covered collar components leaves 828 labelled
  components with identical status maps and identical physical chains on all
  282 completed components. Nevertheless their first missing-component exit
  records are pairwise distinct. The sharper P=761/1015 pair retains the same
  local role, owner, frontier, tie, and continuation packet but has physical
  teeth (761,21)/(1015,28) and different metric exits. Thus neither declared
  quotient determines its named consumer. All four rows have many safe
  translated half-turn clocks; no counterexample-family exclusion, ledger
  decrement, or general quotient theorem is asserted.
source: root + 828-body consumer scout + clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
  - THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow
related:
  - THM-4345-lrc14-halfperiodic-anchor-strip-euclidean-remainder-and-current-envelope
  - THM-4349-lrc14-adaptive-owner-clock-completeness-and-twelve-height-tail-cutoff
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
  - THM-4359-source-normal-row-eight-constructible-response-affine-modification
primary_script: 04-computation/lrc14_828_completed_chain_first_exit_nonfactorization_thm4363.py
primary_output: 05-knowledge/results/lrc14_828_completed_chain_first_exit_nonfactorization_thm4363.out
primary_script_sha256: d7497ce77698db8e99578a0bf184dcb15b27619dba11e9083ee86d966c3e77d8
primary_output_sha256: 28a8697c4f10fe6c832c13af43f1336134dfa42171a758754233fa2634f6b189
independent_referee_script: 04-computation/lrc14_828_completed_chain_first_exit_nonfactorization_independent_referee_thm4363.py
independent_referee_output: 05-knowledge/results/lrc14_828_completed_chain_first_exit_nonfactorization_independent_referee_thm4363.out
independent_referee_script_sha256: febd3555882384b8e4d32591755f396113896ab483eab0b6da3dd1a04a25629f
independent_referee_output_sha256: a93a554177502ae7cae6a41f61330db6646c02f9f8fda745c7b2c910448c5de2
hash_basis: raw LF bytes
audit: >
  PASS WITH SCOPE REPAIR. The primary compares the canonical renewal
  enumerator with a fresh exact interval implementation on all 4*828 retained
  instances and separately audits the twelve collars. An import-free referee
  rebuilds all teeth and cell-audits every wall interval on those 3,312 row-
  component instances. Normal/-O/hash-seeded/frozen LF streams match.
---

# THM-4363 -- LRC14 828-body completed-chain and first-exit nonfactorization

**PROVED FINITE-EXACT FOUR-ROW NONFACTORIZATION + INDEPENDENTLY AUDITED.
LRC(14) REMAINS OPEN. THESE ROWS ARE SAFE CONTROLS, NOT COUNTEREXAMPLES, AND
THE THEOREM LOWERS NO LRC LEDGER.**

## 1. Declared rows and endpoint convention

Use THM-4348's physical anchor components and open danger teeth:

```text
I_k=[(14k+1)/(28h),(14k+13)/(28h)],             0<=k<2h,
T(w,n)=((14n-1)/(14w),(14n+1)/(14w)).                     (1)
```

At a frontier, select the active tooth with farthest right endpoint; at an
equal right endpoint select the wider tooth, equivalently the smaller speed.
Open tooth endpoints are inactive, hence equality at distance `1/14` is safe.

Fix

```text
h=420,             anchor=2h=840,             u=3, 13u=39,
F=(11,1691,3371,5051,6731,8411,10091,525,945),
W_P=(3,39)+F+(P),
P in P0={241,255,761,1015}.                              (2)
```

Each `W_P` has twelve distinct positive odd speeds, and each thirteen-speed
row `{840} union W_P` is primitive.

The inheritance pass was:

- closest proved mechanism: THM-4335/4348's physical tooth address and exact
  farthest-reach renewal;
- canonical hostile: THM-2050's warning that complete local coordinates can
  miss the later global consumer;
- corrected near miss: failure of two fixed clocks and named finite gates
  does not imply failure of the complete translated clock grid;
- least-used sidecar: the selected prefix on a component whose terminal
  status is `missing`.

The live board was

```text
physical component | collar | status | completed chain | missing prefix
first exit | tooth address | Euclidean remainder | endpoint convention.    (3)
```

## 2. Named finite gates and the complete-grid firewall

All four rows represent every denominator `2,...,14`. At the two specific
half-turn clocks

```text
1/2 +/- 1/11760                                             (4)
```

their minimum distance is

```text
829/11760 < 1/14,                                          (5)
```

with unique blocker `5051`. The named doubled-denominator unit-bank
certificates for every `p=8,...,14` do not fire. The exact adaptive MC7
capacity audit over every represented divisor gives

```text
P       represented divisors   minimum capacity   minimizers
241             58                  10/7           7,21,35,105
255             61                  10/7           7,21,35,105
761             58                  10/7           7,21,35,105
1015            61                   9/7           7,35.             (6)
```

These are failures of sufficient certificates, not danger. Indeed, the
complete translated grid `1/2+j/11760`, `0<=j<11760`, contains respectively

```text
214, 222, 218, 170                                         (7)
```

safe points, with first displayed safe points
`121/240,121/240,281/560,281/560`. Thus the examples are explicitly safe;
no complete-grid failure is claimed.

## 3. Twelve strict collars and 828 residual components

The six signed `u=3` walls have numerator residues

```text
r in {1,13,15,27,29,41} mod 42.                          (8)
```

Their two adjacent anchor components have indices `20r-1,20r mod 840`, so

```text
C={19,20,259,260,299,300,539,540,579,580,819,820}.       (9)
```

At every wall, the inward component is strictly contained in a `3`-tooth and
the outward component in a `39`-tooth. The adjacent-side margin is
`1/11760`; the opposite-side margins are `547/11760` and `391/152880`.
Consequently these are strict open-tooth containments, not boundary touches.
Deleting exactly `C` leaves the labelled residual set

```text
R={0,...,839}\C,                       |R|=828.           (10)
```

## 4. The completed-chain quotient

For each `P` and `k in R`, run the deterministic recurrence from `(1)`.
Record its terminal status

```text
missing: a gap occurs before the right endpoint,
span:    one selected tooth covers the component,
renew:   at least two selected teeth cover the component.               (11)
```

Define `Q(P)` to retain the fixed collar set, the labelled status at every
`k in R`, and the complete physical `(speed,tooth-address)` chain whenever
the status is `span` or `renew`. It deliberately erases the selected prefix
when the status is `missing`.

Exact enumeration gives

```text
Q(241)=Q(255)=Q(761)=Q(1015),                            (12)

missing=546,                  span=276,
renew=6,                      completed=282.             (13)
```

All 282 completed chains agree literally, not merely after a role relabeling,
and none contains the varying speed `P`.

For reproducibility, the independent verifier serializes the common status
word as `k|M/S/R` plus LF and the completed map as
`k|status|speed@address,...` plus LF, over increasing `k`. The respective
byte-length/hash pairs are

```text
4860:  4d34447b9eca8c8a9302a0f799a56300b4b96135cd0c3a245a97960975f9a347,
10940: 64704c712a0e3e70ca0ecb3264834b3f61c4d473ddc20ea8c44ccc5c2c616d11. (14)
```

These hashes belong to those declared serializations; no digest is called
representation-independent.

## 5. The first-exit consumer does not factor

For a missing trace, let its exit record be

```text
(k,k mod h,epsilon,first uncovered x,clearance,binding-speed set),         (15)
```

where `epsilon=0,1` is the half-period copy bit (and equals the inherited
selected-tooth owner bit in these four traces). Let `E(P)` be the record at
the least missing `k in R`. For all four rows that component is

```text
k=23,                    I_23=[323/11760,67/2352].       (16)
```

The exact prefixes and exits are

```text
P=241:   945@26 -> 3371@93 -> EXIT,       x=1303/47194;
P=255:   255@7  -> 5051@140 -> EXIT,       x=1961/70714;
P=761:   945@26 -> 761@21  -> EXIT,        x=295/10654;
P=1015:  945@26 -> 1015@28 -> EXIT,        x=393/14210.   (17)
```

The four exits are pairwise distinct. Each has clearance exactly `1/14`,
with its final speed the unique binding speed; because the corresponding
tooth is open at its endpoint, the point is safe. The independent wall-cell
audit finds no other active tooth there.

Therefore `E(241),E(255),E(761),E(1015)` are pairwise distinct while `(12)`
holds. There is no function `f` on this declared four-source quotient such
that

```text
E=f o Q.                                                  (18)
```

This proves nonfactorization through the declared completed-chain quotient.

## 6. The role packet loses the metric tooth

The pair `P=761,1015` has the same partial role word `C1 -> P -> EXIT`, the
same first frontier, singleton right-wall tie sets at both choices, owner bit
zero, and the same missing continuation status. Their numerical addresses
`21,28` and exits nevertheless differ.

Precisely, let `S(P)` retain at each selected step of the first missing trace
the tuple

```text
(k,k mod h,copy bit,step index; selected-tooth owner parity;
 role,orientation; frontier; ordered tie-role tuple;
 next role/EXIT; terminal status),                                  (19)
```

where the owner parity is `(2n-(copy bit)w) mod 2`. Direct evaluation gives

```text
S(761)=S(1015),                    E(761)!=E(1015).       (20)
```

Thus `E` also does not factor through `S` on this two-source set. This local
packet forgets the whole physical tooth `(P,n)`, not merely its address.

The first common frontier is

```text
x0=b(945,26)=73/2646=1/36-1/5292.                       (21)
```

Write the varying odd speed as `P=36n+r`. Literal activity of its address-`n`
tooth at `x0` is

```text
|P*x0-n|<1/14
 iff |73P-2646n|<189
 iff |18n-73r|<189
 iff |P-147r|<378.                                      (22)
```

When active, its outgoing right wall is

```text
b(P,n)=(14n+1)/(14P)=1/36+(18-7r)/(252P).               (23)
```

For `(P,n,r)=(761,21,5),(1015,28,7)`,

```text
21*1015-28*761=7,
(14*21+1)1015-(14*28+1)761=352,                         (24)

b(761,21)-b(1015,28)=352/(14*761*1015)>0.               (25)
```

Thus the declared role packet does not determine the metric consumer. The
physical pair `(P,n)` is one concrete repair for these traces: `n` is the
Euclidean quotient in `P=36n+r`, while the remainder controls the wall
correction in `(23)`. Neither `n` alone nor `(P,n)` is claimed minimal or
sufficient for a general successor or global LRC certificate.

## 7. Audit and scope

The primary compares two exact implementations on all `4*828=3312` retained
instances, with a separate exact audit of the twelve collars. The import-free
referee rebuilds the open intervals, selection order, named gates, collars,
statuses, chains, exits, and arithmetic identities, and independently audits
every tooth wall and intervening cell on the same 3,312 instances. Normal,
optimized, hash-seeded, and frozen outputs agree byte-for-byte.

The proved universe consists only of the four rows in `(2)`. The theorem does
not classify the other 828-body rows, restrict an arbitrary hypothetical
counterexample, show that all consumers need the full missing prefix, prove
the displayed address packet sufficient, lower a residual rank, or advance
the logical LRC(14) bound. It identifies two precise unsafe quotients and a
missing metric sidecar exposed by one precise consumer.

Reproduce from the repository root:

```text
python3 -B 04-computation/lrc14_828_completed_chain_first_exit_nonfactorization_thm4363.py
python3 -B -O 04-computation/lrc14_828_completed_chain_first_exit_nonfactorization_thm4363.py
python3 -B 04-computation/lrc14_828_completed_chain_first_exit_nonfactorization_independent_referee_thm4363.py
python3 -B -O 04-computation/lrc14_828_completed_chain_first_exit_nonfactorization_independent_referee_thm4363.py
```
