---
id: THM-4379
title: "LRC14 sixteen-step labeled-binder observability and phase decoder"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4363/4365/4367 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. On the fixed h=420 odd tail, the
  declared complete labeled first-exit record has sharp global horizon 16
  shifts/17 observations. Its closed marked phase interval is explicit, and
  it resolves every strict-active metric fibre at horizon zero. No safe-time,
  seam-entry, ledger, or LRC(14) consequence is asserted.
source: root + lrc_binder_observer + clean-room referee / next-sharp session, 2026-09-03
depends_on:
  - THM-4363-lrc14-828-body-completed-chain-first-exit-nonfactorization
  - THM-4365-lrc14-cofinite-828-quotient-fibre-and-centered-residue-exit-law
  - THM-4367-lrc14-active-first-exit-scale-collision-classification
related:
  - THM-4374-lrc14-seventeen-step-metric-exit-observability-and-shift-congruence-rigidity
mistake_firewall:
  - MISTAKE-222
primary_script: 04-computation/lrc14_sixteen_step_labeled_binder_observability_thm4379.py
primary_output: 05-knowledge/results/lrc14_sixteen_step_labeled_binder_observability_thm4379.out
primary_script_sha256: cd876cca0b1563fcb08a6c09337677e700b339f39553461f978dada436dd04a2
primary_output_sha256: f4d8ac1206006dfd5ae889556c697e567a8e68f18a69eac6738c82df659f5845
independent_referee_script: 04-computation/lrc14_sixteen_step_labeled_binder_observability_independent_referee_thm4379.py
independent_referee_output: 05-knowledge/results/lrc14_sixteen_step_labeled_binder_observability_independent_referee_thm4379.out
independent_referee_script_sha256: 6bcd40128da170fb6e935653683f92b9b648c75d3f5b003068a18e510061f4fd
independent_referee_output_sha256: 3712cbc3921618ab9ced6e4e639842eaa840705a0dc45bb2ace1499042e9850d
hash_basis: raw LF bytes
audit: >
  PASS WITH DECLARED-RECORD, FIXED-TAIL, U_16, FULL ROLE-PROJECTION, SIDECAR,
  AND BOUNDARY-COUNT WORDING REPAIRS. The 4,374,628-check primary and
  5,109,296-check independent integer-rational referee verify the exact record
  law, sharp decoder horizon, partitions, boundary ties, active-fibre collapse,
  and role-normalized hostile in all audited modes.
---

# THM-4379 -- LRC14 sixteen-step labeled-binder observability and phase decoder

**PROVED ELEMENTARY RELATIVE TO THM-4363/4365/4367 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. ON THE FIXED `h=420` ODD TAIL, THE DECLARED LABELED
FIRST-EXIT RECORD HAS SHARP HORIZON 16 SHIFTS / 17 OBSERVATIONS. EVERY ROW IS
SAFE AND LRC(14) REMAINS OPEN.**

## Scope and inheritance

The source is only the fixed `h=420` family of
[`THM-4365`](../../01-canon/theorems/THM-4365-lrc14-cofinite-828-quotient-fibre-and-centered-residue-exit-law.md),
on odd `P>=11019`. The active metric-fibre normal form is inherited from
[`THM-4367`](../../01-canon/theorems/THM-4367-lrc14-active-first-exit-scale-collision-classification.md),
and the comparison target is the sharp 17-shift metric theorem
[`THM-4374`](../../01-canon/theorems/THM-4374-lrc14-seventeen-step-metric-exit-observability-and-shift-congruence-rigidity.md).
The exit-record type and owner parity are those of
[`THM-4363`](../../01-canon/theorems/THM-4363-lrc14-828-body-completed-chain-first-exit-nonfactorization.md).

Inheritance pass:

- closest mechanism: THM-4365's centered-residue exit law and THM-4374's
  half-phase rotation;
- canonical hostile: THM-4367's equal metric exits erase physical scale;
- corrected near miss: boundary equality is safe but has two physical
  binders, not the singleton fixed binder;
- least-used sidecar: wall side and physical Euclidean tooth address.

Live board:

```text
centered phase | physical scale | Euclidean address | open boundary
binder set | wall side | owner parity | finite output word
```

The Anchor was the sharp physical-record horizon. The Niche was the owner and
boundary-side content. The Wildcard was the role-normalized loss hostile in
Section 6.

The applicable `META-PATTERNS` card was **Controlled forgetting and unlabeled
quotients require a sidecar**: the fibre product collision is precisely scale
erasure, and the physical tooth (or an injective address for it) is the needed
sidecar. No new method-card promotion is proposed.

## 1. The declared complete local record

Put

```text
A=3371,       S=1303,       M=14A=47194,
N=M/2=23597,  X0=S/M.
```

For an odd tail parameter `q`, define

```text
S q = M n(q) + rho(q),       -M/2 < rho(q) <= M/2,
z(q) = (A-rho(q))/2 mod N,   0 <= z(q) < N.              (1)
```

The half-phase evolves exactly by

```text
z(q+2)=z(q)-S mod N.                                      (2)
```

The fullest first-exit record used here is

```text
R(q)=(k,k mod h,copy bit,status;
      selected physical chain with addresses and owner parities;
      selected frontiers; first uncovered point; clearance;
      binding physical teeth with address, wall side, owner parity).       (3)
```

No centered residue is artificially added to the emitted record. The phase
`z` is a proof coordinate. Tooth notation below is
`speed@address[side,owner]`. On component `k=23`, the copy bit is zero, so the
THM-4363 owner formula

```text
(2 address - (copy bit) speed) mod 2
```

is identically zero. Thus owner parity itself supplies no separation.

The fixed selected prefix is

```text
945@26 -> 3371@93,
frontiers=(323/11760,73/2646),
X0=1303/47194.                                           (4)
```

THM-4365 yields the following complete piecewise record law.

```text
1 <= z <= 3370 (strict active):
  chain   = 945@26 -> 3371@93 -> q@n(q) -> EXIT
  exit    = (14n(q)+1)/(14q)
  binders = {q@n(q)[R,0]}

z=0 (rho=+3371, positive boundary):
  chain   = 945@26 -> 3371@93 -> EXIT
  exit    = X0
  binders = {3371@93[R,0], q@n(q)[R,0]}

z=3371 (rho=-3371, negative boundary):
  chain   = 945@26 -> 3371@93 -> EXIT
  exit    = X0
  binders = {3371@93[R,0], q@n(q)[L,0]}

z notin {0,...,3371} (strict inactive):
  chain   = 945@26 -> 3371@93 -> EXIT
  exit    = X0
  binders = {3371@93[R,0]}.                              (5)
```

Every case has status `missing` and clearance exactly `1/14`. At a strict
active point the final varying tooth is selected and the exit is its right
wall. At either boundary the varying open tooth is not selected, but it still
binds metrically at the exit. The sign of `rho` is visible as its right/left
wall side.

For `q_t=P+2t`, the exact output word is therefore

```text
B_H(P)=(R(q_0),R(q_1),...,R(q_H)),
rho_t=<rho_0+2St>_M,
n_t=(S(P+2t)-rho_t)/M,
z_t=z_0-tS mod N.                                       (6)
```

This retains both scale (`q_t`) and phase (`z_t` and the boundary side in the
proof/record respectively), rather than replacing the physical tooth by an
abstract ordinal.

## 2. Sharp global-on-the-fixed-tail horizon

Let

```text
J={0,1,...,3371} subset Z/NZ.                            (7)
```

By (5), the numeric varying speed occurs in `R(q)` if and only if `z(q)` lies
in `J`. Consequently the first occurrence at time `t` exposes the literal
binder `q_t=P+2t`, and recovers

```text
P=q_t-2t.                                                (8)
```

The preimages of `J` at times `0,...,H` are `J+tS`. Because `S=1303` is
shorter than the interval `J`, their integer lifts overlap consecutively:

```text
union_{t=0}^H (J+tS) = [0,3371+1303H]     for 0<=H<=15.  (9)
```

At `H=15` this ends at `22916`, leaving exactly the 680 phases

```text
{22917,...,23596}.                                      (10)
```

At `H=16` the lift reaches `24219>=23596`, so its projection covers all of
`Z/NZ`. Every parameter therefore exposes its numeric varying binder by time
16, and (8) makes `B_16` injective.

This is sharp. Take

```text
P1=11123,       P2=P1+M=58317,       z(P1)=z(P2)=22927. (11)
```

For `t=0,...,15`, both phases stay strictly outside `J`; both words are the
same repeated fixed record with binder `{3371@93[R,0]}`. At time 16 their
phase is `2079`, and the records split as

```text
P1+32=11155:   exit=4313/156170,   binder={11155@308[R,0]},
P2+32=58349:   exit=22555/816886,  binder={58349@1611[R,0]}. (12)
```

Hence the minimal global observability horizon is exactly

```text
16 shifts = 17 observations.                            (13)
```

This is one shift shorter than THM-4374's metric horizon `17` (18
observations).

## 3. Exact finite-horizon equivalence and phase filtration

For `0<=H<=15`, put

```text
U_H={3372+1303H,...,23596},
|U_H|=20225-1303H,                 U_16=empty.            (14)
```

Then the full word equivalence is exactly

```text
B_H(P)=B_H(P')
iff P=P' or both z(P),z(P') lie in U_H.                  (15)
```

Indeed, parameters in `U_H` emit only the common fixed record. Outside
`U_H`, a first varying binder occurs in the word; equality at that position
forces the same numeric binder and hence the same starting parameter by (8).
At `H=16`, equality of words is equality of parameters.

Thus, on the whole infinite tail, for every `H<=15` there is exactly one
nonsingleton word fibre, and it is infinite; every other word fibre is a
singleton. At `H=16` all fibres are singletons. On the finite half-phase
circle the unresolved bucket sizes are

```text
H:     0     1     2     3     4     5     6     7     8
|U|: 20225 18922 17619 16316 15013 13710 12407 11104  9801

H:     9    10    11    12    13    14   15  16
|U|:  8498  7195  5892  4589  3286  1983  680   0.     (16)
```

The exact first marked-hit distribution is

```text
t=0: 3372;       t=1,...,15: 1303 each;       t=16: 680. (17)
```

Compared with THM-4374's strict-metric distribution `3370,1303x15,682`,
closing the target adds the two time-zero boundary phases and gives the
net counts in (17); the intervening orbit slabs reshuffle.

## 4. Active metric-fibre filtration

On a strict-active THM-4367 metric fibre, write

```text
q=bg,       A-rho=ag,       a+Sb=A kappa,
n=(g kappa-1)/14.                                      (18)
```

All scales in the fibre have the same metric exit `kappa/(14b)`, but their
full records bind at the distinct numeric teeth

```text
bg @ ((g kappa-1)/14).                                  (19)
```

Therefore the full record is already injective on every active metric fibre
at horizon zero. The sharp maximum fibre-size comparison is

| output | `H=0` | `H=1,...,15` | `H=16` | `H=17` |
|---|---:|---:|---:|---:|
| metric only (THM-4374) | 241 | 94 | 49 | 1 |
| full physical record | 1 | 1 | 1 | 1 |

The verifier exhausts all `6,106` structural types, `13,427` structural scale
points, and `281,073` scale pairs used in THM-4374's finite reduction. Every
metric fibre is constant metrically and singleton under the full record.

## 5. Hostile boundaries

The least tail representatives are exactly

```text
P=43823, rho=-3371, n=1210:
  {3371@93[R,0],43823@1210[L,0]};

P=50565, rho=+3371, n=1396:
  {3371@93[R,0],50565@1396[R,0]}.                        (20)
```

Both exits equal `1303/47194`; both varying teeth are inactive because teeth
are open. Nevertheless the clearance tie is physical and must remain in the
binding set. Collapsing either row to `{3371}` is the first failed implication
in the naive strict-active/inactive dichotomy.

## 6. Loss ledger and why the horizon improves

The decisive added coordinate is the numeric physical binder, not owner
parity:

- **source:** full record `R(q)`;
- **quotient:** role-normalize every varying-tooth occurrence in both the
  selected chain and binder set, erasing its numeric `q@n` data;
- **preserved:** metric, fixed/varying/tie type, wall side, owner parity;
- **destroyed:** physical scale and Euclidean address;
- **needed sidecar:** on a fixed strict-active metric fibre, the scale `g`
  (equivalently the literal speed, or `n` together with the fibre data); at a
  boundary, an address also needs the tie/wall-side data;
- **hostile:** THM-4374's pair `253031,258645` remains equal through horizon
  16 after this role normalization, and splits only at time 17.

The repaired primary and independent referee both apply this normalization
to the full declared record.  Thus the one-step gain is not available to an
unlabeled owner/tie word. The full record earns it because the first varying
binder is already a lossless
parameter address. The sharp same-phase pair (11) separately shows that phase
alone does not restore scale.

As in THM-4374, equality is consequently the only output-compatible forward
`+2` congruence for the full record: iterating any such congruence through 16
shifts gives equal `B_16` words and therefore equal parameters.

## 7. Exact verification

The primary uses exact `Fraction` geometry and compares the centered law
against direct open-tooth renewal and all thirteen distances on 47,194 rows.
Its role-normalization hostile projects every varying tooth in both the chain
and binder set.  The import-free clean-room referee uses normalized integer
rational pairs, a separately designed renewal, phase inverse/decoder, and a
constructive census of every active structural fibre shape.

```text
python -B 04-computation/lrc14_sixteen_step_labeled_binder_observability_thm4379.py
python -B 04-computation/lrc14_sixteen_step_labeled_binder_observability_independent_referee_thm4379.py
```

The primary passes `4,374,628` checks and the independent referee passes
`5,109,296`.  Normal, optimized, isolated, and fixed-hash-seed runs are
byte-identical to their respective frozen raw-LF outputs.

```text
primary script:    cd876cca0b1563fcb08a6c09337677e700b339f39553461f978dada436dd04a2
primary output:    f4d8ac1206006dfd5ae889556c697e567a8e68f18a69eac6738c82df659f5845
referee script:    6bcd40128da170fb6e935653683f92b9b648c75d3f5b003068a18e510061f4fd
referee output:    3712cbc3921618ab9ced6e4e639842eaa840705a0dc45bb2ace1499042e9850d
hash basis:        raw LF bytes
```

The independent design was frozen before the candidate implementation was
read.  Both routes verify the complete declared record projection, phase
cover, whole-tail equivalence classes, exact boundaries, active scale fibres,
and the role-normalized THM-4374 hostile.

## 8. Honest frontier

This theorem concerns only the first missing component of one
explicit common quotient fibre. It neither changes the safe time `1/14`, nor
enters the unresolved `2+12` seam, nor excludes a counterexample family, nor
decrements an LRC ledger. It shows that a rich physical observer recovers its
own varying speed after a sharp finite delay. **LRC(14) remains OPEN.**
