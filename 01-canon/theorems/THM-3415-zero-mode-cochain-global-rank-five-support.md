---
id: THM-3415
title: "Zero-mode-cochain global rank-five support"
status: >
  PROVED all-q global rank-five support theorem, downstream of THM-3408's
  COMPUTER-ASSISTED PROVED arithmetic cutoff, plus VERIFIED-EXACT Q<=500
  primitive hostile census and recurrence/Farey companion; independently
  proof- and arithmetic-audited.  The exact global rank is five iff (10|q or
  12|q), 8 does not divide q, and 9 does not divide q.  The primitive
  classification beyond Q=500 is not claimed.  No LRC(14) ledger decrement
  follows.
source: codex-2026-08-15-major-frontiers
audit: exact density and order-three-complement proof; THM-3408 additive-order dual gate; disjoint Q10/Q12 witness replay; rare-coordinate Q<=500 census; normal/optimized transcripts byte-identical; independent proof and arithmetic audits PASS, including a selective Q=555 hostile
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3408-fixed-zero-additive-order-duality-and-six-core-corridor
related:
  - THM-3414-fixed-zero-six-owner-base-classification
  - THM-3410-projective-cochain-wedge-ray-tree-tariff-and-residue-scalar-hubs
  - THM-3411-pairwise-independent-toggle-filter-and-sharp-high-minor-norm
verified_by:
  - 04-computation/lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.py
  - 05-knowledge/results/lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.out
  - 04-computation/lrc_zero_mode_cochain_rank5_support_probe_20260815.py
  - 05-knowledge/results/lrc_zero_mode_cochain_rank5_support_probe_20260815.out
script: 04-computation/lrc_zero_mode_cochain_rank5_support_probe_20260815.py
output: 05-knowledge/results/lrc_zero_mode_cochain_rank5_support_probe_20260815.out
script_sha256: b753d8a1f8fd2e8c636ab81815b1561db2a3c451feda5b7dacfef8b57c060306
output_sha256: 7528c4b508cd4a762f22642192c0209f6ba97d9987ad6d7c3157bd7ccfc01832
semantic_sha256: 8c4a9bc6fee718a2c6f6ec282b14bb7f7968aec40ec37b9c8d12a3e13e8528dc
hash_basis: LF-normalized bytes
---

# THM-3415 -- zero-mode-cochain global rank-five support

**PROVED (with the COMPUTER-ASSISTED PROVED THM-3408 cutoff as an explicit
dependency) + VERIFIED-EXACT + INDEPENDENTLY AUDITED; the primitive census
beyond `Q=500` is not claimed.**

## 1. Statement

For `q>=2`, let `rho_ZMC(q)` be the minimum number of distinct positive
transverse owners (`q` does not divide an owner) whose strict danger blocks
cover all `q` sheets at a common THM-3398 mode centre and whose complete mode
cochain vanishes.  If no such family exists, put `rho_ZMC(q)=infinity`.

Then

```text
rho_ZMC(q)=5
iff (10|q or 12|q) and 8 does not divide q and 9 does not divide q.   (1)
```

Together with the already proved rank-four classification,

```text
rho_ZMC(q)=4 iff 8|q or 9|q,                                      (2)
```

this determines the complete global support through rank five.

More locally, every primitive quotient carrying a cover by at most five
owners is divisible by at least one member of

```text
{8,9,10,12}.                                                       (3)
```

The converse to `(3)` is intentionally not stated as a primitive theorem.
The exact independent census through `Q=500` finds primitive rank five at

```text
zero twist: Q=16,18;
half twist: Q=10,12,
            every 8|Q with Q>8,
            every 9|Q with Q>9,                                   (4)
```

and nowhere else in that finite universe.  Statement `(1)`, unlike `(4)`,
is all-q.

The two global atoms are exact disjoint partitions:

```text
Q=10, residues (1,3,4,7,9):
{0,9}, {3,6}, {2,7}, {1,8}, {4,5};

Q=12, residues (1,5,7,8,11):
{0,11}, {2,9}, {3,8}, {1,4,7,10}, {5,6}.             (5)
```

Thus OR and XOR of the five **sheet masks** agree in `(5)`.  No XOR statement
is made for the augmented prime-breaker coordinates.  This is a partition
identity, not a tournament equivalence; prime-breaker coordinates remain
part of primitive realizability.

## 2. Primitive divisor reduction

THM-3405 and the divisor-ancestry reduction write an active family as

```text
U=dV, gcd(V)=1,  g=gcd(q,d),  Q=q/g.                               (6)
```

Every zero-cochain certificate is the `g`-fold pullback of a primitive
`Q`-sheet cover at twist `epsilon in {0,1}`.  For selected primitive residue
types `r_i`, gcd-one positive lifts are encoded by the augmented condition

```text
gcd(M_epsilon,r_1,...,r_s)=1,
M_0=Q, M_1=2Q,                                                      (7)
```

where every full cover has `s>=2`.  In particular primitive gcd one implies

```text
lcm_i(Q/gcd(Q,r_i))=Q.                                             (8)
```

Only this necessary direction is used below.  Formula `(6)` gives the exact
divisor minimum

```text
rho_ZMC(q)=min_(Q|q,Q>=2,epsilon) rho_epsilon^prim(Q).              (9)
```

The previous quotient-order proof gives the universal floor four and `(2)`.
That input is an unnumbered proved result, not a consequence of THM-3405
alone: see the
[divisor-ancestry proof](../../07-reflections/lrc-zero-mode-cochain-divisor-ancestry-and-boolean-realization-codex-20260815.md)
and its independently replayed
[exact companion](../../04-computation/lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.py).
It remains to prove `(3)` for a primitive five-cover.

## 3. Half twist forces the four bases

For an owner of quotient order `m`, let `h(m)` be the largest possible
half-twist block count on `m` quotient sheets.  Exactly

```text
z(m)=1+2 floor((m-1)/14),
a(m)=2 ceil(floor((m-1)/7)/2),
h(m)=a(m)                         if m is even,
h(m)=max(a(m),z(m))               if m is odd.                      (10)
```

Both `a(m)` and `z(m)` are at most `(m+6)/7`.  Hence orders at least 16
have density strictly below `1/5`, and orders at least 37 have density
strictly below `1/6`.  Direct evaluation below those analytic cutoffs gives

```text
5h(m)>=m iff m in {3,5,8,9,10,15},                                 (11)

6h(m)>=m iff
m in {3,5,8,9,10,11,12,15,17,22,23,24,29,36}.                      (12)
```

Suppose first that no selected owner has order three.  If an order-eight or
order-nine owner occurs, `(3)` already holds.  Otherwise `(11)` says every
block has density at most `1/5`, with equality only at orders `5,10,15`.
Five blocks can cover only if all five attain equality.  If order ten occurs,
then `10|Q`.  If it does not, `(8)` gives `Q|15`.  At `Q=5` all nonempty
blocks are the same singleton.  At `Q=15`, every size-three order-five or
order-fifteen block contains sheet seven, so five maximum blocks cannot be a
partition.  Both residual moduli fail.

For `r=(Q/m)s`, dividing numerator and modulus by `Q/m` reduces the strict
half-twist predicate to modulus `2m`; its sheet mask is `m`-periodic.
Therefore its intersection with an order-three mask descends exactly to
`lcm(3,m)`, rather than merely admitting a density estimate.

Now suppose an order-three block occurs.  Every nonempty order-three type has
the same sheet block; call its complement `C`.  It has size `2Q/3`.  Descent
to `lcm(3,m)` and exact evaluation of the finite list `(12)` show that the
largest fraction of all sheets contributed to `C` by an order-`m` block is

```text
greater than 1/6 only for m=9,
equal to   1/6 only for m=8,12,                                   (13)
less than  1/6 for every other m.                                 (14)
```

Outside `(12)`, `(14)` follows already from total block density.  If at least two
order-three types were selected, at most three other blocks contribute to
`C`, so `(13)--(14)` force an order eight or nine.  With exactly one
order-three type, four blocks must cover `C`.  An order-nine block again gives
`9|Q`; otherwise equality is necessary in all four terms, so their orders
belong to `{8,12}`.  The lcm with order three is then divisible by 12 (and by
8 if order eight occurs).  This proves `(3)` at half twist.

## 4. Zero twist reduces to one odd hostile and fails there

THM-3408 associates to each fixed-zero quotient order `m` its exact density
`alpha(m)` on every additive-order sheet stratum.  Call a modulus base-free
when none of

```text
{15,16,18,20,24}                                                    (15)
```

divides it.  This is the same five-base set that the later THM-3414 proves is
equivalent to the existence of a fixed-zero cover by at most six owners.
THM-3414 is therefore a stronger compatible route to the first divisibility
conclusion below, but it is not an independent audit: both arguments inherit
THM-3408's computer-assisted arithmetic cutoff.  The sharper five-owner
argument retained here explains the equality boundary.  THM-3408's complete
arithmetic cutoff implies that every
base-free order satisfies

```text
alpha(m)<=1/5,
alpha(m)=1/5 iff m in {22,33,44,50}.                                (16)
```

If a base-free `Q` had a zero-twist cover by at most five owners, the unit
sheet stratum and `(16)` would force exactly five owners, all at equality.
If order 50 occurs, move to the additive-order `Q/5` stratum.  Its image
order drops from 50 to 10, where `alpha(10)=0`; every other term is at most
`1/5`, so the total is at most `4/5`.  If order 50 does not occur, all five
orders lie in `{22,33,44}`.  On the `Q/11` stratum they drop to orders
`2,3,4`, all of alpha-density zero.  Both alternatives contradict a cover.
Therefore every zero-twist five-cover is divisible by a base in `(15)`.

Assume for contradiction that it nevertheless avoids all four target bases
in `(3)`.  Divisibility by `16,18,20,24` would force divisibility by
`8,9,10,12`, respectively, so only `15|Q` remains.  Since `10` does not
divide `Q`, the modulus is odd; since `9` does not divide it, `v_3(Q)=1`.

At zero twist every block contains sheet zero.  If no order-three owner is
selected, every allowed odd order has density at most `1/5`; even equality
cannot cover because five blocks have a common point.  If an order-three
block is selected, its complement again has size `2Q/3`.  Outside

```text
{3,5,15,17,29},                                                     (17)
```

total zero-block density is below `1/6`.  Exact descent to `lcm(3,m)` gives
the complement fractions for `(17)` as

```text
0, 2/15, 2/15, 2/17, 10/87,                                       (18)
```

all strictly below `1/6`.  Four remaining blocks cannot cover `2Q/3`; extra
order-three duplicates only worsen the bound.  This eliminates the last
hostile and proves `(3)` at zero twist.

## 5. Proof of the exact global classification

The two rows in `(5)` satisfy the augmented gcd gate `(7)`, so they are
primitive half-twist covers of ranks five at `Q=10,12`.  They cannot have
rank four by `(2)`.  Pulling them back along `(6)` gives rank at most five on
every multiple of 10 or 12.  If the multiple is divisible by neither 8 nor
9, `(2)` and the universal floor force exact rank five.  This proves the
reverse implication in `(1)`.

Conversely, suppose `rho_ZMC(q)=5`.  Formula `(9)` supplies a primitive
quotient `Q|q` with a cover by at most five owners.  By `(3)`, one of
`8,9,10,12` divides `Q`.  The first two alternatives would make `q` rank four
by `(2)`, contrary to the assumption.  Hence `10|q` or `12|q`, and the rank
four alternatives are absent.  This proves the forward implication and
therefore `(1)`.

For completeness, the infinite primitive breaker constructions behind the
finite pattern `(4)` are explicit.  For `Q=8k`, `k>=2`, use

```text
(1,k,3k,5k,7k),                                                     (19)
```

and for `Q=9k`, `k>=2`, use

```text
(1,k,5k,6k,7k).                                                     (20)
```

The last four residues pull back the rank-four sheet partition; residue one
breaks every prime divisor of `2Q`.  The zero-twist witnesses at `Q=16,18`
are `(1,3,5,7,8)` and `(1,5,6,7,9)`.  The companion proves only the absence
of further primitive rank-five cases through `Q=500`; no all-q primitive
converse is inferred from that census.

## 6. Harmonic, Farey, Fibonacci, and Berggren transports

The rank-five support consists of 32 residue classes modulo 360:

```text
10,12,20,30,50,60,70,84,100,110,130,132,140,150,156,170,
190,204,210,220,228,230,250,260,276,290,300,310,330,340,348,350.    (21)
```

Consequently

```text
#{q<=N:rho_ZMC(q)=5}=(4/45)N+O(1),
sum_(q<=N,rho_ZMC(q)=5) 1/q=(4/45)log N+O(1).                       (22)
```

The disjoint rank-four and rank-five supports together have natural and
harmonic density `2/9+4/45=14/45`.  Thus the subset-of-the-harmonic-series
view has an exact coefficient, not merely a metaphor.

For reduced fractions, retain the pair `(a,q) mod 360`.  There are `82,944`
primitive residue-pair states and exactly `4,128` have a denominator in
`(21)`, giving

```text
4128/82944=43/864.                                                   (23)
```

An elementary Mobius lattice-point count makes primitive residue states
equidistributed in the Farey triangle, so `43/864` is also the limiting
proportion of Farey fractions whose denominator has exact ZMC rank five.
The Calkin-Wilf or Stern--Brocot child maps act on these residue pairs, making
`(21)` a finite-state accepting condition on the reduced-fraction tree.

The Fibonacci sequence modulo 360 has period 120.  In one period rank four
occurs at indices divisible by six, while rank five occurs at

```text
n=15,45,75,105 mod 120.
```

Equivalently, for `n>=3`,

```text
rho_ZMC(F_n)=5 iff n=15 mod 30.                                    (24)
```

The periodic rank-four list formally contains index zero, but `F_0=0` is
outside the theorem's degree domain `q>=2`; no claim at `n=0` is intended.

Finally the Berggren `U`-spine label

```text
Q_n=4n^2+12n+11=2c_n+1                                             (25)
```

is odd, so it never lies in the rank-five support.  It has rank four exactly
for `n=1,5 mod 9`; all other spine labels have rank at least six or no cover.
More generally every plane label `2c+1` from a primitive Pythagorean triple
is odd, so the same rank-five exclusion holds on the full ternary Berggren
tree.  The special `U`-spine supplies the sharper rank-four residue law.

## 7. Carrier and scope

Each atom in `(5)` is one anchor block plus four disjoint petals.  That is the
precise source of a recurring “size four” shadow: after fixing the order-five
or order-three anchor, four obligations remain.  No intrinsic pairwise
observable or orientation is present, so the carrier is a pointed five-edge
partition clutter, not a tournament.  THM-3411's four-toggle filtration is a
method-level analogue—both retain a four-way high-order sidecar after a lower
quotient—but supplies no object-level map.

THM-3410 gives nonzero projective scalar fibres over related ternary and
five-colour physical partitions.  The present theorem sits on the zero
cochain fibre and classifies owner count; THM-3410 prices information lost
when the physical partition quotient forgets the cochain.  Neither theorem
embeds these covers into an LRC(14) row, supplies a physical current, or
decrements the live ledger.

## 8. Exact companion

The standard-library companion independently builds every augmented bank for
both twists and every `2<=Q<=500`.  Its rare-coordinate branch-and-bound
search examines

```text
371,760 raw types,
184,338 unique augmented types,
183,424 maximal types,
381,717 memoized states,
386,458 branches.                                                   (26)
```

It reproduces the rank-four atoms, the finite primitive pattern `(4)`, all 45
global rank-five degrees through 500, the analytic critical-order lists,
both XOR partitions, 115 explicit breaker lifts, the 360-state arithmetic
transports, and the two hostile zero-complement routes.  Normal and optimized
outputs are byte-identical.  Its event and semantic digests are

```text
28f6ed1fa9322355bbde1b302f27536ac9ceef62cdb7c532cde4e3dd10369798,
8c4a9bc6fee718a2c6f6ec282b14bb7f7968aec40ec37b9c8d12a3e13e8528dc.
```

Two independent audits pass: one reconstructed the proof and descent gates,
and one verified the arithmetic, hashes, and normal/optimized replay.  As a
selective beyond-census hostile, both twists at `Q=555` were checked and have
no primitive cover of rank at most five.  This spot check is evidence only,
not an extension of the finite converse in `(4)`.

**QED.**
