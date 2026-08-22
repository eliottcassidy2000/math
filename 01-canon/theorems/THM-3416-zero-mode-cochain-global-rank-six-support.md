---
id: THM-3416
title: "Zero-mode-cochain global rank-six support"
status: >
  PROVED all-q global rank-six theorem, downstream of THM-3414's
  COMPUTER-ASSISTED PROVED fixed-zero atlas and THM-3415's rank-four/five
  divisor minimum, with an elementary all-q half-twist anchored-complement
  proof; VERIFIED-EXACT primitive hostile census for 2<=Q<=300 and exact
  harmonic/Farey/Fibonacci/Berggren transports; INDEPENDENTLY PROOF- AND
  ARITHMETIC-AUDITED, including a clean-room target-free Q<=1000 scan and five
  larger hostile moduli.  LRC(14) remains open.
source: codex-2026-08-15-major-frontiers
audit: half-twist capacity tail and four anchored complement gates; exact 17/29 lcm-493 mixed control; THM-3414 independently audited fixed-zero input; literal Q11/Q15/Q23/Q25 witness replay; rare-coordinate Q<=300 hostile census; normal/optimized transcripts byte-identical; independent proof/arithmetic audit PASS with clean-room target-free Q<=1000 and Q=1479,1972,2465,3451,4913 hostiles
referee_audit: independent formula, complement-colour, reflection/CRT, divisor-ancestry, atom, arithmetic, hash, routing, scope, and normal/-O replay audit CLEAN
depends_on:
  - THM-3414-fixed-zero-six-owner-base-classification
  - THM-3415-zero-mode-cochain-global-rank-five-support
related:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3408-fixed-zero-additive-order-duality-and-six-core-corridor
  - THM-3410-projective-cochain-wedge-ray-tree-tariff-and-residue-scalar-hubs
script: 04-computation/lrc_zero_mode_cochain_rank6_support_probe_20260815.py
output: 05-knowledge/results/lrc_zero_mode_cochain_rank6_support_probe_20260815.out
script_sha256: ed7672fd75b7c5ede23c0b7752e06849faa9c693d18e0639eee5b39f17d03a21
output_sha256: 3eb932643e76f2ac7836d8d435bf259183e2c3778216b0458bf4487b5894102a
semantic_sha256: 18abbdc82a40b0299d9dc59cd745d52a09e9706d3cf0cb0379ab5ed2f064df22
independent_script: 04-computation/lrc_zero_mode_cochain_rank6_support_thm3416.py
independent_output: 05-knowledge/results/lrc_zero_mode_cochain_rank6_support_thm3416.out
independent_script_sha256: c946cf9a66daa13da790cc9b1129993f9b5c1a8a2bdc6dec1bcc07b644989122
independent_output_sha256: 733a3ed02910348b87b561df7f1c79eef2dc431a8d702f2c0f807db33c1298cc
independent_semantic_sha256: 99892baf39b3d2b1b6a802bf21e0fe4164f155030d8ad051bc7ae26513b01ca3
hash_basis: LF-normalized bytes
---

# THM-3416 -- zero-mode-cochain global rank-six support

**PROVED (with THM-3414's COMPUTER-ASSISTED PROVED fixed-zero atlas as an
explicit dependency) + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

For `q>=2`, let `rho_ZMC(q)` be the minimum number of distinct positive
transverse owners whose strict danger blocks cover all `q` sheets at a common
THM-3398 mode centre and whose complete mode cochain vanishes.  Then

```text
rho_ZMC(q)=6
iff
(11|q or 15|q or 23|q or 25|q)
and 8,9,10,12 all do not divide q.                                 (1)
```

Together with THM-3415, the complete support through rank six is

```text
rho_ZMC(q)=4 iff 8|q or 9|q;

rho_ZMC(q)=5 iff (10|q or 12|q) and 8,9 do not divide q;

rho_ZMC(q)=6 iff (11|q or 15|q or 23|q or 25|q)
                    and 8,9,10,12 do not divide q.                 (2)
```

At half twist, let `r_(1/2)(Q)` be the minimum number of literal sheet blocks
covering `Z/QZ`, before imposing the augmented primitive gcd gate.  Then

```text
r_(1/2)(Q)<=6
iff some member of {8,9,10,11,12,15,23,25} divides Q.              (3)
```

Every primitive augmented half-twist cover by at most six owners therefore
has modulus divisible by one of the eight bases in `(3)`.  The converse is
not a primitive statement: pullback can destroy the breaker coordinate.  For
example, the `Q=33` pullback of the `Q=11` atom covers all sheets with residues
`(3,6,9,15,21,27)` but has augmented gcd three; the companion finds no
primitive rank-at-most-six cover at `Q=33`.

The new minimal positive atoms are

```text
Q=11, half twist: (1,2,3,5,7,9);
Q=15, zero twist: (1,2,3,4,5,7);
Q=15, half twist: (1,4,6,7,8,10);
Q=23, half twist: (1,4,5,7,9,11);
Q=25, half twist: (1,9,10,11,19,21).                              (4)
```

Each row satisfies the augmented gcd gate.  The `Q=11,23,25` rows are exact
sheet partitions, so OR equals XOR there.  Both `Q=15` rows instead have one
distinguished multiply covered sheet.  No augmented-mask XOR claim follows.

## 2. Divisor ancestry and the half-twist order model

THM-3415 records the proved divisor-minimum formula

```text
rho_ZMC(q)=min_(Q|q,Q>=2,epsilon in {0,1}) rho_epsilon^prim(Q).     (5)
```

For a half-twist owner residue `r`, put

```text
m=Q/gcd(Q,r),              r=(Q/m)s,              gcd(m,s)=1.      (6)
```

Dividing the phase word by `Q/m` reduces the strict predicate to modulus
`2m`.  Its sheet mask is `m`-periodic, and intersections of orders `a,m`
descend exactly to `lcm(a,m)`.  The maximum block size per `m` sheets is

```text
z(m)=1+2 floor((m-1)/14),
a(m)=2 ceil(floor((m-1)/7)/2),
h(m)=a(m)                         if m is even,
h(m)=max(a(m),z(m))               if m is odd.                     (7)
```

Both branches satisfy

```text
h(m)<=(m+6)/7.                                                    (8)
```

Consequently `h(m)/m<1/6` for every `m>=37`; direct evaluation below that
cutoff gives the all-order list

```text
6h(m)>=m iff
m in {3,5,8,9,10,11,12,15,17,22,23,24,29,36}.                    (9)
```

Assume for contradiction that `Q` avoids all eight bases in `(3)` but has a
cover by at most six blocks.  Total block mass forces at least one selected
block of density greater than `1/6`.  Removing from `(9)` all orders already
carrying a forbidden base leaves exactly

```text
E={3,5,17,29}.                                                     (10)
```

An actual order-17 or order-29 block above density `1/6` is necessarily a
maximal `z`-branch block, of sizes `3/17` and `5/29` respectively.

## 3. The four anchored-complement gates

For a maximal order-`a` anchor in `(10)`, let `C_a` be its complement.  Five
remaining slots must contribute on average

```text
theta_3 = 2/15,
theta_5 = 4/25,
theta_17=14/85,
theta_29=24/145                                                   (11)
```

of all `Q` sheets to `C_a`.  We eliminate the possible anchors in the order
`3,5,17,29`.

### Order three

The order-three block is one residue class modulo three, so `|C_3|=2Q/3`.
If `3` does not divide a companion order `m`, CRT leaves exactly two thirds
of that block in `C_3`.  The all-order density classification

```text
5h(m)>=m iff m in {3,5,8,9,10,15}                                 (12)
```

shows that the contribution is below `2Q/15` unless `m=5` or a forbidden
base divides `m`.  But an order-three anchor and an order-five companion
would force `15|Q`, also forbidden.

If `3|m`, the dangerous phase words run through residue classes modulo three
in a balanced interval.  At most `ceil(2h(m)/3)` of the `m` phase positions
avoid the anchor.  From `(8)`, for `m>=33`,

```text
ceil(2h(m)/3)/m <= (2m+26)/(21m) < 2/15.                           (13)
```

The only target-free multiples of three below 33 are `3,6,21`; their exact
complement contributions are `0,0,2/21`.  Thus every companion contributes
strictly less than `2Q/15`, and five cannot cover `C_3`.

### Order five

Now suppose no order-three block occurs.  The maximal order-five block is a
single residue class modulo five, with complement `4Q/5`.  By `(8)`, orders
above 50 have density below `4/25`.  Exact evaluation through 50 shows that,
after deleting forbidden bases and order three, only

```text
m in {5,17,29,31,37,43}                                           (14)
```

can even meet the required total density.  Exact descent to `lcm(5,m)` gives

```text
m:                    5,     17,     29,    31,     37,      43
max |B_m meet C_5|/Q: 0,   12/85,   4/29,  4/31,  24/185,  28/215,             (15)
```

all strictly below `4/25`.  Five companions again fall short.

### Order seventeen

Suppose neither order three nor order five occurs and select an order-17
anchor if one is present.  Its complement has size `14Q/17`.  The analytic
tail at the quota `14/85`, followed by exact evaluation below 40, leaves only
orders `17,29` after the inherited exclusions.  Their exact contributions are

```text
order 17 into C_17: 2/17,
order 29 into C_17: 70/493,                                       (16)
```

both strictly below `14/85`.  Equation `(16)` is the decisive mixed
`17/29` control: it is computed on `lcm(17,29)=493`, not inferred from raw
density.  Five companions cannot cover `C_17`.

### Order twenty-nine

Finally suppose only order 29 remains from `(10)`.  Its complement has size
`24Q/29`.  At quota `24/145`, the analytic tail begins at order 38; after the
same exclusions only order 29 survives the finite critical list.  Two maximal
order-29 blocks overlap enough that a second contributes at most

```text
4/29 < 24/145                                                       (17)
```

to the first complement.  This last contradiction proves the necessary
direction of `(3)`.

### Independent reflection/CRT closure of the 17/29 tail

Concurrent work found a second all-order mechanism for the last two anchors.
After orders three and five and all target-bearing orders are removed, every
order other than 17 and 29 has density at most

```text
7/43.                                                              (17a)
```

The bound follows from `(8)` for `m>=43`; exact evaluation below 43 has
target-free maximum `6/37`, while order 43 attains `7/43`.

Every half-twist mask is invariant under the reflection

```text
sigma(ell)=-1-ell.                                                  (17b)
```

On `Z/17Z` and `Z/29Z`, this reflection has one fixed point and all other
orbits are pairs.  Thus `a>=1` order-17 blocks cover at most `1+2a` quotient
points, and `b>=1` order-29 blocks cover at most `1+4b`.  Their missed
fractions are at least

```text
A_a=(16-2a)/17,          B_b=(28-4b)/29.                           (17c)
```

When both orders occur, `493|Q` and the equal-fibre CRT map to
`Z/17Z x Z/29Z` makes the jointly missed set a product cylinder of density at
least `A_a B_b`.  Put `c=6-a-b`; the remaining blocks cover at most `7c/43`.
The three exact gaps are

```text
b=0:       A_a-7(6-a)/43       =(33a-26)/731      >=7/731,
a=0:       B_b-7(6-b)/43       =(31b-14)/1247     >=17/1247,
a,b>=1:    A_a B_b-7c/43
             =(344ab+1043a+699b-1442)/21199       >=644/21199.     (17d)
```

All are positive.  This independently removes every pure or mixed 17/29
profile, with the product cylinder supplying the coordinate that scalar
density forgets.  It agrees with the anchored `70/493` gate but explains the
same obstruction through symmetry and CRT rather than pairwise complement
pricing.  The independently written `lrc_zero_mode_cochain_rank6_support_thm3416.py`
freezes `(17a)--(17d)` and a separate primitive census through `Q=200`.

## 4. Positive half-twist bases and the primitive boundary

The complete minimal witness bank is

| `Q` | residues | quotient-order profile |
|---:|---|---|
| 8 | `(1,3,5,7)` | `8^4` |
| 9 | `(1,5,6,7)` | `(3,9^3)` |
| 10 | `(1,3,4,7,9)` | `(5,10^4)` |
| 11 | `(1,2,3,5,7,9)` | `11^6` |
| 12 | `(1,5,7,8,11)` | `(3,12^4)` |
| 15 | `(1,4,6,7,8,10)` | `(3,5,15^4)` |
| 23 | `(1,4,5,7,9,11)` | `23^6` |
| 25 | `(1,9,10,11,19,21)` | `(5,25^5)` |

Every row covers literally and has augmented gcd one.  Scaling a row by
`Q/a` pulls the literal sheet cover from base `a` back to every multiple `Q`
of `a`, proving the reverse direction of `(3)`.  Such scaling need not retain
augmented gcd one; this is why `(3)` is deliberately a literal-cover theorem
and only its necessary direction is used for primitive quotients.

The new partition atoms have block-size profiles

```text
Q=11: (2,1,2,2,2,2),
Q=23: (4,3,4,4,4,4),
Q=25: (4,4,5,4,4,4).                                               (18)
```

Each is one anchor plus five disjoint petals.  The `Q=15` zero-twist row has
fourteen sheets of multiplicity one and its fixed sheet of multiplicity six;
the half-twist row has multiplicities `1^14,4^1`.

## 5. Zero twist and proof of the global classification

THM-3414 proves, for every `Q>=2`,

```text
a fixed-zero cover by at most six owners exists
iff one of {15,16,18,20,24} divides Q.                              (19)
```

This is a COMPUTER-ASSISTED PROVED all-`Q` finite-profile atlas, independently
audited; it is not a bounded-modulus scan.  Each base in `(19)` contains one
of the eight half-twist/global bases:

```text
15;       16 -> 8;       18 -> 9;       20 -> 10;       24 -> 12. (20)
```

Hence any primitive cover by at most six owners, at either twist, has quotient
divisible by one of the eight bases in `(3)`.

The four rows in `(4)` give rank at most six on every multiple of
`11,15,23,25`.  If none of `8,9,10,12` divides `q`, the rank-four and rank-five
classifications `(2)` exclude every lower rank, proving the reverse direction
of `(1)`.

Conversely, if `rho_ZMC(q)=6`, the divisor minimum `(5)` supplies a primitive
quotient `Q|q` covered by at most six owners.  Sections 3 and `(19)--(20)`
force one of the eight bases to divide `Q` and hence `q`.  A lower base
`8,9,10,12` would make the global rank at most five, contradiction.  Therefore
one of `11,15,23,25` divides `q`, and all lower bases are absent.  This proves
`(1)`.

## 6. The exceptional-order generalized tournament

There is one lawful tournament-like object here, but it is a proof sidecar,
not the cover carrier.  On vertices `E={3,5,17,29}`, define the observable

```text
c(a,b)=maximum fraction of an order-b block lying in C_a.          (21)
```

The exact matrix, with rows indexed by anchors, is

| `a\b` | 3 | 5 | 17 | 29 |
|---:|---:|---:|---:|---:|
| 3 | `0` | `2/15` | `2/17` | `10/87` |
| 5 | `4/15` | `0` | `12/85` | `4/29` |
| 17 | `14/51` | `14/85` | `2/17` | `70/493` |
| 29 | `8/29` | `24/145` | `72/493` | `4/29` |

Orient `a -> b` when `c(a,b)<theta_a`.  The strict arcs are

```text
3 -> 17,29;
5 -> 17,29;
17 <-> 29.                                                         (22)
```

The pair `{3,5}` is a missing edge at equality, but its lcm is the forbidden
base 15.  The pair `{17,29}` is genuinely bidirected.  Thus the intrinsic
sidecar is exactly a four-vertex generalized tournament with one missing
target-killed edge and one both-way edge.  Choosing the elimination order
`3,5,17,29` turns it into a proof scheduler.  It forgets residue orientation,
prime breakers, and literal sheet union, so it is not equivalent to a cover.
The reflection/CRT proof restores a different sidecar for the bidirected pair:
the product of the two missed cylinders.  Neither representation subsumes the
other; together they explain both the pairwise obstruction and its all-mixture
stability.

The actual positive objects remain six-edge pointed partition clutters.  This
distinction explains why tournament language is useful for the obstruction
grammar but not for reconstructing the atoms.

## 7. Harmonic, Farey, Fibonacci, and Berggren transports

The exact rank-six support is periodic modulo

```text
M=lcm(8,9,10,11,12,15,23,25)=455400.                              (23)
```

Exactly `55,000` residues are accepted, so

```text
#{q<=N:rho_ZMC(q)=6}=(25/207)N+O(1),
sum_(q<=N,rho_ZMC(q)=6) 1/q=(25/207)log N+O(1).                    (24)
```

The disjoint support through rank six has density and harmonic coefficient

```text
2/9 + 4/45 + 25/207 = 149/345.                                   (25)
```

For reduced fractions, retain `(numerator,denominator) mod M`.  The exact
primitive state count and accepting count are

```text
J_2(M)=131,383,296,000,
N_6(M)=16,682,793,600,                                             (26)
```

so the limiting Farey/Stern--Brocot denominator proportion is

```text
N_6(M)/J_2(M)=157981/1244160.                                     (27)
```

The child maps act on this finite residue-pair state space; no ordering of
fractions is lost, but only denominator acceptance is read out.

The Fibonacci sequence modulo `M` has period 1200.  For `n>=3`,

```text
rho_ZMC(F_n)=6
iff (10|n or 25|n) and 6 does not divide n and 15 does not divide n. (28)
```

Equivalently the accepted indices modulo 150 are

```text
10,20,25,40,50,70,80,100,110,125,130,140,                         (29)
```

of density `2/25` among indices.

For the Berggren `U`-spine label

```text
Q_n=4n^2+12n+11=2c_n+1,                                           (30)
```

the candidates 15,23,25 never divide `Q_n`, while `11|Q_n` exactly at
`n=0,8 mod 11`; the lower odd obstruction `9|Q_n` occurs at `n=1,5 mod 9`.
Therefore

```text
rho_ZMC(Q_n)=6
iff n mod 11 is 0 or 8 and n mod 9 is neither 1 nor 5.             (31)
```

This gives 14 classes modulo 99, of density `14/99`.  On the full ternary
Berggren tree every label `2c+1` is odd, so exact rank six is the finite-state
condition

```text
(11 or 15 or 23 or 25) divides 2c+1,          9 does not.          (32)
```

Tracking the Pythagorean triple modulo `lcm(9,11,15,23,25)=56925` under the
three Berggren matrices gives a literal ternary automaton.  The companion
counts accepted nodes through depth ten as

```text
1,0,2,5,10,41,133,378,1210,3519,10634.                            (33)
```

This is a finite exact prefix, not a claimed limiting tree density.

## 8. Exact companion and scope

The standard-library companion independently constructs augmented banks and
runs a rare-uncovered-coordinate branch-and-bound solver for both twists and
every `2<=Q<=300`.  It examines

```text
133,058 raw types,
65,604 unique augmented types,
65,063 maximal types,
734,163 memoized states,
751,364 branches.                                                  (34)
```

It finds 64 primitive rank-six twist instances and exactly 35 global rank-six
degrees through 300, reconstructs every finite critical-order and complement
gate, verifies all atoms and transports, and freezes event and semantic
digests

```text
354654b682bb9f9796e11e6f67cc4511bcb4d403a5ba3604e2471abbe14e0706,
18abbdc82a40b0299d9dc59cd745d52a09e9706d3cf0cb0379ab5ed2f064df22.
```

Normal and optimized Python outputs are byte-identical.  The finite primitive
census is not promoted beyond `Q=300`; the all-`q` global theorem instead uses
the analytic tails, exact order descents, divisor ancestry, and THM-3414.

The separate reflection/CRT companion has script/output/semantic hashes

```text
c946cf9a66daa13da790cc9b1129993f9b5c1a8a2bdc6dec1bcc07b644989122,
733a3ed02910348b87b561df7f1c79eef2dc431a8d702f2c0f807db33c1298cc,
99892baf39b3d2b1b6a802bf21e0fe4164f155030d8ad051bc7ae26513b01ca3.
```

It independently freezes the all-order order-three/order-five classifiers,
reflection union maxima, CRT deficits, positive atoms, arithmetic transports,
and a primitive cap-six census through `Q=200`.

An independent clean-room audit reconstructed the capacity and all four
anchored-complement gates, searched every target-free half-twist modulus
through `Q=1000`, and checked the larger hostiles
`1479,1972,2465,3451,4913` together with their proper divisors.  It also
replayed direct factor-997 pullbacks and THM-3414's fixed-zero transcript.
No counterexample was found.

Nothing here constructs a physical current, embeds an atom into an LRC(14)
row, or decrements the live ledger.  LRC(14) remains open.

**QED.**
