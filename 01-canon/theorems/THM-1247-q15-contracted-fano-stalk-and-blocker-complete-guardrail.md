---
id: THM-1247
title: THE q=15 SUNFLOWER IS A CONTRACTED FANO STALK, NOT A CLOSING OBSTRUCTION — Fano incidence, the Kakeya cut, and complete third-blocking coexist with a lonely phase
status: PROVED (unique sign-address contraction; two invariant Fano planes with common degenerate image; CRT independence of speed chi7; q15/cut/carrier-spoke/all-sum-beat blocker-complete lonely guardrail; dependency-free exact referee; sorry-free Lean finite core)
source: codex-2026-07-19-S78 continuation with sunflower-fano-probe agent
depends_on: [THM-1156, THM-1238, THM-1240, THM-1241, THM-1242]
related: [THM-1244, HYP-7870]
script: 04-computation/lrc14_q15_contracted_fano_guardrail_thm1247.py
output: 05-knowledge/results/lrc14_q15_contracted_fano_guardrail_thm1247.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCQ15ContractedFanoGuardrail.lean
script_sha256: e446f579f11378514f32449469a5496f0d4a67268ae71c5b73f704001da4c19b
output_sha256: 90a76fcfe439bb2ee6925c203672813ca5998ecb5be62623f0f13270dd866068
formalization_sha256: 4f4756f1a6d89009a6261fcbfcae6091105c45da6d1e8bdcacde35d8bc870195
---

# THM-1247 — q=15 contracted Fano stalk

## 1. Exact seven-point sign-address refinement

Let

```text
P_r={r,15-r},                   r=1,...,7,             (1)
```

be the seven nonzero points of
`(Z/15Z minus {0})/{+/-1}`.  For a speed `s` with nontrivial quotient period,
its strict `q=15` danger mask is

```text
M_s={p:sp==0,+1,-1 mod 15}.                            (2)
```

There are exactly six distinct masks:

```text
{0} union P_r,             r=1,2,4,7,
{0} union P_5,
{0} union P_3 union P_6.                               (3)
```

Thus THM-1242's six-mask full clock has the canonical seven-address
refinement with the unique contraction

```text
P_3~P_6.                                               (4)
```

Writing its six petals as

```text
A=P1, B=P2, C=P3 union P6, D=P4, E=P5, F=P7,          (5)
```

the cyclic nonzero word is

```text
A B C D E C F F C E D C B A.                          (6)
```

It is palindromic, so its signed consecutive-transition current vanishes.

## 2. The two invariant Fano planes

Multiplication by two acts on the sign addresses as

```text
T=(1 2 4 7)(3 6)(5).                                  (7)
```

Exactly two labelled Fano planes are `T`-invariant:

```text
F0={123,145,167,246,257,347,356},
F1={126,137,145,234,257,356,467}.                     (8)
```

Here is a short classification proof.  The only `T`-fixed three-subset is
`356`, so it is the fixed line.  The other two lines through `5` are forced
to be `145,257`.  The line through `1,2` is then either `123` or `126`; its
four-element `T`-orbit supplies the remaining lines.  This gives precisely
(8).

After the contraction (4), both planes become the same six-vertex
multihypergraph:

```text
{12C,145,17C,24C,257,47C,C5C}.                        (9)
```

The intrinsic negative `chi_7` address line is exactly `{3,5,6}`.  Under
(4) it degenerates to the repeated-vertex line `C5C`.  The sunflower is
therefore genuinely a **contracted Fano stalk**, but it is not a simple Fano
plane on six owner obligations and its negative line cannot force three
distinct supports.

## 3. Speed chi7 is an independent bundle

THM-1156 colours owner speeds after removing their seven-adic valuation.
That colour does not transport through the address plane above.

Fix the defining difference pair `a,a+15` with `a==1 mod 15`.  The four
choices

```text
a mod 7 =1,2,3,5                                      (10)
```

realize respectively the ordered speed-colour pairs

```text
(+,+),(+,-),(-,+),(-,-).                              (11)
```

For each of the other mask representatives `5,3,2,4,7 mod 15`, impose
independently residue `1` or `3 mod 7` to choose either colour.  The Chinese
remainder theorem solves every pair of congruences, and adding distinct
multiples of `105` separates the speeds without changing masks or colours.
Consequently all `2^7=128` owner-colour words occur on the same sunflower.

The address `chi_7` and speed `chi_7` are therefore independent bundles
unless a new metric transport law couples them.

## 4. A cut-compatible full-clock packet

Take

```text
V={1,2,7,8,18,19,20},
c=1,
D=(2,7,8,18,19,20),
G=[1/14,13/14].                                        (12)
```

The defining pair `7+8=15` realizes all six masks:

```text
M_1 ={0,1,14},
M_7=M_8={0,2,13},
M_20={0,3,6,9,12},
M_19={0,4,11},
M_18={0,5,10},
M_2 ={0,7,8}.                                         (13)
```

They cover `Z/15Z` and meet pairwise only at zero after the defining pair is
identified.  On the consecutive beat block `p=2,...,13` inside `G`, the
third blockers are

```text
p=2,13:       defining pair dangerous,
p=3,6,9,12:  blocker 20,
p=4,11:      blocker 19,
p=5,10:      blocker 18,
p=7,8:       blocker 2.                               (14)
```

Every pair-safe `q=15` beat is therefore consumed.

The THM-1241 drift constraints also hold.  At pivots
`(2,7,8,18,19,20)`, the absolute invoices are

```text
(62,42,40,40,42,46),                                  (15)
```

and the fastest invoice is `46>20/14`.  Every adjacent speed gap is greater
than `20/210`, so `7->8` is itself a legal macroscopic cut emitting (13).

## 5. All carrier spokes and positive sum beats are blocked

The six centered carrier-spoke phases and complete blocker sets are

```text
d    t_d      blocker set
2    2/3      {18}
7    1/2      {2,8,18,20}
8    5/9      {18}
18   10/19    {2,19}
19   1/2      {2,8,18,20}
20   11/21    {2,19}.                                 (16)
```

For example, selecting `2->18->19->20->2` gives a literal four-cycle in
THM-1240's blocker functional graph.

More strongly, every one of the fifteen fast-fast pairs has a positive
pair-sum beat in `G` with a strict third blocker.  The table gives the beat
numerator `p`, one blocker `b`, and curvature `14D-q`:

```text
pair       p   b    14D-q       pair       p   b    14D-q
(2,7)      1   18      19       (7,20)      3   18      57
(2,8)      1   20      18       (8,18)     11    7     114
(2,18)     2   20      36       (8,19)      3   18      15
(2,19)     3    7      63       (8,20)      3   19      28
(2,20)     3    7      62       (18,19)     5    7     187
(7,8)      3   20      69       (18,20)     4   19      18
(7,18)     3    8      31       (19,20)     5    8     199
(7,19)    13    2     156
```

Thus common-clock saturation, a valid Kakeya cut, carrier-spoke cycling, and
complete positive sum-beat blocking all coexist.

## 6. Yet the packet is lonely

At `t=1/12`, the distances in the order (12) are

```text
1/12,1/6,5/12,1/3,1/2,5/12,1/3.                     (17)
```

Hence

```text
min_(v in V)||vt||=1/12>1/14.                         (18)
```

The exact-seam graph is empty as well: none of the twenty-one reduced pair
sums `(a+b)/gcd(a,b)` is divisible by `14`.  THM-1156's zero-seam
third-support mechanism is therefore not activated.

## 7. What the Fano and tournament quotients lose

On sunflower petals, first-occurrence tie gauge produces the transitive path
`A->B->C->D->E->F`; the palindromic current (6) shows why this orientation
has no metric content.

There is nontrivial tournament telemetry on carrier-spoke obligations.  For
labels `i,j`, take the cross-clearance observable

```text
c_ij=||i t_j||-||j t_i||,                             (19)
```

orient `i->j` when `c_ij>0`, and break equality by increasing speed.  Its
fingerprint is

```text
score multiset (1,1,2,3,3,5),
3 directed triangles,
SCC sizes (5,1),
9 Hamiltonian paths.                                  (20)
```

This is a genuine nontransitive tournament, yet (18) shows it still does not
encode cover truth.

The missing coordinate is off-grid metric continuation.  If
`sp==+1 mod 15`, the owning tooth extends from `p/15` by

```text
29/(210s) to the left,             1/(210s) to the right; (21)
```

for residue `-1`, the widths reverse.  The contracted Fano stalk, mask
sunflower, and blocker tournament all lose this orientation, the intervening
teeth, lift height, and owner reuse across the cyclic word.

We challenged runners, masks, sign addresses, Fano points/lines, pair beats,
carrier spokes, blocker obligations, off-grid germs, and proof obligations as
vertices.  The smallest faithful next carrier is

```text
(contracted Fano address stalk; cyclic mask word;
 oriented tooth germs and lift heights; off-grid owner continuation).      (22)
```

## 8. Verification and scope

The dependency-free referee enumerates all `30` labelled Fano planes and
finds exactly the two invariant planes (8), checks their common contraction,
all `128` independent speed-colour words, the six masks, twelve in-gap q15
beats, six carrier spokes, and all fifteen positive fast-fast sum beats.  It
also verifies the lonely and seam ledgers and the tournament fingerprint.
Normal and optimized outputs are byte-identical.

The Lean module kernel-checks the six mask classes, palindromic petal word,
pair uniqueness and invariance of both Fano planes, their common contraction,
degeneration of the negative line, full witness clock, lonely phase, all
carrier-spoke blockers, and all fifteen sum-beat blockers using ordinary
kernel reduction.  The short invariant-plane classification and CRT colour
freedom remain explicit paper proofs; there are no placeholders or
`native_decide` calls.

Frozen hashes are

```text
source         e446f579f11378514f32449469a5496f0d4a67268ae71c5b73f704001da4c19b
output         90a76fcfe439bb2ee6925c203672813ca5998ecb5be62623f0f13270dd866068
formalization  4f4756f1a6d89009a6261fcbfcae6091105c45da6d1e8bdcacde35d8bc870195
```

THM-1247 proves that the current Fano/`chi_7`, cut, mask, and sampled-blocker
constraints cannot exclude the mixed `q=15` branch.  It does not prove that
this branch embeds in an actual six-comb cover.  A viable closing theorem must
charge oriented off-grid continuation or reuse of an owner across several
address cells; projective incidence alone is exhausted.

