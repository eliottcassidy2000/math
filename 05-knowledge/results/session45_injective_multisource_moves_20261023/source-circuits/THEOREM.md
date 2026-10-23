# A native three-source dyadic cancellation circuit

**Status: PROVED algebraic wiring identity; FINITE-EXACT original order-eight
certificate and sharply scoped minimum-source census.** No global complete
atlas rewiring, all-per-port guard preservation, or all-order existence claim.

## Inheritance and connection contract

The inherited single-source-pair rectangles keep a chosen mod-eight character
but fail a different harmonic character and can be blocked by one-port
capacity. The closest positive mechanism is a common-target signed sum of
rectangles: source charges add before the harmonic guard is tested. The least
used sidecar is the actual four-word black `F/R` resource, together with its
literal Hamiltonian path and a single shared source-owner frame.

The active concepts are: (1) mixed harmonic labels; (2) integer source charge;
(3) resource-disjoint native paths; (4) shared indexed pins and source frames;
(5) the distinction between nodewise harmonic currents and individual-port
phases. The map is from four actual pinned source-path substitutions to a
signed paired-black incidence perturbation. It preserves the listed owners,
source matrices, original selected mod-eight phases, source/target margins,
and every original nodewise dyadic harmonic bank. Quotient labels alone forget
path availability, resource collisions, and owner agreement; the packet
restores those coordinates literally.

Dependencies are the original frozen order-eight owner/edge atlas, the audited
complete two-primary electrical filtration in
`metagraph_order8_electrical_20261023`, and the Session45 mixed-label compiler.
The checker directly reauthenticates the selected harmonics and every native
path. Completeness of the *entire* dyadic harmonic bank uses that independently
audited filtration, rather than inferring completeness from several examples.

## 1. An all-size conditional common-target wiring identity

Let a symmetric integer Laplacian `L` be fixed. Choose distinct targets `U,V`
and a finite source support disjoint from the targets. Suppose source `s` has
integer signed multiplicity `c_s`, with `sum_s c_s = 0`. A positive unit means
replace one old paired-black resource from `s` to `U` by one from `s` to `V`;
a negative unit means the reverse. Each black resource contributes conductance
2. Write

`a = sum_s c_s e_s`,  `b = e_U - e_V`.

Then the exact Laplacian change is

`Delta L = 2 (a b^T + b a^T)`.

Indeed the source diagonal changes cancel within every slot, the two target
diagonal changes cancel because `sum c_s=0`, and the four signed cross-entries
are precisely those of the displayed matrix. Thus source and target weighted
margins are retained. Every harmonic `x modulo m` is retained exactly when

`2 a (x_U-x_V) + 2 b (sum_s c_s x_s) = 0 modulo m`.

If some `c_s` is a unit modulo `m` (in particular, if one is `+1` or `-1`),
disjoint supports imply the equivalent two scalar conditions

`2(x_U-x_V)=0` and `2 sum_s c_s x_s=0 modulo m`.

For all harmonic `x` simultaneously, finite-ring annihilator/image duality
makes these equivalent to `2b,2a` lying in `im L modulo m`. This is an all-size
criterion, not an all-size existence proof. A physical wiring realization
additionally needs all requested native ports, one globally compatible owner
frame on every shared source, exact shared pins, and an injective selection
of the full unavailable-resource-complement inventory.

For three sources `A,B,C`, take `a=e_A-2e_B+e_C`. There are two constituent
rectangles `(A,B;U,V)` and `(C,B;U,V)`. They use TWO distinct resources at `B`;
repeating one port twice does not realize this vector.

## 2. Original order-eight packet

Targets: `U=2625`, `V=2628`.
Sources: `A=1981`, `B=2458`, `C=2659`.
All local-to-physical maps are the identity on the **same** eight owners
`0,...,7`. The two `B` slots use one identical actual source matrix and one
identical owner map, not two isomorphic source copies. Full source matrices,
all masks, both orientations, and every indexed pin are in
`three_source_packet.json`.

The following paths give the actual old-to-new orientation. Clocks start at0.

- Source1981, positive: `74135260 -> 72134650`.
  Fixed clocks `{0:7,2:1,3:3,7:0}`; indexed arc `(2,3)=(1,3)`.
- Source2659, positive: `32517640 -> 31247650`.
  Fixed clocks `{0:3,4:7,5:6,7:0}`; indexed arc `(4,5)=(7,6)`.
- Source2458, negative first: `73451620 -> 74251630`.
  Fixed clocks `{0:7,3:5,4:1,5:6,7:0}`; indexed arc `(3,4)=(5,1)`.
- Source2458, negative second: `72564130 -> 74563210`.
  Fixed clocks `{0:7,2:5,3:6,7:0}`; indexed arc `(2,3)=(5,6)`.

Each path is Hamiltonian in its stated actual source. Each side consists of
four distinct normalized black resources and16 distinct literal
`{port,Kport,Cport,KCport}` tuples; the union has eight distinct normalized
resources and32 distinct literal tuples. The two endpoint resource sets are
disjoint. The actual source path counts are493,585,575 for sources1981,2458,
2659, respectively.

Every selected per-slot original mod-eight phase is exactly2 before and after.
This is a *selected single-character* phase statement; it does not claim every
individual-port harmonic phase is retained.

## 3. Genuine cancellation, not independently safe rectangles

The original mixed module is

`Q = (Z/2)^4 direct-sum Z/128`,
`ell(s)=(Y0(s),...,Y3(s) modulo2; X(s) modulo128)`.

The concrete labels are

- `ell(A)=(0,1,1,0;23)`, with `X(A)=151`;
- `ell(B)=(1,1,0,0;127)`, with `X(B)=127`;
- `ell(C)=(0,1,1,0;103)`, with `X(C)=231`;
- `ell(U)=ell(V)=(0,0,1,1;13)`, with `X(U)=X(V)=141`.

Therefore `ell(A)-2ell(B)+ell(C)=0 in Q`, whereas both
`ell(A)-ell(B)` and `ell(C)-ell(B)` are nonzero. The two separate rectangles
have exact `X` debts respectively

`48(e_U-e_V)` and `208(e_U-e_V)`.

They individually fail modulo256, though each keeps the chosen original
mod-eight character. They also have nonzero debts against some authenticated
mod-four `Y` characters. Their joint `Y` debts vanish modulo4 and their joint
`X` debt is

`256(e_U-e_V)`,

which vanishes modulo256 but **is not zero over the integers**. The audited
mixed-module criterion therefore gives preservation of every original
nodewise dyadic bank.

The literal raw integer vector `X` is NOT harmonic modulo512 and the raw debt
is NOT zero modulo512. Beyond the terminal exponent256, the genuine harmonic
is scaled, for example `2X modulo512`, whose joint debt is512 times the target
difference. The checker verifies both sides of this boundary directly. Do
not claim an unscaled primitive lift of `X`, exact integer preservation, or
an infinite fixed-integer tower.

## 4. A sharp finite source-support boundary in a declared universe

The census fixes the following universe:

1. Every one of the3528 canonical original source matrices, including
   self-converse and nontrivial-automorphism sources.
2. All1,102,392 literal Hamiltonian paths in those matrices.
3. Actual size-four normalized black `F/R` resource orbits only; blue size-two
   resources are excluded.
4. Moves between two distinct common target types with equal mixed label,
   neither target in the selected source support.
5. Each old/new path pair keeps both endpoint owners and at least one entire
   indexed internal arc at clocks `k,k+1`, `1<=k<=5`.
6. Four unit slots: two positive and two negative. A repeated source must use
   distinct resource keys in each endpoint state. Shared physical source
   frames are the fixed canonical frames, not independently relabeled slots.

Original metagraph tilings are normalized words, equivalently `Aut(T)`-orbits,
so normalized `F/R` no-reuse is faithful to that original unit-capacity carrier.
This is not a minimum for the distinct alternative carrier that assigns
separate capacity to every hidden owner-labelled lift. For example, source43
with canonical code17301538 has paths `(7,6,5,3,2,1,0,4)` and
`(7,6,5,2,1,3,0,4)`, both normalized word344062: one normalized `F/R` tile
but two disjoint raw labelled `K/C` orbits. The selected positive packet is
unaffected: its three source automorphism groups are trivial. Keeping both
resource meanings explicit prevents hidden lifts from manufacturing extra
original tiling capacity.

The independent endpoint-fibre enumeration agrees exactly with the generator's
indexed-arc-fibre enumeration:700 eligible source/target-pair menus. Only eight
menus admit two resource-disjoint transitions. They are

`(901;649,1833)`, `(1886;2625,2628)`, `(2203;69,72)`,
`(2210;69,72)`, `(2333;69,72)`, `(2458;2625,2628)`,
`(2769;493,727)`, `(3298;69,72)`.

There is no genuine two-source doubled-rectangle repair, and exactly one
three-source arithmetic-progression repair, namely this packet (up to reversing
all signs and interchanging the two outer sources). “Exactly one” refers to
the source/target support pattern, not a claim that its detailed path/pin
presentation is unique.

Why this proves minimum source support in this universe: a nonzero balanced
source charge cannot have one source. On two sources, four slots with no
cancelling zero slots either form an independently safe single rectangle or
have source charge `2e_A-2e_B`; the latter requires capacity2 at both sources
and `2(ell(A)-ell(B))=0`, with a nonzero label difference for genuine
cancellation. The complete capacity/label test finds none. A genuine
three-source four-slot charge must be `e_A-2e_B+e_C` up to signs, and the
complete test finds exactly the stated one. Therefore three distinct source
types is the minimum for a genuine four-slot repair in this universe.

The boundary does NOT classify arbitrary target networks, different pins,
more than four slots, independently changed source-owner frames, or arbitrary
occupied remainders. More source labels can increase capacity without improving
this minimum. Nodewise guards alone never prove a globally admissible port
allocation.

## 5. Reproduction and controls

Run from this packet directory, supplying the original pinned repository:

`python3 check.py --repo /path/to/repo --out receipt.normal.json`

`python3 -O check.py --repo /path/to/repo --out receipt.optimized.json`

Both receipts are byte-identical:1,102,470 checks. The checker uses explicit
exceptions rather than validation-critical assertions. It authenticates the
original source bytes, all complete-source masks, all selected actual paths,
full black orbit keys, literal K/converse closures, exact shared source frame,
all pin values, each individually nonzero debt, the zero joint modular debt,
and the raw-X/scaled-X mod512 hostile. It independently regenerates all
eligible resource-edge menus by endpoint fibres rather than the searcher's
indexed-arc fibres.

The independent Session45 injective-completion audit also checks the frozen
packet directly at the actual raw source-matrix/path tuple level. This rules
out inventing capacity by changing labels on the repeated source.
