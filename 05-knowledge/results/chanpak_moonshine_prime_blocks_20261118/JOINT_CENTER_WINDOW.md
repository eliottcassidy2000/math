# Select the center and grow the protected choices in the same construction

**Status:** all-level physical windows, exact inverse, and complete
odd-modulus phase control for the padded family. This does not prove
universal absolute placement.

## Working backward from the two demands

If the replacement choices grow but keep the continuation, a later reader
must see one fixed arm signature while the Hamiltonian count ranges over
a palette. If the absolute center is selectable, a separate control must
move the count without magnifying the repair step. The two requirements
can be combined: stack protected arm pairs, then use a transitive deletion
as the reader. That reader has deletion count one, so it does not scale
the shifts. Its binary contacts move the center.

The starting mechanism is the [protected arm law](HALF_GRID_REPAIR.md).
An arm has a marked root; in a word all roots are identified and all
edges between successive deletions point forward. Write
(h,g,B,U,V,Ubar,Vbar) for the Hamiltonian count, deletion count, deletion
two-cover count, and signed/unsigned endpoint reads, as in that note.
Changing only h in one arm contributes its change multiplied by the
product of all the other deletion counts. Different edits add independently.

## Why not iterate 17 and 30 literally?

The old 17 and 30 half-steps came from arms with g=15,17,
and their paired two-cover count B=53760. Both paired g=255 and B
are divisible by 15. For m repeated pairs, all h variants are
divisible by 15^(m−2) when m≥2. Indeed all four endpoint reads are
divisible by 15^(m−1): the U transport sums products of m−1
factors, each a g or B, and the signed transport is a pure g power.
In the recurrence for h, the term involving h at length m−1 has a
fresh factor g; the other cross terms already contain at least
m−2 factors. Local repairs preserve the divisibility. The shared-root
transfer also keeps this divisor under arbitrary exterior arm contexts.
Thus this repetition cannot serve primes at m≥3. Coprime *local steps*
alone did not prevent a common factor in the center.

## Change the second core

Use the old g=17 pair H=65,69 with (g,B,U,V,Ubar,Vbar)
=(17,224,72,−6,135,−11). For the second pair use the following g=29
arms, independently counted by subset DP:

| H | outgoing masks | root |
|---:|---|---:|
| 91 | (38,108,56,113,99,64,5) | 4 |
| 95 | (96,1,99,7,15,26,58) | 4 |

Their common signature is (29,320,198,−20,93,−9).
One cell is the g=17 arm followed by the g=29 arm. Its all-low
profile is
```
(h_1,G,B,U,V,Ubar,Vbar) = (9689,493,71680,26406,-340,24747,-319).
```
The two upgrades add 4·29 and 4·17 respectively. This time
both B and U are coprime to G=17·29.

At depth m, repeat the cell m times. Let h_m be the count with
all low arms (iterate the seven-coordinate arm transfer), and put
```
g_m = G^m,
U_m = U (B^m−G^m)/(B−G),
V_m = V G^(m−1).
```
If j of the g=29 arms and k of the g=17 arms are upgraded,
the total is h_m+4G^(m−1)(17j+29k) with 0≤j,k≤m.
Choose, for example, the first j or k arms of each kind. All choices
have the same full signature apart from h.

The normalized palette has no gaps in
```
[476, 46m−476]          (m≥28).
```
To see this, for a target r in that interval let j_0 in [0,28]
be 12r mod29, since 17·12≡1 mod29, and put k_0=(r−17j_0)/29.
Then k_0≥0 since 476=17·28. If k_0>m, advance j by 29
until k, which drops by 17 at each step, first reaches ≤m.
It does so before j exceeds m: the largest j≤m in the correct
residue is at least m−28, at which k≤m by the upper bound on r.
The first such k is at least m−16. This also gives an exact test
for the entire palette by checking 0≤j,k≤m or rejecting.

## Invert the center, then invert the repair

The right reader is one transitive-deletion tape. Encode its incontacts
by a nonnegative integer I in binary, least significant position first.
Its signed end read is −(I mod2); its own Hamiltonian count is 1+C(I),
where C is the inversion checksum from the half-grid repair:
```
C(0)=0,
C(2J+1)=C(J),
C(2J)=C(J)+J+(J mod2)  (J>0).
```
Also 0≤C(I)≤I: an incontact in position b can contribute
at most 2^b, by summing weights over all earlier positions,
and position zero contributes nothing. The complete physical count is
```
N = h_m + U_m I + g_m C(I) − V_m (I mod2) + 4G^(m−1)(17j+29k).
```
For fixed m,I every repair in the interval above is physical at the
same continuation. If a future retained-root move needs the ordinary
end count as well, append an ordinal sink, which sets that count to
zero without changing N and transports the shared weighted starts.

For a true absolute inverse restrict to the high address annulus
```
L_m = floor(U_m/(4g_m)),        ceil(L_m/2) ≤ I ≤ L_m,   m≥2.
```
The sum of the checksum, signed correction, and repair is nonnegative
and strictly less than U_m/2. The checksum is ≤U_m/4.
For the other terms, their ratio to U_m is bounded by
```
(|V|+184m)/U · (G/B)^(m−1) ≤ (|V|+184)/U < 1/4.
```
The middle expression is largest at m=1 (the numerator grows
by less than a factor two per step, and G/B<1/2). Consequently
the absolute address is uniquely **I=floor((N−h_m)/U_m)**.
Subtract its full center, divide exactly by 4G^(m−1), and use
the palette test above. No rounding approximation is left in the answer.
The address also has to agree with an independent arithmetic read:
```
I ≡ U_m^(−1) (N−h_m)  mod G^(m−1).
```
The inverse exists since U_m≡U B^(m−1) mod G. Meeting that
congruence is not by itself a match.

The depth is also uniquely bounded by the target. Put Y_m=U_m²/g_m.
All the heights allowed here lie in [Y_m/16,Y_m]. The lower bound
uses I≥L_m/2 and U_m/g_m≥8; the upper bound uses
h_m≤g_m+U_m (the full Hamiltonian contribution is among weighted
starts), the address bound, and the residual bound above.
The envelopes are disjoint since Y_(m+1)/Y_m ≥ B²/G > 16.
Read the depth by iteration, or reject in the intervening gap.
The replay includes prime 1,870,331,984,921, decoded to
(m,I,j,k)=(2,981,1,2), and a physically built graph at depth 30
whose exact count follows by transfer. The full subset checks are
at small orders; the large-depth code checks the compiled arcs.

## Two neutral sinks solve the growing modular supply problem

Here is a knob that changes the center capacity without spending the
repair. Add two ordinary sinks inside **each** cell after its two
arms. They form a shared-root transitive arm with signature
(g,B,U,V,Ubar,Vbar)=(1,4,3,−1,0,0). The cell count and both repair
increments stay fixed. The new cell parameters are
```
(G,B',U',V') = (493,286720,107103,-493).
```
Both B',U' are still units mod G, and the exact inverse formulas
and envelope estimates hold with them. The change crosses a
load-bearing threshold:
```
B  < G²,    B' > G².
```
For the unpadded family the address interval is only
Theta((B/G)^m) long, so at that depth it supplies a vanishing
fraction of the classes mod G^(m−1). With padding the high
annulus contains **at least G^(m−1) consecutive addresses** for
every m≥2: its size is at least U_m/(8g_m)−1, which is at least
(U'/(8G))(B'/G)^(m−1)−1 > G^(m−1).
Therefore every residue mod G^(m−1) has a real center at each
depth. For m≥28 each such center also has the full growing
repair interval on its own lattice. This is joint supply, not
an inference from separate modular and count experiments.
It does not claim the full modulus 4G^(m−1).

## What still blocks an arbitrary absolute target

- The palette fills a progression, not consecutive odds. At
  fixed m,I its stride is 4G^(m−1), and the floor address of a
  particular target can still miss its required residue.
- Each layer supplies Theta(m(B/G)^m) distinct counts while
  height is Theta((B²/G)^m), replacing B by B' for padding.
  The lower bound uses the gap-free central palette at each
  address; the quotient read keeps different addresses disjoint.
  The bounded ratio of successive height scales then gives
  Theta(log X · X^s) counts up to X, with
  s=log(B/G)/log(B²/G)<1. Both lanes are provably sparse.
- A fixed finite repertoire of full-signature repairable arms
  cannot simply keep adding independent unit directions. A
  g=1 deletion is transitive (any directed triangle forces
  H≥3), and B fixes its width while U fixes its root contacts
  as a binary index. Thus no two g=1 arms with the same full
  signature have different H. For the other arms, all raw count
  shifts in a fixed frame are divisible by
  (product g_i)/(lcm of the repairable g_i), even before
  considering the local changes. Reusing a finite repertoire
  unboundedly makes this common divisor grow.

The next structural task is to allow several absolute addresses over
the *same* coarse height, or to let a continuation-preserving move
change the deletion modulus itself. Two neutral sinks show how to
increase modular supply without widening the hole in the local
palette; they do not make the absolute placements overlap.

Replay: `python3 joint_decode.py` (the single-core precursor)
and `python3 paired_window.py` (two choices, exact center, padding).
