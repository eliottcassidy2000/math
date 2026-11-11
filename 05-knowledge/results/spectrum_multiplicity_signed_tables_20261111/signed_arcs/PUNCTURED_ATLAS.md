# The two missing arcs are missing modulus support

Date: 2026-11-11. **PROVED:** the support dictionary, finite-rank surgery,
and terminating controller on the two-puncture target atlas. **OPEN:** that
this target atlas equals the actual Hamiltonian-path-supported atlas.
The intended interpretation is reconstructed from the user's pointers; the
mathematical statements do not depend on claiming this was the only historical
convention. No Collatz convergence assertion is used.

## Inheritance and the recovered carrier

The October 3 source
`05-knowledge/results/tiling_modular_atoms_20261003.md`, §3, supplies the
central blocks. The October 5
`05-knowledge/results/collatz_spectrum_equality_and_stitched_sheets_20261005.md`,
§§2–4, retains the actual mod9/mod11 examples and the alternating signs.
The adjacent-input shortcut-output convention is in
`05-knowledge/results/paired_prime_collatz_20261003.md`, §3. The explicit
raw-edge interpretation and four mod8 classes are inherited from
`05-knowledge/results/session255_modular_collatz_20261111/collatz_dictionary/README.md`,
§2. These different source coordinates are kept distinct.

The closest mechanism is the raw-edge bijection. The canonical hostile is
reading the mod9 center pair 2,7 as the raw positive edge 2→7. The corrected
near miss is using the entire atlas when the intended graph is supported on
attainable tournament counts. The least-used sidecar is the **table modulus
itself**. Live board: support/multiplicity, signed edge sources, first-hit
ports, cycle feedback, and gluing. Anchor = the HP spectrum; niche = finite
rank graph surgery; wildcard = a complete mixed-sign descent controller.

For positive x and sign σ∈{+1,−1}, let

    Cσ(x)=3x+σ  (x odd),       Cσ(x)=x/2  (x even),
    mσ(x)=4x+σ.

The nonzero-index multiplication table modulo m has central entries
`x,3x+σ`. For odd x read these as the raw edge. For even x=2h, decode the
second entry by `(6h+σ−σ)/6=h` and take the forward halving edge 2h→h.
Thus m≡5,1 mod8 encode the plus odd/even branches, and m≡3,7 mod8 encode the
minus odd/even branches. The two sheets are **disjoint labelled copies**
until an explicit vertex identification is made.

Let S be the set of positive Hamiltonian path totals of finite tournaments.
Define the HP-supported atlas A(S) to retain the edge `(x,σ)→(Cσ(x),σ)`
exactly when `4x+σ∈S`. Vertex destinations remain present even when their
own outgoing edge is missing. The edge has the corresponding table modulus
as an immutable label. This is an edge-supported graph, not the subgraph
induced on vertices whose integers themselves belong to S.

**Support theorem.** The map from odd moduli m≥3 to labelled atlas edges is
a bijection. Consequently the desired odd-spectrum theorem is equivalent to

    A(S)=A(all odd m≥3) minus { (5,+)→(16,+), (2,−)→(1,−) }.

Indeed m≡1 mod4 has the unique representation 4x+1 with x≥1; m≡3 mod4
has the unique representation 4x−1 with x≥1. The two omitted moduli are

    21=4·5+1,                 7=4·2−1.

The attainable count 1 has no nonzero multiplication table and is a separate
base case. This explains both requested missing arrows simultaneously.
It is an exact reformulation, rather than a proof that no further edges are
missing. In particular, deleting a vertex instead of its outgoing edge
would lose valid incoming edges and is the wrong operation.

The parent session supplies another genuine operation on these labels:
ordinal multiplication by C3 sends m→3m and therefore

    (x,σ) → (3x+σ,−σ).

At even x this operation is not a raw Collatz step. The sequence of missing
moduli 7→21→63 has labelled sources `(2,−)→(5,+)→(16,−)`. This is the
precise multiplicative reading of 2,5,16. More generally multiplying labels
`4x+σ` and `4y+τ` gives source/sign
`(4xy+τx+σy,στ)`. The parent packet proves the tournament multiplicity
version; only this elementary coordinate identity is used here.

## Exact surgery: one puncture changes feedback, the other changes flow

Take a finite forward-closed set V on one signed sheet containing the
specified edge g→t. Let A be its adjacency matrix, with rows as sources,
and delete this one edge to obtain A0. Work over formal power series in z,
so R0=(I−zA0)^−1 always exists. Adding the edge back gives

    det(I−zA)=det(I−zA0) [1−z R0(t,g)],                  (1)
    R(x,y)=R0(x,y)
       + z R0(x,g) R0(t,y)/[1−z R0(t,g)].              (2)

These are respectively the rank-one determinant identity and the geometric
series for walks that cross the restored edge. They can also be proved by
splitting a walk at its first crossing, its later returns, and its final
suffix. No assumption about unexamined infinite orbits occurs.

For the plus puncture g=5,t=16, the known suffix is

    16→8→4→2→1→4→2→1→… .

It never visits 5, hence `R0(16,5)=0`. The determinant is unchanged.
Restoration nevertheless changes every route entering 5, through the
nonzero feed-through term `z R0(x,5)R0(16,y)`. Thus the first-hit channel
detects the missing 21 even though the full cycle characteristic polynomial
does not. This is an all-finite-forward-closed-set statement, not just one
small numerical example.

For the minus puncture g=2,t=1, the deleted graph has exactly the path
`1→2` from 1 to 2, hence `R0(1,2)=z`. Restoration multiplies the determinant
by `1−z²`, and (2) has this denominator. It closes the 1↔2 cycle. The other
known minus cycles are untouched. Unknown cycles, if any, are likewise
untouched unless they used the removed edge; a functional graph cycle using
2 must be the displayed two-cycle.

The two disjoint sheets together therefore have a diagonal two-gate return
matrix `diag(0,z)`. One gate attaches a downstream orbit, the other closes a
cycle. Equal support codimension two conceals these different operations.
The first-hit numerator is a necessary additional observer for the 21 hole.

## The punctured union still has a finite terminating controller

Now explicitly identify the equal integer vertices of the two signed sheets,
keeping each edge's sign and modulus. Use the **candidate target support**
of every odd m≥3 except7,21. Stop at vertex1. There is a deterministic policy
reaching1 from every positive integer without either removed edge:

* At even n use the plus halving edge.
* At odd n≡1 mod4 other than1,5 use plus, then two plus halvings.
* At odd n≡3 mod4 use minus, then two plus halvings.
* On reaching5 use the fixed word

      5→14→7→22→11→32→16→8→4→2→1.                 (3)

The odd signs in (3) are minus at5, plus at7, minus at11. Its modulus labels
are exactly

      19,57,29,89,43,129,65,33,17,9.

These are actual attainable small HP totals; the accompanying receipt gives
literal tournament witnesses, inherited where convenient and checked by a
small subset DP. In particular the forbidden minus edge2→1 is replaced by
the plus edge with label9, while the forbidden plus edge5→16 is replaced by
the whole detour. Edge labels are therefore indispensable even when endpoints
coincide.

**Proof of termination.** Away from1,5, the selected odd numerator is
divisible by4. Its three-edge block ends at `(3n+σ)/4<n`; every even block
halves. All such blocks end at most `7n/9`: the plus case has n≥9; the minus
case is below3n/4. Thus block endpoints decrease until1 or5 is reached.
The exceptional detour ends at1 in ten edges. A bound of
`3 ceil(log_(9/7)(n))+10` edges is sufficient (with the n=1 case immediate).
No arithmetic logarithm is needed for the algorithm; it only tests residue
classes and the exceptional state5.

This theorem shows the two specified punctures are compatible with a small
complete mixed-sign route grammar. It does **not** provide tournament
witnesses for every modulus along arbitrary routes: that is precisely the
remaining support theorem. Nor is it a policy restricted to one Collatz sign.
Its useful finite port is the concrete detour (3), not a claim that numerical
descent manufactures an unspecified HP source.

## Reproduction and transfer contract

Run `python3 -B check.py` in this directory; `--write` regenerates the receipt.
The code checks exact raw decoding, determinant polynomials by an independent
cycle decomposition, first-hit walk counts, the finite controller, and the
small literal source witnesses. No broad Collatz orbit census is used.

Source→target: attainable odd total→one labelled signed edge. Preserved:
support, with exact exceptional labels. Lost: the tournament, order,
multiplicity, and extension port response. Needed sidecar: a source witness
and its actual constructor. Cheapest decisive test: the same integer edge
2→1 has labels7 and9 on opposite sheets, only the latter allowed.

Source→target: punctured functional graph→its cycle determinant. Preserved:
cycle feedback. Lost: the whole transient attachment at5. Needed sidecar:
the first-hit numerator in (2). Decisive test: restoring5→16 changes the
1-hitting series while leaving the determinant exactly unchanged.
