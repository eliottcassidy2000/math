# Eight increment lanes: exact mod-8 coverage for one span-3 extension

A literal span-3 word is encoded by letter masks 0 through 7. At vertex i, bit j−1 reverses the edge to i−j, for j=1,2,3; all other edges point forward. Let h and h' be the Hamiltonian path counts before and after the final letter. Both are odd.

**Result.** Every one of the 64 ordered pairs (u,v) of odd residues modulo 16 is realized by a literal one-letter extension, with **both tournaments strongly connected** and the new order at most 7. In particular, at each fixed old odd residue u modulo 16 there are **8/8** possible residues of the *half-increment* D=(h'−h)/2 modulo 8. Equivalently, the normalized ratio h'/h modulo 16 runs over all eight units. One should not ask for arbitrary H residues mod 8: only four of the eight are possible by parity. Nor is 8/8 a factor-eight bound on numerical growth: the real sharp one-letter growth factor is 39.

The 64-row exact certificate is `STRONG_TRANSITIONS.json`. Each row specifies one short literal word, exact h and h', and the resulting residue pair. A standalone verifier uses subset Hamilton path counting and forward/reverse reachability, independent of the profile transfer used to find the words. Coverage is verified by direct integer arithmetic on these witnesses, so there is no inference from a length cutoff to all-order nonexistence in this assertion.

## Simple algebraic witnesses and a sharp boundary

If strongness of both sides is not required, a small transparent bank suffices:

| h'/h mod 16 | word | h → h' |
|---:|:---:|---:|
|1|00|1 → 1|
|3|002|1 → 3|
|5|0004|1 → 5|
|7|0021|3 → 5|
|9|00044|5 → 13|
|11|00363|1 → 11|
|13|000404|5 → 17|
|15|00046|5 → 11|

For any odd u, prepend a ribbon P with H(P)≡u/h mod 16 by ordinal sum, leaving all cross edges forward. Then the old and new counts are multiplied by the same H(P), and the quotient class stays fixed. The short H=3 and H=5 ribbons generate the full unit group mod 16 under multiplication: the eight numbers 3^i 5^j, 0≤i<4, 0≤j<2 are distinct modulo 16. This gives arbitrary old state without requiring a separate search.

The transient phenomenon is real: if the old tournament is transitive (h=1), all possible h' are
```
1,3,5,9,11,15,17,19,23,25,27,31,33,39.
```
Thus only **7/8** half-increments occur in that subfamily: D≡6 is missing, equivalently h'≡13 mod 16. The extra 5→17 row realizes precisely D=6 and ratio class 13. This exclusion is all-order, because every old label earlier than the final six is forced to precede every possible out-neighbor of the new vertex in the old transitive order. Every cycle through the new vertex starts its old portion at such an out-neighbor and then increases. The removed labels lie on no directed cycle and therefore do not affect the Hamiltonian count, by the weighted odd-cycle packing identity. The remaining six old labels preserve span 3 when restricted. Enumerating topological orders through six old labels and testing the final letter therefore suffices, and is included in the verifier.

Finally, the real ratio **39** remains sharp if one fixes *any* old odd residue mod 16 and allows arbitrary prefixes: ordinal-prepend an H(P) residue as above to the sharp 1→39 example. This gives quotient class 7 modulo 16 in every old-residue row, and hence all four corresponding pairs modulo 8. It does not prove a sharp adjusted growth constant in the other quotient classes.

## Meaning for selection

This is a one-step local coverage result, not a proof of exact integer selection, nothing about power savings for missed absolute values, and not a spectral estimate inherited from the quasi-Riemann hypothesis. The 7/8→8/8 distinction is structural: a restriction to transitive old slices loses one entire local lane, and interacting nontransitive prefixes recover it. Even the all-strong 64-row bank does not imply that one chosen prefix admits all eight extensions. The append word has only eight letters, and their consequences depend on finer state than h modulo 16.

## One strong prefix, all eight lanes with a single letter choice

The single-prefix question itself has a positive witness. The strong old word `013325` has h=25. Its eight append choices give

| append mask | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| h' |25|61|43|53|103|115|95|81|
| D mod 8 |0|2|1|6|7|5|3|4|

This gives a unique local choice for each requested half-increment, including the previously hidden D≡6. The zero letter appends a sink and preserves the count, so that extension is not strong. Every nonzero letter extends the strong prefix to a strong tournament: the new vertex has an out-neighbor by the letter, an in-neighbor by the forced long forward edge, and the old vertices form one strong component.

The three-bit map D is quadratic over the integers, and in this example is

    D = 18a + 9b + 39c − 13ab − 12ac − 13bc,

where a,b,c are the three mask bits for offsets 1,2,3. The general form follows from the odd-cycle packing formula: each new packing includes at most one cycle through the new vertex. Such a cycle has successor i and predecessor j, and its contribution to D is a nonnegative coefficient times a_i(1−a_j), with a_j=0 if the predecessor is out of range. Thus the interactions are physically forced, not independent arbitrary signs. All coefficients of this example are checked directly. The resulting permutation of eight half-increment residues is an exact small-scale steering mechanism, not a consequence of equidistribution.

Ordinally prepending an arbitrary odd-residue H(P) multiplies every D in this table by the same unit, so there is a (possibly nonstrong) single-prefix fan for **any specified odd old residue mod 16**, without using a separate extension witness per class. The strong 64-row table above proves the distinct stronger-type conclusion in which the extension also has to be strong, at the cost of allowing the prefix to depend on the target pair.


In fact all eight old residue rows have such a fan with a strong prefix: `STRONG_FANS.json` supplies the eight prefixes and the complete lists of eight exact output values, checked by subset DP. In each row masks 1 through 7 give strong extensions and the inevitable mask 0 gives a sink. Thus a single strong prefix in each old residue can steer to any **different** residue while preserving strongness. The identical residue can also be obtained strongly in the 64-transition certificate, but there the prefix must change. This qualification is necessary: an eight-way fan has exactly one choice per residue, and its zero letter always consumes its unchanged residue.

### Why the fan is not an infinite steering strategy; two-letter repair

A one-letter fan cannot support arbitrary residue requests forever while remaining a fan: a requested unchanged residue forces letter 0, and three successive unchanged requests force 000. Its three vertices form a transitive terminal block; every further appended vertex can beat only these three, so the next tournament is the ordinal sum of the older core and a tournament on four vertices with a three-vertex transitive prefix. The latter has H in {1,3,5} by the exact small count (or the transitive-prefix table restricted to three old labels). Thus only three next residue choices remain. This is a failure of the one-letter universal policy, not of general online control with variable blocks.

For the same strong H=25 fan, **two** nonzero appended letters suffice for all eight lanes, even the unchanged lane, with both intermediate and final tournaments still strong. The following columns give (ordered append masks; final H) for D mod 8 in order:
```
(2,4;73) (1,6;155) (1,1;109) (1,5;143)
(2,1;81) (4,2;275) (1,4;133) (1,2;119)
```
These too are exact checked choices, not a claim of an invariant family of two-letter fans.

### Exact finite Fourier criterion

For a three-bit extension map define P(z)=Σ_mask z^{D(mask)}. It is 8/8 iff 1+z+⋯+z⁷=(z+1)(z²+1)(z⁴+1) divides P in Z[z]. Reduction of P modulo z⁸−1 gives the eight residue multiplicities; vanishing at all nontrivial eighth roots of unity makes that remainder a multiple of the displayed degree-seven polynomial, and its value at z=1 fixes the multiple to one. This is a finite cyclotomic zero test, not a use of the quasi-Riemann hypothesis. It does not constrain where in the integers the eight exponents sit.
