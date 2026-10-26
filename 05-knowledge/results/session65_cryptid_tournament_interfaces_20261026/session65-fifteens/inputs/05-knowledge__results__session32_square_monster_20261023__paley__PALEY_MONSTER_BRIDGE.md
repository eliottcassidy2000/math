# Paley7 and the Monster 7A normalizer: an exact shared quotient action

**Status.** CITED Monster local-group structure; PROVED elementary transport to a seven-point affine action; FINITE-EXACT Paley7 counts. This is a concrete group-action bridge, not an explanation of forbidden Hamiltonian counts 7 and 21. Research date: 2026-10-23. No Git operations were performed by this worker.

## 1. Inheritance and scope

Closest proved mechanism: the square-affine ordinary group and its converse extension, already separated in `01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`. Canonical hostile: Paley7 has ordinary regular vertex subsets but no regular subsets under the signed group. Corrected near miss: `01-canon/MISTAKES.md`, MISTAKE-378, and the symmetry repair in `01-canon/theorems/THM-309-five-cycle-design-paley.md`, distinguish square slopes from nonsquare orientation reversals. Only that repaired symmetry statement is inherited here, not unrelated cycle claims in THM-309. Relevant sidecar: `05-knowledge/results/deficit_six_centers_20261021b/README.md`, section 5, and its `check.py` lines 534–543 explicitly identify the order-seven center as Paley7 with `(H,Aut,h,b)=(189,21,9,3)`.

Concept board: quadratic-residue orientation; normalizer quotients; coset action versus conjugation; signed path orbits; deficit-six centers; parafermion proof boundary. The new connection concerns the first four, validates the fifth at this one center, and does not import a conjectural assertion from the sixth.

## 2. The exact seven-point action

On `V=F7`, put `Q={1,2,4}` and orient `x -> y` exactly when `y-x` belongs to Q. The maps

`G={x -> ax+b: a in Q, b in F7}`

preserve every arc, and form `C7:C3`, of order 21. Generators are `t:x -> x+1` and `s:x -> 2x`, with `s t s^(-1)=t^2`.

This is the full ordinary automorphism group, not just a subgroup. An automorphism fixing 0 preserves its out-neighborhood `{1,2,4}`, which is a directed triangle, so its restriction has at most three choices. This restriction determines the remaining vertices: the out-neighborhoods in Q of 3, 5, 6 are respectively `{4}`, `{2}`, `{1}`. Thus the point stabilizer has order at most 3; translations are transitive, giving at most 21 automorphisms. The displayed 21 maps attain the bound.

The map `x -> -x` reverses all arcs. Composing any anti-automorphism with it gives an ordinary automorphism. Consequently

`Aut(P7)=7:3`, and `Aut±(P7)=AGL(1,7)=7:6`.

Square multipliers preserve orientation and nonsquare multipliers reverse it. The ordinary group is sharply transitive on the 21 directed arcs, equivalently on the 21 unordered pairs. The full affine group is sharply two-transitive, but is not an ordinary tournament automorphism group. The orientation choice is an actual sidecar: the two orbital tournaments of the order 21 action are converse to one another.

## 3. The Monster source and quotient

For a Monster element g in class 7A, set `C=<g>` and `N=N_M(C)`. The cited local structure is

`C_M(C)=C7 × He`, `N=(7:3 × He).2`.

A primary reference is Jianbei An and R. A. Wilson, *The Alperin weight conjecture and Uno's conjecture for the Monster M, p odd*, LMS Journal of Computation and Mathematics 13 (2010), 320–356, Proposition 4.1, proof Case 1, p. 324; DOI [10.1112/S1461157009000059](https://doi.org/10.1112/S1461157009000059), [author-hosted PDF](https://webspace.maths.qmul.ac.uk/r.a.wilson/pubs_files/S1461157009000059a.pdf). The date in the DOI is not its publication year.

The action of the outer factor is not inferred merely from dot notation. Thomas Breuer, Kay Magaard and Robert A. Wilson, *Some steps in the verification of the ordinary character table of the Monster group*, dated 12 December 2024, explicitly describe N as the subdirect product of `7:6` and `He.2`, over their common quotient C2. Section 4.13, p. 20 states this; section 5, p. 22 gives the corresponding GAP construction. [Author-hosted computational supplement](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf). This is imported group structure, not a new construction of the Monster.

Here is the precise deduction. Since He is perfect,

`[C_M(C),C_M(C)]=He`.

The centralizer is normal in N; its derived subgroup is characteristic, hence He is normal in N. The cited subdirect description then supplies

`N/He ≅ 7:6 ≅ AGL(1,7)`.

Equivalently, the quotient has a normal C7 and faithful C6 conjugation action on it. Choosing a complement C6 gives the degree-seven action on its left cosets; in affine coordinates `(b,a)` acts by `x -> b+ax`. Its index-two square-multiplier subgroup is exactly `Aut(P7)`. Thus the normalizer has a seven-point signed action with kernelHe, and the index-two subgroup `(7:3 × He)` acts through the ordinary Paley group. The order 21 factor also embeds as a subgroup in N.

**Do not confuse the two actions.** Conjugation of N on C itself fixes the identity and has image C6, not7:6; C acts trivially by conjugation on itself. Translations in the faithful quotient action arise from the coset/affine action, not from conjugation. A coordinate, complement, and orientation choice are needed to identify this quotient action with the displayed labeled tournament. No action of the entire Monster on seven points is asserted.

The common carrier therefore preserves affine incidence, quadratic-residue orientation up to its sign, and the corresponding ordinary/signed permutation actions. It discards the Held kernel and almost all Monster data. In particular, neither a Monster element's full conjugacy information nor a VOA representation can be recovered from this quotient.

## 4. Exact Hamilton words and presentations

H means the number of ordered vertex words following tournament arcs; no rotation quotient is taken. An ordinary presentation is an `Aut(P7)` orbit of these words. The ordinary group acts freely because a word lists every vertex.

Every Hamilton word has a unique square-affine normalization whose first two entries are 0, 1: subtract the first vertex and divide by the first step, which is a nonzero square. Exhaustive enumeration, independently checked by subset dynamic programming, gives

`H=189`, `|Aut|=21`, `h=H/|Aut|=9`.

A word p is blue precisely when its positional reversal `r_p(p_i)=p_(6-i)` is an anti-automorphism. Exactly 63 words, or 3 ordinary presentations, are blue. Hence `b=3`, `h-b=6`, and the six nonblue presentations form three converse pairs. This recovers the previously verified Paley7 deficit-six center, not a classification of all deficit-six centers.

An anti-automorphism acts on Hamilton words only after also reversing their positions. Using `sigma(x)=-x`, the induced involution is `p -> normalize(sigma(reverse(p)))`. Its three fixed points are exactly the blue presentations; its other six points form three pairs. The signed action on 189 words has orbit sizes `21,21,21,42,42,42`. The full signed group has no regular vertex subset, whereas the ordinary group has 84; all 128 subsets are checked.

This is a genuine transport of the global converse operation from the quotient's nonsquare coset. It is **not** the fixed-Hamilton-path operation that complements only nonpath arcs in the repository's metagraph wiring construction. A vertex permutation/global converse and that presentation-dependent edge operation have different types. The quotient does not predict metagraph wiring.

## 5. Parafermion route: citation and stopping boundary

Jin-Beom Bae, Jeffrey A. Harvey, Kimyeong Lee, Sungjay Lee and Brandon C. Rayhaun, *Conformal Field Theories with Sporadic Group Symmetry*, Communications in Mathematical Physics 388 (2021), 1–105, DOI [10.1007/s00220-021-04207-7](https://doi.org/10.1007/s00220-021-04207-7), [arXiv:2002.02970v2](https://arxiv.org/pdf/2002.02970), section 3.3, equations 3.173–3.175, distinguishes established inputs from expectations. Equation 3.173 gives the local group structure above. Equation 3.174 supplies the subalgebra chain

`P(7) -> V^+_(sqrt(2) A6) -> V^+_(sqrt(2) N(A6^4)) -> V^+_Leech -> V^natural`.

The lattice inclusions and the cited parafermion subalgebra theorem establish existence of this embedded parafermion VOA. Here `P(7)` is the parafermion algebra, not the Paley tournament P7. The text only expects a Held action on its commutant and supplies evidence for the fusion automorphism lifting into class 7A. Its proposed Held commutant/M-com construction has additional conjectural qualifications. None is used to prove our quotient bridge. In particular, this paper supplies no map from Paley Hamilton words or presentations into VOA states.

## 6. What this does and does not establish

The bridge is `N_M(<7A>)/He -> AGL(1,7)=Aut±(P7)`, with index-two orientation-preserving subgroup `7:3`. It has an explicit degree-seven action and exact finite tests. It is stronger than a numerical coincidence and weaker than a Monster construction or a tournament-to-VOA functor.

The numbers 7 and 21 occur here as the vertex count and ordinary group/arc count. The actual Paley Hamilton count is 189, not 7 or 21. Nothing here proves why 7 and 21 are forbidden values of H over all tournaments, characterizes those gaps, identifies the Baby Monster by a Hamilton statistic, or promotes scalar coincidences into a moonshine theorem. That explanation remains **OPEN / UNSUPPLIED by this bridge**.

## Reproduction

Run `python3 -B check.py > exact.json` and `python3 -B -O check.py > optimized.json`, then compare the files. The script uses explicit runtime checks, not removable assertions. It enumerates all 5040 permutations and all 128 subsets, compares permutation-discovered groups to affine formulas, verifies the independent subset DP, retains all nine normalized presentations and converse pairings, and includes its own SHA-256 in the result. This is finite-exact verification of P7, not a proof by computation of any unbounded classification.
