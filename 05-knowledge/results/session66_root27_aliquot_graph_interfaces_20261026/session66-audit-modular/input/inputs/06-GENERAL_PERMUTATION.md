# General signed pair-orbit defect polynomial

**2026-10-26. Status:** elementary proof, with complete FINITE-EXACT checks described below; independently hostile audit pending. Proposed by the Session65 root while reviewing the involution specialization. No novelty claim: the mechanism is the classical signed-cycle/parity projection already used in repository signed-monomial work. It adds a quantitative form to the odd-order automorphism obstruction.

Let `σ` be any permutation of a finite labelled vertex set `V`. For each tournament `T` define `D_σ(T)` as in `REPORT.md`, counting unordered pairs whose direction is not preserved by `σ`. Let `O` run over the orbits of `σ` on unordered pairs. Put `L_O=|O|`, and let `epsilon_O=+1` if `σ^{L_O}` fixes the two endpoints individually, and `epsilon_O=−1` if it swaps them.

## Theorem

`F_σ(z):=Σ_T z^{D_σ(T)} = Π_O [(1+z)^{L_O}+epsilon_O(1−z)^{L_O}]`.

**Proof.** Choose one arbitrary ordered representative for each unordered pair in a fixed pair-orbit. Each edge orientation becomes a bit. Transport by σ gives an affine cyclic system of comparisons; its total twist is zero or one according as `epsilon_O` is positive or negative. If `d_i` records the failure bit at each step, telescoping around the orbit yields

`Σ_i d_i ≡ (1−epsilon_O)/2 (mod2)`.

Conversely every failure-bit pattern satisfying that parity constraint reconstructs exactly two edge assignments, obtained by selecting the first edge's orientation and propagating. Thus the orbit's generating factor is twice the even or odd part of `(1+z)^L`, namely the displayed factor. Different pair-orbits have disjoint edge variables, giving the product. QED.

## Which pair-orbits have negative twist?

Exactly one occurs for each **even vertex-cycle** of σ. If the endpoints belong to different vertex-cycles they cannot be swapped. For endpoints in the same cycle of length `d`, a swap after `L` steps implies `L=d/2` and the endpoints are antipodal, forcing `d` even. Conversely the antipodal pairs of an even `d`-cycle form one pair-orbit of length `d/2`, and its transport swaps endpoints.

Let `e(σ)` be the number of even vertex-cycles, and `q(σ)` the number of all unordered-pair orbits. Then

`min_T D_σ(T)=e(σ)`,

`D_σ(T) ≡ e(σ) (mod2)`,

`#{T:D_σ(T)=e(σ)}=2^{q(σ)} Π_(even vertex-cycles C) (|C|/2)`.

An invariant tournament exists iff every vertex-cycle is odd, equivalently iff σ has odd order. The sign of σ alone retains only the defect parity: `(12)(34)` is an even permutation but has minimum defect two. This is an explicit hostile against replacing the full cycle sidecar by permutation parity.

The involution theorem follows with `C(f,2)` positive length-one pair-orbits, `b` negative length-one pair-orbits, and `m=b(f+b−1)` positive length-two pair-orbits.

## Controls and connection contract

`check_general_permutation.py` builds the actual signed pair action from each vertex permutation, enumerates every tournament, and compares the complete histograms with the formula. It also checks the negative-orbit classification, minimum, minimizer count, and parity.

Finite universe: one representative of **every permutation conjugacy type for n=0,…,6**, with every tournament on its vertex set: 30 cases and 367,966 tournament/permutation combinations. Conjugating vertex labels transports tournaments bijectively and preserves the defect, so cycle-type coverage is complete for this invariant. Normal and `-O` executions are byte-identical.

- A 3-cycle has `F(z)=2+6z²` and minimum zero; the two minimizing tournaments are its two cyclic orientations.
- A 4-cycle has `F(z)=8z+48z³+8z⁵` and minimum one.
- Two transpositions on four vertices have `F(z)=16z²+32z⁴+16z⁶`, despite the permutation having positive sign.

Source: a literal vertex permutation plus complete pair incidence. Target: signed pair-orbit lengths, then defect polynomial. Preserved: complete defect histogram; ordinary-invariance existence; minimum. Lost: labels, arithmetic realization and any original tournament's specific edges once one sums over T. Sidecar: the actual action and its signed pair transport; an arithmetic application must separately justify that this action is the relevant one. Hostile: retaining only `sgn(σ)` fails at two transpositions. The polynomial is an enumerator over orientations, not the Session63 subset-orbit polynomial, and its defect order is not a legal-toggle rank or an arithmetic orbit-realization theorem.

Reproduce:

```
python3 check_general_permutation.py > general-receipt.json
python3 -O check_general_permutation.py > general-receipt-optimized.json
cmp general-receipt.json general-receipt-optimized.json
```
