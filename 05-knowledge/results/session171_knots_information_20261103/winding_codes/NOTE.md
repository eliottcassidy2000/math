# Winding codes: sharp hidden-return lengths on a constant-spectrum tournament carrier

**Status: PROVED elementary all-width statements, conditional only on the inherited six-vertex seed certificate; FINITE-EXACT independent seed and shortest-path checks.** The coding parameters are standard linear-code facts proved below. No priority claim is made. The odd Hamiltonian-count spectrum conjecture remains OPEN.

## Inheritance and connection contract

The closest mechanism is [S170, monodromy, §§3–5](https://github.com/eliottcassidy2000/math/blob/2d137c00eaf36a6fbe9f06044ac56d2c3bf7c7d6/05-knowledge/results/session170_holonomy_homotopy_20261103/monodromy/NOTE.md): the physical product of six-state tournament hexagons covers the product of three-state logical circles, and winding modulo two is precisely the physical return. S168 supplies the literal six masks and full-profile certificate. THM-4094, `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`, supplies the strong-prime-lane spectrum boundary.

The canonical hostile is one coordinate triangle: its class returns after three flips, but its pinned tournament changes by three arcs. The corrected near miss is treating class return or a quotient return as literal source return. The retained sidecar here is a parity-check syndrome of the winding vector. Unlike S170's nonabelian all-or-nothing sign quotient, the present commuting deck group admits an entire linear-code lattice of compatible observations.

**Source → target:** a binary linear code C in F₂ᴺ → an intermediate cover of the inherited torus cover. The map quotients deck translations by C. It preserves logical class states, literal edge projection, and constant path profile; it destroys exactly the deck returns in C. The missing information is an element of C, recoverable only with additional observations or a promised bound on the possible error. The decisive test is the shortest nonzero codeword, realized by three actual arc flips per nonzero coordinate.

Concept board: Hamiltonian spectrum / pinned versus class return / winding lattice / coding distance / observation capacity. This is an exact bridge, not an assertion that coding or topology creates the missing prime Hamilton counts.

## 1. Physical setup

In lexicographic pair-bit order on six vertices, use

    B = (153, 8345, 8337, 8369, 177, 185).

Consecutive masks, including the last to the first, differ in one arc. Antipodal masks B_h and B_(h+3) differ in exactly three arcs. Each seed is strong and has H=29 and labelled-slot polynomial

    f(r)=r⁶+12r⁴+16r².

On N fixed ordered blocks let T_h be their ordinal sum, h∈(Z/6)ᴺ. The inherited physical graph is C₆^□N, logical class graph C₃^□N, and each state has H=29ᴺ and full profile f(r)ᴺ. Taking the full cubical products of the polygon circles, including all mixed squares and higher cubes, yields the cover

    E=Rᴺ/(6Zᴺ) → B=Rᴺ/(3Zᴺ).

The deck group is F₂ᴺ, acting by h↦h+3c. A closed logical path of winding w∈Zᴺ returns to deck vector e=w mod2.

The codeword zero always means the original pinned source. A nonzero codeword is not the same labelled tournament, even though it is isomorphic to it. Intermediate sheets below are equivalence classes of pinned tournaments, not new tournament isomorphism classes.

## 2. Linear codes are exactly the connected intermediate winding covers

Let C≤F₂ᴺ have dimension k and let A be a full-row-rank (N−k)×N binary matrix with kernel C. Define

    L_C={w∈Zᴺ : w mod2∈C}.

**Theorem 1.** Quotienting E by the deck subgroup C gives the connected cover

    E_C=Rᴺ/(3L_C) → B,

of degree 2^(N−k). Its monodromy on a loop w is exactly the syndrome A(w mod2). A logical loop closes in E_C iff its deck return lies in C. It closes at the original physical source iff its deck return is zero. Every connected intermediate cover of E→B arises uniquely this way.

**Proof.** The lattice 3L_C contains 6Zᴺ, and its quotient by 6Zᴺ is C. Its index in 3Zᴺ is [F₂ᴺ:C]. Lifting a loop of displacement 3w returns modulo 3L_C exactly when w∈L_C. The map F₂ᴺ/C→im A given by c+C↦Ac is a bijection, proving the syndrome description. An intermediate connected cover determines a subgroup between 2Zᴺ and Zᴺ; reducing it modulo 2 gives C and recovers L_C. This last claim can also be obtained directly by identifying deck translations whose endpoints become identified. ∎

Thus N−k syndrome bits describe the quotient sheet, while k further bits are needed in the worst case to recover an arbitrary physical source sheet from a known quotient sheet. There are exactly 2ᵏ possibilities in that fiber. Nonlinear encodings cannot avoid this cardinality lower bound.

## 3. Code distance is an exact physical hidden-return threshold

Write d(C)=min{weight(c):0≠c∈C}, with d({0})=∞. Count a path by its number of literal single-arc flips.

**Theorem 2 (sharp hidden-return length).** Among paths starting at a fixed pinned tournament, projecting to a closed logical path and closing in E_C but not at their pinned source, the shortest length is exactly

    3 d(C).

The same number is the minimum Hamming distance between arc masks of distinct pinned tournaments identified by the code quotient within one logical fiber.

**Proof.** A logical loop changes coordinate i by an integer multiple of three signed unit steps. If its parity return c_i is one, that sum is an odd multiple of three, requiring at least three steps in coordinate i. Consequently any hidden nonzero return costs at least 3 weight(c)≥3d(C). For a minimum codeword c, make one positive three-edge triangle in every coordinate with c_i=1. The resulting path has exactly 3d(C) steps, returns logically and in the quotient, and ends at h+3c≠h. Each changed block is antipodal and differs in exactly three independent arcs. Therefore the final arc masks differ in exactly 3 weight(c), proving also the second assertion. ∎

This lower bound survives paths through other labelled tournaments when measured only by endpoint arc-flip distance: one flip changes one arc, so no route can beat the endpoint Hamming distance. The cover interpretation and constant-profile intermediate vertices are asserted only for the displayed carrier.

**Corollary.** For paths of length strictly less than 3d(C), quotient closure certifies actual source closure. The strict inequality is sharp. A syndrome observes every nonzero deck error of support at most t iff d(C)>t. It uniquely determines any promised deck error of support at most t iff d(C)>2t.

For the last assertion, two such errors have the same syndrome iff their difference is a codeword of weight at most 2t. Conversely any codeword of weight at most 2t can be partitioned into two supports of size at most t. Thus this is an iff, not only a sufficient decoding bound.

An r-bit syndrome correcting t deck errors must satisfy

    2ʳ ≥ Σ_(j=0)^t binom(N,j),

because all errors of support at most t must have distinct syndromes. The guarantee assumes a bound on the actual parity-return support. It does not correct wrong coordinate labels, an incorrect circuit description, or arbitrary errors in the syndrome readout itself.

## 4. Concrete order-42 controls

For N=7 all displayed tournaments have order 42 and H=29⁷, independent of the chosen code or observation. Five choices illustrate the exact tradeoff:

- C=F₂⁷: no observation; first hidden return is 3 flips.
- C=even-weight code [7,6,2]: one syndrome bit; first hidden return is 6 flips.
- C=Hamming code [7,4,3]: three syndrome bits; first hidden return is 9 flips, and all single-coordinate deck errors are uniquely correctable.
- C=repetition code {0,1111111}: six syndrome bits; first hidden return is 21 flips, and every deck error on at most three coordinates is uniquely correctable.
- C={0}: seven bits; there is no hidden return.

For the Hamming case take the seven columns of A to be the seven nonzero vectors of F₂³. There is no zero column or repeated column, so no codeword has weight one or two. Columns 001,010,011 sum to zero, giving distance exactly three. The seven nonzero single-error syndromes and the zero syndrome exhaust F₂³, meeting the counting bound for t=1.

The repetition case is implemented by checks c_i+c_7 for i=1,...,6. Its only nonzero codeword has weight seven. A quotient can therefore retain only 64 of 128 physical sheets yet certify every quotient-closing loop shorter than 21 flips as a true source return.

## 5. Hamiltonian spectrum firewall and amplifier

**Proposition 3.** If a fixed tournament R has H(R)=q, ordinal adjoining the N-block carrier produces pinned states R⊕T_h with

    H(R⊕T_h)=q·29ᴺ,
    F_(R⊕T_h)=F_R fᴺ.

Every code, winding syndrome, and hidden-return threshold above persists inside this same scalar H value. If R is nonempty, the tournament order is |R|+6N.

**Proof.** Every directed Hamilton path in an ordinal sum traverses its factors in order and restricts to a Hamilton path in each. Concatenation is a bijection, giving multiplicativity. The inherited labelled-slot identity gives the full-profile statement. The flips occur inside the designated blocks, so all metric and winding arguments are unchanged. ∎

This provides arbitrarily large, tunably protected return information inside single attainable odd Hamilton counts. It does not provide any new value outside the multiplicative family q·29ᴺ. In particular, for N≥1 it cannot construct a prime value other than 29. The strong-prime realization lanes of THM-4094 remain untouched. Taking N=0 removes the coding carrier and merely recovers the pre-existing witness R.

Any observation computed only from H, or even from the full path-cover profile, is constant on this carrier. It cannot certify source return. The lower bound is not vague information loss: the code identifies exactly which 2ᵏ pinned return states are invisible, and its distance gives the exact first physical failure. Conversely, arithmetic spectrum realization and return-to-source certification are independent obligations; solving the latter on this carrier does not solve the former.

## 6. Exact checks and scope

Run `python3 check.py` and `python3 -O check.py`. Outputs agree byte for byte. The script uses explicit exceptions rather than assertions.

The checker independently counts Hamilton paths by all 720 permutations in every seed, checks actual arc differences, enumerates all 16 linear subspaces of F₂³, and for each uses BFS in the actual 216-vertex C₆³ graph to find the first hidden quotient-closing return. BFS is independent of the distance formula. It agrees with 3d(C) in every case, including no hidden return for the zero code. It separately enumerates the seven-bit control codes and directly measures endpoint tournament arc differences.

The theorem does not assert that source membership in an arbitrary tournament construction is linear, that a quotient sheet is an actual new tournament class, or that topology forces a requested odd Hamilton count. The next useful targeted experiment is to identify a naturally occurring source-language constraint whose return syndromes already form a code, rather than declare an arbitrary code and treat that choice as an intrinsic property of the unmarked tournament graph.
