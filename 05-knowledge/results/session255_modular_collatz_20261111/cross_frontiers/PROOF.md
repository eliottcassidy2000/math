# Finite repair compilers and the exact location of modular memory

**PROVED IN NOTES:** the general support/residue compiler and affine fibre theorem below. **FINITE-EXACT:** eleven-state minimization of the concrete four-coordinate equality language and the saved controls. Inherited tournament formulas are used at their stated scope, not audited anew. No literature-priority claim. Collatz, Chan–Pak, JC(2), LRC(14), and the targeted HN improvement remain open.

## Inheritance, board, and actual maps

Starting source: main `6bb050cb3a0d`. The closest mechanism is the nonnegative common frame in `05-knowledge/results/h13_paired_descent_20261111/dynamics/PROOF.md`. Its actual 401 repair is in that packet's `dynamics/INCOMING_REPAIR.md`. Session254's `SYNTHESIS.md` and `chan_pak/PROOF.md` supply the distinction between a legal repair and a formal modular replacement. The new all-size birth theorem `infinite_birth_families_20261111/BIRTH_THEOREMS.md`, B4–B7, says why scalar transition edges do not compose without their common source. The bounded-defect boxes and growing reset family in `uniform_return_metagraph_20261111/defect/BOUNDED_TEMPLATES.md` and `reset/COMPOSITION.md` give the other boundary: finite classification at each fixed budget is not a uniform unbounded reset theorem.

Canonical hostile: two equal H values need not have equal future responses. Corrected near miss: smaller response does not preserve exact future target values. Least-used sidecar here: the **support of the difference in one common positive coordinate frame**, rather than its norm or the individual count.

Concept board: (1) modular multiplication kernels, (2) signed Collatz affine words, (3) positive rooted response matrices, (4) actual count-neutral repair languages, (5) fixed-defect versus unbounded-state finite models. Anchor is a native target-preserving repair compiler; niche is exact 3-primary memory; wildcard is modular access without integer access.

## 1. General finite compiler theorem

Let a finite alphabet A act by nonnegative integer d-by-d matrices M_a. Let x,y be integer states with delta=x−y coordinatewise nonnegative. Fix an output coordinate i, an odd or even modulus m>=1, and any predicate R on residues (Z/m)^d. The following language is regular:

    words w with (M_w x)_i = (M_w y)_i
    and R(M_w y mod m).

It is accepted by at most 2^d m^d states. Consequently, if such a word exists, one exists of length at most 2^d m^d−1. A shortest word can be found by finite breadth-first search. The theorem also permits a predicate on any specified finite-state sidecar, by taking the product and multiplying the bound by its number of states.

**Proof.** For S=support(delta), define

    B_a(S)={j: there is k in S with (M_a)_(j,k)>0}.

No cancellation is possible, so support(M_a delta)=B_a(S) exactly. The second coordinate updates by v↦M_a v modulo m. Start at (support(delta), y mod m); accept exactly when i is absent from S and R(v) holds. Induction proves the asserted equivalence for every word. Any repeated state can be removed from an accepting path, proving the length bound. No bound on the growing integer vector is required. QED.

For several moduli m_j one should use their lcm, or the compatible residue fibre product, not the unrestricted Cartesian product when factors overlap. Coprime prime powers give exactly the m^d residue states by CRT. Every factor must read the **same operation letter** at each step. Independently chosen local words do not glue.

**What this manages:** an infinite family of actual compositions of the supplied operations, subject to exact equality of the two outputs and a modular condition. **What it does not manage:** arbitrary equal-target repairs, signed cancellations, or the exact integer value of the final common output. When x,y are realized by actual starting tournaments and all letters are the inherited actual O/Z/J births, the returned word is a literal construction on both graphs. The full starting graphs remain the realization sidecar; a residue or response vector alone is not a graph.

The assumption delta>=0 is essential. In dimension two with output row (1,1), differences (1,−1) and (1,−2) have the same signed support but output respectively zero and −1. A sign pattern does not replace magnitudes outside the positive cone.

## 2. An exact eleven-state equality machine for O/Z/J

In the inherited invertible coordinates

    q(H,E,U,d)=(h,u,v,t)=(H,U,U−d,H+U−E),

O, Z, J act as

    O(h,u,v,t)=(h,2u,u+v,h+2u),
    Z(h,u,v,t)=(t,u+v,2v,t+2v),
    J(h,u,v,t)=(t,h+u,u,h+2u).

Thus all count-neutral future words of every comparable pair of sources depend on only four support bits. This **does not require delta to be bounded**. Number bits h,u,v,t as 1,2,4,8. Starting from all sixteen formal supports, ordinary deterministic partition refinement, with acceptance h=0, yields exactly these eleven language-equivalence classes:

    {0}, {2}, {4}, {6}, {8}, {12}, {10,14},
    {1}, {3}, {5,7}, {9,11,13,15}.

This eleven-state claim is FINITE-EXACT, certified by the complete transition table and refinement in probe.py/receipt.json. It is not a claim that every formal support is attained by a pair of tournament states. The generic proof needs only sixteen states. Using the displayed quotient improves the response-congruence witness bound to **11m^4−1**.

Two examples recover and extend concrete repairs:

- The inherited H13→H5 ray has delta=(0,8,0,0), support2. Its neutral language is Z*O* union Z*JO*, as already proved in the H13 note. This is a positive control, not a new attribution.
- The actual two-edge 401 repair has delta=(0,12,12,12), support14. Its exact neutral language is O*. Indeed O preserves h=0 and makes t positive; Z or J sends the positive t into h, and thereafter h and t remain positive. Thus the repair can be placed before arbitrarily many universal sinks at no count cost, but not before any Z/J birth at exact equality.

The **new compiler** additionally asks for a downstream modular response while staying in the exact neutral language. On the cheaper actual 401 source q=(401,573,532,933), the neutral O^k words have U=573·2^k. Hence requesting U≡1 modulo m reduces to an explicit discrete cyclic-orbit condition, not a vague modular repair promise. For example modulo5 the shortest word is O; modulo3,7,9 there is no solution. Modulo7 the source is a unit but lies in the wrong doubling coset: 573·2^k visits only6,5,3. Modulo11 the empty word already succeeds. This separates a nonunit obstruction from a unit-orbit obstruction. All statements concern the native response state of actual appended-sink tournaments, not only abstract affine maps.

This is useful for the bounded-template frontier: every finite library of actual comparable reset pairs can be compiled separately and then queried under arbitrary future O/Z/J words. There is no need to expand graph order to a guessed cutoff. It cannot make the fixed-defect library cover unbounded defect, and it cannot identify two scalar birth edges whose middle graphs differ.

## 3. Exact fibre theorem for signed Collatz words

Fix an odd modulus m and a prescribed word of r signed odd steps

    x_(j+1)=(3x_j+epsilon_j)/2^(k_j),
    epsilon_j in {−1,+1}, k_j>=1.

The signs and exponents are part of the word. This defines an affine modular map even when some ordinary integer inputs would not follow that word. Put K=sum k_j. Its multiplier is a=3^r 2^(−K) mod m, and its additive term b is obtained by ordinary affine composition.

**Theorem.** Put d=gcd(3^r,m). The modular word has exactly m/d outputs, each with exactly d inputs. Its outputs are precisely b+d(Z/m). Its input ambiguity has a lossless sidecar of d symbols, and no smaller sidecar can recover every input from its output.

**Proof.** Multiplication by a has kernel of size d since 2 is a unit. Its image is the subgroup d(Z/m). Translating by b gives the claim. All solutions in a fibre have a fixed residue c modulo m/d. Write each input uniquely as c+j(m/d), 0<=c<m/d, 0<=j<d. The output together with j recovers c, since multiplication by a/d is invertible modulo m/d. A fibre of d elements requires at least d distinct sidecar values by injectivity. QED.

Writing m=3^e q with gcd(3,q)=1, this becomes d=3^min(r,e). Thus **only the 3-primary component loses starting-state memory under a fixed signed word**. Every modulus coprime to3, including composites such as25 or35, gives a permutation. The prime/composite multiplication-table split is real for arbitrary nonzero rows, but it does not justify treating all composite Collatz channels as lossy: loss depends on the multiplier and its local valuation. Prime3 itself is already a completely collapsing multiplier3 channel.

For a general fixed affine word with numerator multipliers a_j and unit denominators, replace 3^r by product a_j. This is the precise operation-sensitive version of the prime/composite information picture. At a nonzero row of a prime table, kernel size is1; at composite modulus a zero-divisor row has kernel size gcd(a,m). Statements about bits mean log2(d) Hartley ambiguity for this observable, not a Shannon entropy assertion without a probability law.

**Signed-central-table scope.** A central-table encoding can propose epsilon=±1, but orientation, k_j, and divisibility/valuation conditions are not recovered from the reduced output. The theorem takes these as retained word data. For a genuine integer orbit restrict inputs to its admissible 2-adic cylinder; the displayed full-residue fibre count is not automatically the number of admissible positive integer predecessors. It neither proves descent nor constructs an infinite ordinary Collatz trajectory.

## 4. The concrete cross-frontier transfer and boundary

The common object is an **operation-indexed finite observer**, not a claimed equality of Collatz and tournament trajectories.

- Source: an actual comparable rooted tournament pair with a legal birth word. Map: (pair, word)↦(Boolean support of q-difference, q-residue). Preserves exact final H equality and complete modular q response. Destroys graph adjacency and integer magnitudes. Sidecar: initial graphs plus word. Decisive test: the literal 401 repair and its O* language.
- Source: a fixed signed Collatz affine word. Map: input residue↦output residue. Preserves that word's affine response. Destroys exactly a kernel of size3^min(r,v3(m)). Sidecar: kernel index plus sign/exponent/orientation data. Decisive test: m=9 loses memory while m=25 does not.
- Source: Session254's neutral padding and the birth-family atlas. Proposed use: attach their complete modular interfaces or verified finite-state sidecars to a native repair automaton. Preserved predicate is only the chosen interface equality; marked endpoints must remain in the interface. The existing 401 four-coordinate carrier does not contain the full substitution jet, so it cannot simply be declared a universal substitution repair.

This last distinction is operationally important. Finite-state completeness for a **fixed observer** is a theorem; completeness for every useful observer, every source graph, or every odd target remains a separate mathematical obligation. The compiler is a new finite search instrument with honest equality gates, not a scalar coverage argument.

## Reproduction

    python3 -B 05-knowledge/results/session255_modular_collatz_20261111/cross_frontiers/probe.py

The exact script checks all16 support starts through word length8 (157456 state instances), prints the complete transition and minimized-state partitions, performs finite product searches on the literal cheaper401 response for m=3,5,7,9,11, and checks the affine fibre/sidecar formulas for every odd3<=m<50, 0<=r<=4, and K in {r,2r+1} (240 cases). The theorem proofs handle all word lengths and moduli; the bounded controls test only the new claims, not inherited conjectures.
