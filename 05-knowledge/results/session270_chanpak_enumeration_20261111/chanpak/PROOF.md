# A million-target compiler and averaged-core path-cover tomography

**Status:** PROVED transfer, all-orders moment transform and exact linear-rank certificate; FINITE-EXACT million-target coverage. Chan–Pak remains OPEN. No new-family or literature-priority claim is made.

## Inheritance and actual advance

The inherited full-cover substitution mechanism is Session268 `spectrum/PROOF.md`, with hidden state q_r, sharp homogeneous-probe cost 2r−2, and sparse fixed exponential atom surfaces. The relevant older niche is **THM-4878, `01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md`**: its G(2,0) is exactly the two-apex family used here, its finite-state principle is known, its 200-state nonnegative realization and empirical Hankel rank17 are inherited. That source certified coverage through2047, not an all-target theorem. No novelty is claimed for the family, rational-series method, or number17.

What this packet adds is an explicit36-state endpoint realization, a rational17-state realization with finite *invariant-subspace* certificate valid for all words, a reusable direct-address compiler for499,998 targets, and a general fixed-core average identity that converts arbitrary path-cover data into a small rational ordinary GF. The minimal-rank certificate upgrades the inherited numerical rank observation for this language to an exact all-word proof. The average identity supplies a different physical tomography tradeoff, not another audit of prior theorems.

Concept board: returning endpoint state; core-run allocations; weighted rational languages; target range versus moments; homogeneous versus averaged sensors. The hostile is treating a scalar H or an average of H as enough to reconstruct exact range. The corrected near miss is iterating scalar seam formulas without retaining root data. The present state does return under every letter.

## 1. Physical language and a36-state nonnegative representation

Fix a labelled tournament C on k≥1 core vertices. Add an ordered transitive spine v_1→...→v_m. A letter c_i in{0,...,2^k−1} records cross edges: bit j is1 iff core vertex j→v_i. Core and ordered spine are designated. Distinct words are distinct labelled presentations; they need not be distinct unmarked isomorphism classes. For k=2 fix core0→1, so the alphabet is{0,1,2,3}.

In any Hamilton path, order the core vertices by their occurrence, a permutation p. They divide the spine vertices into k+1 slots. Within each slot the vertices appear in increasing spine order. Conversely, assigning each spine vertex to a slot determines at most one path for a fixed p.

Scan the spine left to right. Each of slots0,...,k−1 has state0(empty),1(nonempty, latest vertex beats the next core vertex), or2(nonempty, latest vertex loses to that core vertex). Slot k has state0 or1, since it has no next core vertex. When a first vertex enters a slot j>0, require that the previous core vertex p[j−1] beats it. A later vertex can enter regardless of its color relative to that previous core. Update the slot status according to its relation with the next core vertex. Each letter assigns its vertex to exactly one slot, so each old state has at most k+1 transitions, counted with multiplicity.

Initially all slots are empty, once for each p. Finally reject status2 at any slot. An empty internal slot j=1,...,k−1 is accepted only when C has p[j−1]→p[j]; empty exterior slots are legal. These are exactly the missing core-core adjacencies. Thus accepted runs are in bijection with directed Hamilton paths.

The resulting state count is D=k!·2·3^k, giving D=36 for k=2. With column vectors, v(wc)=M_c v(w), and H(w)=a^T v(w). All M_c have nonnegative integer entries. The full state, not merely H, is preserved under arbitrary future letters. `spine.py` constructs the matrices directly from these rules. Reading one word costs O(mkD) integer additions and O(D) stored integers, excluding integer bit cost; writing the tournament uses O((m+k)^2) bits.

This is a useful specialization of inherited frontier dynamic programming, not a claim to a new algorithmic paradigm.

## 2. Exact17-dimensional quotient and its lower bound

`rank.py` starts with the36-state initial vector and closes its rational span under the four matrices. Independently it closes the acceptance vector under the transposed matrices. It produces reachable dimension28 and observable dimension20. Each closure is explicitly checked invariant under *every* letter, rather than inferred from stabilization of a sampled Hankel table.

For reachable basis B and observable basis C, the pairing C^T B has rank17. Its stored17×17 minor has determinant150994944≠0. Every Hankel entry H(uv) factors through these invariant spaces, giving an upper bound17 for the entire infinite word series; the nonzero minor gives the matching lower bound. Thus the minimum rational linear representation has dimension17. This does not say the smallest nonnegative integer automaton has17 states.

The program selects17 reachable columns V and17 observable rows W, sets K=W^T V, and constructs R_c=K^{-1}W^T M_c V. It checks C^T(M_c V−V R_c)=0 for each c and the corresponding initial-vector identity. This proves that the exported initial vector, output row and four rational17×17 matrices in `rank.json` evaluate *all words*, not just the finite controls. Negative rational coordinates are allowed in this quotient; the36-state positive realization remains the transparent compiler.

## 3. General fixed-core mass transform

Let q_r(C) be the number of ordered lists of r nonempty disjoint directed paths covering the k core vertices. In particular q_1=H(C) and q_k=k!. Define

A_m(C)=sum over all 2^{km} spine color words w of length m of H(C,w).

This is a sum over labelled cross-edge choices with the core fixed, not over isomorphism classes. Set onto(m,t)=t! S(m,t), with onto(0,0)=1 and onto(m,0)=0 for m>0. Then

A_m(C)=2^{km} sum_{r=1}^k q_r(C) sum_{b=0}^2 binom(2,b) onto(m,r−1+b)/2^{2r−2+b}.   (1)

**Proof.** Delete the spine from a path but retain the maximal consecutive core runs. These are an ordered path cover of C, with r runs. The r−1 internal spine slots must be nonempty; exactly b of the two exterior slots are nonempty. There are onto(m,r−1+b) assignments of ordered spine vertices onto those named slots, and each assignment fixes the increasing order in each slot. Exactly2r−2+b cross edges are required. Over uniform independent core-spine orientations, the probability of these constraints is2^{−(2r−2+b)}. The required edges are distinct, even for a singleton internal slot, because the adjacent core endpoints are distinct. Multiply the expectation by2^{km} and sum all r and boundary choices. Every path has a unique decomposition, so there is no overcount. ∎

With y=exp(2^k z), the exponential generating function is

sum_{m≥0} A_m(C) z^m/m! = ((y+1)^2/4) · sum_{r=1}^k q_r(C)((y−1)/4)^{r−1}.   (2)

Indeed sum_m onto(m,t)x^m/m!=(e^x−1)^t. This is a polynomial in y of degree k+1. If its coefficients are d_0,...,d_{k+1}, the ordinary GF is

sum_{m≥0} A_m(C)z^m = sum_{j=0}^{k+1} d_j/(1−j2^k z).                (3)

The j=0 term is just the constant d_0; it contributes only at m=0. Consequently there is a constant-coefficient recurrence of order at most k+1 for m≥k+2, with characteristic polynomial product_{j=1}^{k+1}(X−j2^k). Some poles can cancel. Equivalently include the zero root and use an order k+2 recurrence starting at m=k+2. These are arithmetic-operation bounds for fixed core order, not constant bit complexity for growing exact integers.

### Triangular tomography and its honest cost

For m=r−1, the coefficient of q_r in(1) is

m!·2^{(k−2)m}.

Terms q_j with j>m+1 vanish. Thus A_0,...,A_{k−1} recover the entire cover vector by triangular elimination. Two same-order cores first differing at q_r have identical earlier masses, and their mass difference at spine length r−1 is exactly(r−1)!2^{(k−2)(r−1)} Δq_r.

This uses only r−1 additional vertices *per probe*, compared with2r−2 in Session268's single homogeneous-module context. It also sums2^{k(r−1)} different orientations, so it is not a free query improvement and does not violate the earlier module lower bound. Nonhomogeneous cross edges intentionally break that hypothesis. The m≥r−1 requirement is sharp for this particular averaged-spine observation because a path with m spine vertices has at most m+1 core runs. No lower bound for arbitrary rooted, unaveraged observations is claimed.

### Specialization to the two-core compiler

For core0→1, q=(1,2), and(2) becomes(e^{4z}+1)^3/8. Therefore A_0=1 and

A_m=(12^m+3·8^m+3·4^m)/8 for m≥1.

The ordinary GF is

(1−18z+80z²−48z³)/((1−4z)(1−8z)(1−12z)).

The exact recurrence is A_m=24A_{m−1}−176A_{m−2}+384A_{m−3} for m≥4, with initial1,6,48,432. This supplies instantaneous fixed-core mass evaluation by exponentiation, not instantaneous arbitrary target inversion. Tensor powers of the36-state matrices similarly give rational GFs for every fixed integer moment sum_w H(w)^s, but their dimensions may grow as36^s and support is not recovered merely by quoting those moments.

## 4. Constructive finite coverage and witness costs

`census.cpp` exhaustively traverses every four-letter word at lengths0 through14, using iterative deepening. It records a shortest word for every positive odd target≤1,000,000. No isomorphism pruning, heuristic search, scalar-state merging, or product closure is used. For this run all entries fit uint64: a partial slot assignment is at most2·3^14<10^7, and every Hamilton count is bounded by that same count of assignments. The code should not be used at arbitrary larger depth without revisiting overflow.

**FINITE-EXACT result:**499,998 positive odd targets≤1,000,000 have recorded witnesses; the only absent targets are7 and21. Every such witness has at most16 vertices. All38,764 eligible odds below77,533 already have witnesses on at most13 vertices; all115,258 eligible odds below230,521 on at most14. The exact interval convention is strict upper endpoint in those last two statements. At15 vertices the first additional miss is776,161; at16 vertices none remains below the million cutoff.

Because all shorter words are exhausted first, the recorded sizes are genuinely minimal within this fixed two-apex language, not just minimal among a random sample. They are **not** asserted globally minimal among tournaments. For target1, the language returns its two-vertex core, although a singleton is globally smaller. The graph need not be strong; this is the ordinary Hamiltonian-path spectrum requested by Chan–Pak, not a strong-spectrum theorem. The full conjecture, all-target completeness of G(2,0), and an asymptotic O(log N) exact-target witness bound remain OPEN. A finite interval with16-vertex witnesses proves none of them.

`atlas.b85` stores one little-endian32-bit record for each positive odd index: high4 bits are spine length m; low28 bits are its base4 word, first letter most significant, padded with leading zeros to length m. FFFFFFFF means missing. Zlib compression plus Base85 wrapping preserves the packed2,000,000-byte array; the raw SHA256 is in `atlas.json`. After one decompression, target lookup is direct-address, and emitting its graph costs O((m+2)^2). This is a practical finite compiler, not a globally constant-time solution.

## 5. New-artifact controls and reproducibility

The all-word proof above justifies the transfer. The supplied tests target the new implementation:

* `controls.json`:2205 independent subset-Hamilton-DP comparisons, exhaustively over one-core words through spine7, two-core words through spine5, and cyclic three-core words through spine3.
* `atlas_checks.json`: replay of every one of the499,998 stored words through the transfer, plus independent subset DP on all attained targets≤1001 and deterministic evenly spaced atlas entries. The receipt gives the exact number.
* `moments.json`: direct enumeration of colorings for five fixed cores of orders1–4, triangular recovery, and the displayed two-core recurrence through m=19.
* `rank.json`: invariant-span closure and nonzero determinant certificate, with explicit17-state matrices.

Reproduce from this directory:

    python3 spine.py
    python3 moments.py
    python3 rank.py
    python3 compile_target.py 999999
    python3 compile_target.py --decode /tmp/s270-atlas.bin
    g++ -O3 -std=c++17 verify_atlas.cpp -o /tmp/s270-verify
    /tmp/s270-verify /tmp/s270-atlas.bin

Regenerate the exhaustive census and packed certificate:

    python3 spine.py --cpp > /tmp/s270-census.cpp
    g++ -O3 -std=c++17 /tmp/s270-census.cpp -o /tmp/s270-census
    /tmp/s270-census 14 1000000 > /tmp/s270-census14.json
    python3 pack.py /tmp/s270-census14.json

The census traverses477,218,583 prefix visits including iterative-deepening repeats. Target distribution is stored compactly; raw multi-megabyte JSON need not be published. Prior coverage figures are inherited as context, not reverified here. The next target-directed frontier is exact inversion of the17-state semigroup or a returning cone with a provable digit algorithm; a rational mass GF alone is insufficient.
