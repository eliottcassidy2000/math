# Observer compression defects and native Jacobsthal addresses

Session44, 23 October2026. Status: **PROVED** for the finite weighted-incidence identities, all-order converse statement, and native phase-transfer identities below; **FINITE-EXACT** for the independently reconstructed order-four/order-five tables and their exact algebraic spectra. The computations are certificates, not an all-order native spectral classification. No external novelty claim is made.

## What is new relative to the incoming packet

The incoming [Farey/conductance theorem](../../native_guard_moves_farey_20261023/farey/FAREY_MASSES_AND_NATIVE_CONDUCTANCE.md) already explains `7/15` and `8/15` by the same raw conductor `4/5`, masses `(3,4)` versus `(2,6)`, and the abstract Berggren tree of rational two-state gaps. We do not claim those results anew.

This packet adds:

1. An exact trace-defect criterion for deciding when a source compression loses only zero modes, with the complete order-four eigenspaces as its native example.
2. The first higher-order hostile: both unpaired order-five gaps have algebraic degree seven, and the compressed-source gap is strictly smaller than the full-source gap. Complement-pairing compressed sources produces a degree-five gap, not a new complementary rational pair.
3. An all-order converse-merging invariance on the reflection-invariant fixed-path carrier, proved by equality of target columns. It is not a claim about untransported pins or arbitrary guards.
4. A genuine native fraction mechanism: fixed-endpoint cyclic-cell continuation probabilities follow the Jacobsthal denominators, a sparse subsequence of the inherited Berggren parabolic spine. Their generalized continued-fraction transfer has determinant `-2`, not the determinant-one Farey property. Neither `7/15` nor `8/15` occurs as one of these unguarded continuation probabilities.

The inheritance board is **incidence mass / conditional expectation / kernel contrasts / literal phase / Farey determinant / actual clock**. The anchor is faithful observer compression; the niche is the exact native continuation law; the wildcard is whether its continued fraction can support a Farey transfer. The hostile tests below settle the first nontrivial cases rather than comparing decimal appearances.

## 1. Exact zero-mode compression theorem

Let `W` be any nonnegative finite incidence matrix with strictly positive row masses `r_i` and column masses `c_t`. Put

    K = D_r^(-1) W D_c^(-1) W^T.

Partition its source rows into blocks `B`, and sum rows within each block to form `Wbar`; let `r_B=sum_(i in B) r_i` and `Kbar=D_rbar^(-1) Wbar D_c^(-1) Wbar^T`. Define normalized target profiles

    p_it=W_it/r_i,       p_Bt=Wbar_Bt/r_B.

Let `J` lift a block function to the fine rows, and let `E` take the `r`-weighted average within each block. Then `EJ=I` and, always,

    Kbar = E K J,                                             (1)
    tr(K)-tr(Kbar)
      = sum_B sum_(i in B) r_i sum_t (p_it-p_Bt)^2/c_t.         (2)

Consequently the following are equivalent:

- Every normalized target profile is constant within its source block.
- The right side of (2) is zero.
- `K=J Kbar E`, so the weighted within-block contrast space `ker E` is killed and the nonzero spectrum is unchanged, with multiplicities.
- `K` and `Kbar` have the same nonzero eigenvalues, with multiplicities.

Proof. Equation (1) is direct multiplication. Expanding the squares in (2) and using `sum_(i in B) r_i p_it=r_B p_Bt` gives the trace formula. If profiles agree, every incidence row is its row mass times the common profile, giving the factorization. That factorization splits the source space into `im J` and `ker E`, with operators `Kbar` and zero. Conversely equal nonzero spectra imply equal traces, and every summand in (2) is nonnegative, forcing profile equality. QED.

This is stronger than merely proving a quotient is lumpable: it tests whether *all* forgotten modes are zero. Equal rank alone does not suffice. No assertion that `K-JKbarE` is positive semidefinite is made. With incomparable source partitions, apply (2) to their common refinement; there is no direct coarsening inequality between them.

The dual statement applies to merging target columns: identical normalized source profiles ensure the source operator is unchanged. This uses the entire weighted profile, not a coincidence of support or column sums.

## 2. Full versus compressed order-four eigenspaces

Use the native Hamiltonian path `0 -> 1 -> 2 -> 3`, and free backward chords `(0,2),(0,3),(1,3)`. The full source is the ordinary graph of those actual chords on four vertices; the compressed source replaces `(i,j+1)` by `(i,j)` on three vertices. We delete only the zero/transitive word.

In full-source order `(single edge, P3, two disjoint edges, P4)`, the masses are `(3,2,1,1)` and

    K_full = [[11/15,2/15,1/15,1/15],
              [1/5,  2/5, 1/5, 1/5 ],
              [1/5,  2/5, 1/5, 1/5 ],
              [1/5,  2/5, 1/5, 1/5 ]].

The compressed map merges `P3` and the two-disjoint-edge class. Both normalized target profiles are the point mass at the strong four-vertex tournament. Thus (2) is zero, and the compressed masses become `(3,3,1)` without changing the positive spectrum.

The full eigenspaces are exactly

    eigenvalue 1:     span{(1,1,1,1)},
    eigenvalue 8/15:  span{(4,-3,-3,-3)},
    eigenvalue 0:     {f: f_1=0, 2f_2+f_3+f_4=0}.

A zero-mode basis is `(0,1,-2,0)` and `(0,0,1,-1)`. The first is the specific forgotten contrast. The compressed zero space is `span{(0,1,-3)}`. These are right eigenspaces, with orthogonality measured by the retained row masses, not by uniform class measure.

The further merge of all strong-only source rows leaves masses `(3,4)` and the inherited two-state chain with gap `7/15`. Keeping all eight words and complement-pairing the compressed source instead gives masses `(2,6)`, gap `8/15`, and nonconstant eigenvector `(3,-1)` with eigenvalue `7/15`. These are different weighted carriers. Their complementary gap values do not establish a common full spectrum.

## 3. Exact order-five irrationality and the first failure of compression

Now use the same definitions at native order five, with all `2^6-1=63` nonzero path words. Canonical graph and tournament labels are the minimum lexicographic edge-bit masks under all vertex permutations. The compressed free-word code uses pairs `01,02,03,12,13,23`; the native code uses pairs `01,02,03,04,12,13,14,23,24,34`.

The source partitions are genuinely incomparable. Here are explicit witnesses from the complete enumeration:

- Free words `3` and `33` have the same full graph label `3`, but compressed labels `3` and `12`.
- Free words `3` and `9` have the same compressed graph label `3`, but full labels `3` and `20`.

The exact dimensions and traces are:

| carrier | rows | rank of incidence | trace of coupling |
|---|---:|---:|---:|
| full graph |17|8|2665427/675675|
| compressed graph |10|8|22433/6435|
| joint refinement |22|8|691/143|

The refinement-to-full trace defect is `599548/675675`; the refinement-to-compressed defect is `8662/6435`. Both are strictly positive. In particular, equal incidence rank eight does not make either forgotten contrast space invisible. The full-minus-compressed trace difference is `309962/675675`.

For `lambda` the spectral variable, the full characteristic polynomial is `lambda^9(lambda-1)F(lambda)/34479695250`, where

    F = 34479695250 lambda^7 -101537044560 lambda^6
        +124233482025 lambda^5 -81343863297 lambda^4
        +30498542418 lambda^3 -6453737759 lambda^2
        +695751836 lambda -27925620.

The compressed characteristic polynomial is `lambda^2(lambda-1)C(lambda)/3377602800`, where

    C = 3377602800 lambda^7 -8397030240 lambda^6
        +8429380020 lambda^5 -4360639860 lambda^4
        +1225690398 lambda^3 -180039440 lambda^2
        +12012273 lambda -267391.

Both degree-seven factors are irreducible over the rationals. The checker supplies finite-field certificates: the monic reductions are

    F mod31: x^7-14x^6+6x^5-6x^4+8x^3+3x-12,
    C mod23: x^7-2x^6-10x^5+4x^4+9x^3+5x^2-5x+2.

For each reduction of prime degree seven, it checks `x^(p^7)=x mod f` and `gcd(f,x^p-x)=1`. These two identities imply irreducibility: every irreducible factor has degree dividing seven, and no linear factor is possible. Nonvanishing leading coefficients and Gauss's lemma transfer irreducibility to the rational factors.

Exact rational Sturm counts locate the largest residual eigenvalues:

    635789/10^6 < lambda_full < 635790/10^6,
    642470/10^6 < lambda_compressed < 642471/10^6.

All other roots lie below these intervals and all seven roots lie in `(0,1)`. Thus both gaps have algebraic degree seven, and

    gap_compressed < gap_full.

Calling one source "compressed" does not make its class partition a coarsening of the other. The reversed gap inequality is therefore not a contradiction to coarsening monotonicity.

For completeness, retain all64 words and complement-pair the compressed sources. There are six source rows with masses `(2,12,24,8,6,12)`. The characteristic polynomial is `(lambda-1)P(lambda)/562933800`, where

    P = 562933800 lambda^5 -914443020 lambda^4
        +553227084 lambda^3 -153488811 lambda^2
        +19184700 lambda -874496.

The same irreducibility test works modulo19. Its largest residual root lies in `(573377/10^6,573378/10^6)`, so its gap has degree five. Target converse-merging leaves this operator unchanged by the next theorem.

These are explicit counterexamples to extending the order-four rational gap picture to general native observer incidence. The inherited Berggren tree correctly parameterizes *rational two-state* gaps; it cannot contain these higher-degree exact gaps. Rational matrix entries, positive semidefiniteness, and common incidence rank do not imply rational eigenvalues.

## 4. All-order converse-merging invariance, with its guard boundary

Let `rho(i)=n-1-i`. Reverse a native path word's vertex order and reverse every tournament arc. The result still has the same forward Hamiltonian path, and its backward chord graph is the ordinary relabeling by `rho`. On compressed positions the same operation is ordinary reflection `i -> n-2-i`. Therefore both full and compressed source classes are unchanged, while the native tournament class changes to its converse.

For every source class `g`, this is a bijection between words counted by `W_(g,t)` and `W_(g,t^op)`. Hence these columns are identical, with equal masses. Merging a converse pair replaces two contributions `ww^T/c` by `(2w)(2w)^T/(2c)`, the same matrix. Thus native target converse-merging preserves the entire source coupling at every order, for either source observation and for any further source coarsening.

The proof also works with weights or a restricted word set that are invariant under this *actual reflected-converse map*. A set with indexed ports fixed in the printed positions, or an arbitrary electrical phase filter, need not be invariant. Such a guard requires its own transported-coordinate check. The theorem does not say that complementing the free bits is converse, nor does it connect states in a fixed-source native metagraph.

## 5. Native ternary transfer yields Jacobsthal, not arbitrary Farey fractions

Use the actual regular cyclic-cell host of the inherited [native branch-width theorem](../../session43_observers_phase_guards_20261023/native-branch-width/THEOREM.md). Each cell is a directed `C3`; successive cell phases satisfy `a_(i+1)-a_i in {0,1}`. The exact state is the shared current phase, and an initial actual owner retains the common-phase lift.

Let `A_r(delta)` count binary difference words of length `r` whose integer weight is `delta mod3`. Equivalently it counts the native paths with fixed initial phase and endpoint phase difference `delta`. The transfer matrix is `M=I+Q`, where `Q` cyclically shifts the three phases. The inherited formula is

    A_r(delta)=[2^r+c_(r-2delta mod6)]/3,
    c=(2,1,-1,-2,-1,1).                                  (3)

A complete odd-cell regular host has an even number of differences. Every finite remaining length r used below occurs as a continuation suffix in a sufficiently large such host, so r is not restricted to be even.

The new observable is the exact first-step conditional probability, under the uniform law on these accepted native continuations:

    Pr(d_0=0 | r,delta)= A_(r-1)(delta)/A_r(delta),
    Pr(d_0=1 | r,delta)= A_(r-1)(delta-1)/A_r(delta).       (4)

The two numerator counts differ by at most one. For `r>=3`, both are positive. If they agree, the probability pair is `(1/2,1/2)`. Otherwise, put

    J_r=(2^r-(-1)^r)/3.

Then the pair in (4), in one order, is exactly

    ((J_r-1)/(2J_r), (J_r+1)/(2J_r)).                     (5)

Proof. Insert the six periodic correction values in (3). For each parity of `r`, one endpoint residue has equal child counts; the other two have consecutive counts. Their sum is the odd integer `J_r`. After canceling the displayed common factor two in each fraction, the numerator counts are consecutive and coprime, so each fraction is reduced. QED.

Thus the unequal pairs begin `1/3,2/3`; `2/5,3/5`; `5/11,6/11`; `10/21,11/21`; `21/43,22/43`. Each fraction is a Farey neighbor of `1/2`, but the two fractions of the pair have determinant magnitude `J_r`, not one. They are not Farey neighbors of each other.

Neither `7/15` nor `8/15` occurs in (4): a reduced unequal pair would require `J_r=15`, equivalently `2^r=45+(-1)^r`, impossible; the equal pair is `1/2`. The occurrence of those fractions in the order-four *spectral* calculation is not the first-step probability law of this native three-state transfer.

There is an exact inherited-tree address. For `r>=3` put `k_r=(J_r-1)/2`. The smaller probability is `k_r/(2k_r+1)`, so in the spectral-coordinate convention of the incoming Farey packet it has Berggren word

    A^(k_r-1),         depth=(J_r-3)/2.                  (6)

The actual remaining native clock is `r`; the rational-tree depth grows exponentially in `r`. In the canonical ternary heap convention `h(wA)=3h(w)+1`, the address is `(3^((J_r-3)/2)-1)/2`, which grows doubly exponentially in `r`. This is an explicit clock-distortion theorem, not a native execution-time bound. The source, continuation probability and address are retained; the address alone loses the native path, endpoint owners and integer displacement.

### A genuine continued fraction, but determinant two

The Jacobsthal counts obey

    J_0=0, J_1=1, J_(r+2)=J_(r+1)+2J_r.

Consecutive positive terms are coprime. Their ratios satisfy

    R_r=J_(r+1)/J_r,     R_1=1,
    R_(r+1)=1+2/R_r.                                    (7)

This is a finite generalized continued fraction with numerator2 at each stage. Its transfer matrix is `[[1,2],[1,0]]`, of determinant `-2`. The exact determinant identity is

    J_(r+1)J_(r-1)-J_r^2=(-1)^r 2^(r-1).                (8)

It follows by induction from the recurrence, with `r=1` as base case. Hence consecutive ratio approximants cease to be Farey neighbors once `r>=2`. Equations(7),(8), rather than an unsupported Fibonacci analogy, explain the arithmetic of this native transfer.

### Retain the integer lift and the modulus

The polynomial matrix `M(z)=I+zQ` gives

    (M(z)^r)_(alpha,beta)
       = sum_(k=beta-alpha mod3) binom(r,k) z^k,
    det M(z)=1+z^3,       (M(z)-I)^3=z^3 I.              (9)

The exponent is the integer displacement erased by the phase residue. The inherited legal `000 <-> 111` coupled move changes this exponent by3 while retaining the endpoint phase; it is not an invisible kernel lift. Equation(9) is a genuine ternary return identity and does not authenticate an arbitrary electrical guard as a weight predicate.

All of these counts can be computed by integer transfer modulo any modulus. In the closed formula(3), division by3 is valid modulo8 because3 is a unit. To recover the answer modulo9 from the numerator alone, retain that numerator modulo27 before dividing, not only modulo9. For example the numerator at `(r,delta)=(2,0)` is3 and at `(6,0)` is66: both are3 mod9, but the counts are1 and22, distinct mod9. The lack of invertible third-root Fourier averaging modulo9 does not obstruct integer transfer; it forbids the lossy division step. This says nothing about the uncomputed dependence of an original mod8 electrical phase.

## 6. Reproduction and finite scope

Run:

    python3 -B check.py --write
    python3 -B check.py
    python3 -O -B check.py

Only the first command replaces the frozen receipt. The checker uses Python's standard library and exact integers/Fractions only, with explicit exceptions rather than optimization-sensitive assertions. It reconstructs every order-four/order-five word from its literal arcs, computes canonical source/target classes, independently verifies Hamiltonian counts and automorphism normalization, forms all incidence matrices, and computes characteristic polynomials by Newton trace identities. Finite-field polynomial arithmetic and rational Sturm sequences certify irreducibility and isolations independently of an external symbolic algebra package. All source tables, eigenspace equations, trace defects, converse counts and native transfer controls are included in `receipt.json`.

`SOURCE_PINS.json` records inherited source paths and immutable blob hashes. Its sources do not prove the new degree-seven/degree-five spectral certificates or conditional Jacobsthal law. General native spectral gaps, guarded source connectivity, and an all-size electrical-guard factorization remain OPEN.
