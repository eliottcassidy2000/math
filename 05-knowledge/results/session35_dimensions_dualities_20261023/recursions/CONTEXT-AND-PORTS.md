# Second-round theorem: contextual spectral equivalence, and the squarefree repair

Status: PROVED elementary algebra/combinatorics; witness FINITE-EXACT. This addendum does not change the previously audited RESULT.md.

## 1. Equality persists in every E/D context but Hamiltonian counts do not

Define two same-size tournaments to be marked-spectrally equivalent when their pairs (p_M,F_M) are equal. The formulas in RESULT.md prove:

**Context theorem.** If T and U are marked-spectrally equivalent, then W(T) and W(U) are marked-spectrally equivalent for every finite word W in the source–sink extension E and skew-Sylvester doubling D.

Proof: induction on word length; each operation updates (n,p,F) by the same deterministic formulas, and equal input sizes stay equal. Consequently every observable that factors through this pair after ANY such context is also equal. This is stronger than equality of a few moments, but weaker than tournament equivalence.

The root's six-vertex witness was independently verified. With bit j recording an outgoing arc to vertex j, outneighbor rows are

T: [58,4,9,18,38,14],
U: [22,12,8,17,38,15].

Both have
p(z)=z^6+15z^4+47z^2+1,
q(z)=1^T adj(zI−M)1=2z(3z^4+34z^2+59),
F=q/p.

But H(T)=29 and H(U)=25 by independent subset/end-vertex dynamic programming. E preserves H, so H(E^k T)=29 and H(E^k U)=25 for every k≥0. Thus even infinitely many E/D-context spectral observations cannot recover H universally. We make NO claim that their Hamiltonian counts remain different after every word containing D; that is unnecessary and not established here.

## 2. Exact path preservation requires exclusion memory

Let A be the 0/1 adjacency matrix, and work over the commutative square-zero ring
R_n=Q[x_1,...,x_n]/(x_1^2,...,x_n^2).
Write X=diag(x_i). The matrix AX is nilpotent, so

P_A(X)=X(I−AX)^(-1)

is a finite polynomial matrix. Its (i,j) entry enumerates all directed vertex-simple paths from i to j, weighted by the product of visited vertex variables, with the trivial one-vertex path included. Indeed each matrix-product monomial corresponds to a walk; revisiting a vertex supplies x_v^2=0 and kills it. All surviving coefficients are positive counts, so this is exclusion rather than cancellation.

Define
Φ_A(X)=1+1^T P_A(X)1.

Then for every nonempty subset S,
[x_S]Φ_A=H(T[S]).
In particular [x_1...x_n]Φ_A=H(T). This is an exact path-invariant lift of an ordinary resolvent, with a transparent preserved predicate and no orientation-blind moment truncation.

For the canonical +2 extension, with new source a and sink b,

**Factorization theorem:** Φ_{E(T)}=(1+x_a)(1+x_b)Φ_T.

Proof: a simple path involving old vertices has an arbitrary old simple path, optionally preceded by a and optionally followed by b. The cases containing no old vertices are exactly empty path, a, b, a→b, supplied by the factor times the constant 1. Equality holds for every subset coefficient, not just the full Hamiltonian coefficient.

For the collision, Φ_T and Φ_U already differ at degree three. On vertices {0,1,2}, counts are 3 and 1. Overall numbers of differing labelled subset coefficients are 8 at degree3, 9 at degree4, 3 at degree5, and 1 at degree6. Thus the lost coordinate can be seen as WHICH labels participate in cyclic path structure, not merely a missing high-order scalar moment.

## 3. Exact doubling via a Schur port calculation

For skew-Sylvester doubling, the 0/1 adjacency is exactly

A_D = [[A,A+I],[A,A^T]].

Introduce independent vertex variables in the two sheets, X_0 and X_1; do NOT identify twin variables. The block resolvent is

I−A_D diag(X_0,X_1)
 = [[I−A X_0, −(A+I)X_1],
    [−A X_0, I−A^T X_1]].

Eliminating the second block gives the exact Schur complement

S_0=I−A X_0−(A+I)X_1(I−A^T X_1)^(-1)A X_0.

All inverses exist because their constant term is I and the positive-degree ideal is nilpotent. Standard block inversion reconstructs the full port matrix and hence Φ_D and its Hamiltonian coefficient. This is a precise boundary-composition rule that retains vertex exclusion.

Why twins need separate labels: a doubled Hamiltonian path visits BOTH copies of each old vertex. Identifying x_{i,0}=x_{i,1}=x_i in a square-zero ring would erase every such path. The multiplicative mode therefore doubles not only size but exclusion-memory labels.

Why full ports retain orientation: for i≠j, the degree-two coefficient of (P_A)_{ij} is A_ij x_i x_j. Thus the full matrix is lossless for labelled adjacency. The scalar Φ is not claimed lossless (for example it is invariant under global arc reversal).

## 4. Efficiency and scope

The exact lift does not deliver a free polynomial-time Hamiltonian invariant. R_n has dimension 2^n; usual subset/end-vertex dynamic programming uses O(n2^n) stored values and O(n^2 2^n) arithmetic. Schur port algebra may reorganize structured computations, but naïve multiplication over the squarefree ring can be more expensive than the standard DP. Any claimed compression must retain the exclusion boundary information and needs its own proof.

The concrete unification is now three levels:
- scalar marked spectral pair: very small, closes exactly under E/D, but provably blind to H;
- labelled squarefree scalar path polynomial: exactly encodes all induced Hamiltonian counts;
- full squarefree port matrix: supports exact Schur composition and retains orientation.

The first level's failure is rigorous even under unlimited E/D probing. The second and third levels repair the relevant lost information at an explicit exponential cost. A useful next research target is a restricted tournament family where ports admit a small exact representation stable under both operations, rather than another scalar spectral invariant.

## Verification

Run python3 session35-recursions/collision.py. Output: collision.out.
- Independent exact characteristic and marked-adjugate polynomial calculation for both six-vertex matrices.
- Independent subset/end-vertex Hamiltonian count: 29 versus25.
- All 256 subset coefficients of each eight-vertex E-extension checked against the factorization.
- No claim of polynomial-time complexity, minimal carrier, or full tournament reconstruction from scalar Φ.

## 5. Important identification: the marked pair is exactly the adjacency spectrum

This closes the inheritance loop and prevents overstating the new carrier. Since 2A=M+J−I, define

g_A(z)=det(zI−(M+J))=2^n p_A((z−1)/2).

The determinant lemma gives g_A=p_M−q_M, where q_M=p_M F_M. Skew symmetry implies p_M(−z)=(-1)^n p_M(z) and q_M(−z)=(-1)^(n−1)q_M(z). Therefore

p_M(z)=[g_A(z)+(-1)^n g_A(−z)]/2,
q_M(z)=[(-1)^n g_A(−z)−g_A(z)]/2.

Hence for same-size tournaments:

(p_M,F_M)=(p_N,F_N) iff A_T and A_U are cospectral.

The all-word result can therefore be stated cleanly: adjacency cospectrality is a congruence under canonical source–sink extension and skew-Sylvester doubling. The skew spectrum alone lacks the boundary mark, but the ordinary adjacency spectrum already encodes it through this parity decomposition.

This identifies the correct predecessor: `01-canon/theorems/THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md`, with its cited THM-499 size-six spectral boundary, already separates adjacency spectrum from disjoint odd-cycle packing. No novelty claim is made for non-spectrality of H itself. The new synthesis is its exact preservation under the two selected recursive operations.

For the current witness, direct cycle enumeration gives
T: c3=6, c5=6, D33=1;
U: c3=6, c5=6, D33=0.
At six vertices OCF reads H=1+2(c3+c5)+4D33. Thus the entire 29−25 discrepancy is exactly one missing disjoint triangle pair, weighted by 2². The individual odd-cycle totals are not the missing coordinate; their compatibility is.

## 6. Determinant identity recovers known OCF, not a new theorem claim

In the same square-zero ring, the matrix determinant lemma proves

Φ_A=det(I+(J−A)X)/det(I−AX).

For a tournament J−A=I+A^T. Since X²=0,
I+(I+A^T)X=(I+A^T X)(I+X).
Using det(I+A^T X)=det(I+AX), we obtain

Φ_A=prod_i(1+x_i) det(I+AX)/det(I−AX).

Finite logarithms are valid in the nilpotent positive-degree ideal over Q:

log Φ_A = sum_i x_i + 2 sum_{k≥3 odd} Tr((AX)^k)/k
         = sum_i x_i + 2 sum_{directed simple odd cycles C} x_{V(C)}.

Repeated-vertex terms vanish. Each oriented simple cycle has k rooted representatives, accounting for division by k. Exponentiation selects vertex-disjoint cycles, each with weight2, together with optional singleton vertices. Extracting the full coefficient gives precisely the known odd-cycle collection formula H=I(Ω(T),2).

This is a transparent rederivation of the existing Rédei/Grinberg–Stanley/Irving–Omar mechanism, not a claimed new OCF. It explains how a vertex-marked resolvent differs from the unmarked adjacency spectrum: the square-zero algebra retains disjointness before taking the trace/determinant, while an ordinary scalar spectral projection loses it.

Primary literature checked:
- Grinberg and Stanley, The Rédei–Berge symmetric function of a directed graph, https://arxiv.org/abs/2307.05569 .
- Irving and Omar, Revisiting The Rédei–Berge Symmetric Functions via Matrix Algebra, https://arxiv.org/abs/2412.10572 (v3, October23 2025).
Repository predecessor: `01-canon/theorems/THM-002-ocf.md`, whose all-n external proof routes through Irving–Omar Corollary20 and Grinberg–Stanley Theorem1.39/Lemma6.5. The identities above are proved here directly to avoid reliance on version-sensitive theorem numbering.

## 7. Root scan provenance and exact scope

Root artifacts, inspected but not independently rerun in full here:
`session35-root/scan_marked_spectra.py`, `session35-root/marked_spectrum_scan.json`.

For n=3,4,5,6,7, the enumerated fixed-Hamiltonian-path presentation counts are respectively 2,8,64,1024,32768. Distinct marked-spectral classes: 2,3,9,28,168. H-ambiguous classes: 0,0,0,3,47. These are presentation scans, not counts of isomorphism classes. The class key uses traces Tr(M^k) and all-ones moments 1^T M^k1 for k=1..n; Newton identities plus the resolvent recurrence recover p and q, so the key is exact. n≤7 and k≤7 make int64 overflow impossible (a crude bound on each entry of M^k is n^(k−1), and on the all-entry sum n^(k+1), far below 2^63).

Every tournament admits a Hamiltonian path, so after relabelling the fixed-path presentations meet every tournament isomorphism class. Therefore the negative n≤5 result covers all tournaments of those orders, not merely random examples. The independent six-vertex witness verifies the first failure directly. The all-word conclusion is proved algebraically and does not extrapolate the finite census.
