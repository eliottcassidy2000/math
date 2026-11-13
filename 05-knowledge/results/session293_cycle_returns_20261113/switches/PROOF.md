# Two physical port switches and a complete carry classification

Status: PROOF CANDIDATE pending independent review. This is a scoped family, not Chan–Pak.

## Inheritance and carrier

The one-apex formula is inherited from THM-4876 §6 (`01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md`). Let T(s) have transitive vertices 0,...,n−1 and one apex z, with z→i iff bit s_i=1. Then

H(T(s))=1+2 sum_{i<j} s_i(1−s_j)c(j−i),
c(1)=c(2)=1, c(d)=2^(d−2) for d≥2.

It follows either by counting all odd cycles through z, whose complements are transitive, or the inherited quadratic ear-cut law THM-4097/4114. Session292's sparse guard theorem supplies the physical interpretation: reordering adjacent vertices in the transitive interior preserves all interior cyclic triples but changes selected quadratic port observations. We now pay for two such changes exactly.

## 1. The full-count two-switch invoice

For a binary word w of length m≥2 compare

X(w)=T(10 w 01),   Y(w)=T(01 w 10).

On named vertices this is swapping two disjoint adjacent pairs of the interior total order while leaving apex incidences unchanged; after relabelling the order, the displayed words result. Each of the two interior edge reversals is a genuine tournament operation. No formal negative-weight graph is used.

Let B(w)=the binary integer w_0...w_{m−2}, and C(w)=the binary integer (1−w_{m−1})...(1−w_1). Then

**Full-count identity:** H(Y(w))−H(X(w))=2(B(w)−C(w)).

Proof. A single adjacent10→01 switch at positions i,i+1 changes H/2 by
−1−sum_{p≤i−2}s_p2^(i−p−2)−sum_{q≥i+3}(1−s_q)2^(q−i−3).
An opposite switch reverses this sign, evaluated at its actual current word. For the two separated switches the cross term in H/2 is −2^(m−1). At the initial word the left change is −1−sum_{r=1}^{m−1}(1−w_r)2^(r−1)−2^(m−1), and the right change is 1+2^m+sum_{r=0}^{m−2}w_r2^(m−r−2). Their sum with the cross term is B−C. This computes the full H, not merely the named opening sector. □

## 2. Exact zero and unit repairs, at every length

For m≥2, B=C iff w_r=1−w_{m−1−r} for every r. Thus zero repairs exist iff m is even, and there are exactly2^(m/2) such words: the anti-palindromes.

For m≥4,

B−C=1 iff w=1 0^(m−2) 1;
B−C=−1 iff w=0 1^(m−2) 0.

At m=2 the positive solution is11. At m=3 the positive solutions are011,101,110. Negative solutions are their bitwise complements.

Proof for +1. In binary addition C+1=B, let j be the index of the carry-stopping digit, indexing from the left. Overflow cannot occur because B has the same length. At r<j we have w_r=1−w_{m−1−r}; at r=j both w_j and w_{m−1−j} are1; at r>j the digits w_r and w_{m−1−r} are0. If j=0 this is the stated word. If j>0, zeros in the two suffixes force j≥(m−1)/2. If j>(m−1)/2, the opposite carry digit lies at an index less than j and contradicts the anti-pair relation there. Hence m is odd and j=(m−1)/2. If j≥2, positions1 and m−2 are both0 but must be opposite; impossible. The only alternative is m=3,j=1, giving011 or110. Bitwise complement negates B−C. □

The zero repairs are not hidden new cores: deleting X's final singleton sink and Y's initial singleton source leaves converse cores when w is anti-palindromic. This explains their equality by the inherited SCC/converse mechanisms. The theorem identifies the exact obstruction instead of claiming that symmetric cancellation is a new universal seam.

## 3. An explicit dyadic +2 family and its renewal limit

For m≥3 take w=1 0^(m−2)1. Direct substitution in the inherited formula gives

H(X)=15·2^(m−2)+1,
H(Y)=15·2^(m−2)+3.

Both are actual tournaments on m+5 vertices (one redundant exterior singleton may be removed). Indeed inversion weights of X sum to15·2^(m−3), proving the first formula, and the preceding invoice gives the second. For example m=3 gives31→33 and m=4 gives61→63.

Lengthening the zero run is a graph-defined rearming move on each side. At the level of full counts it sends

h_X→2h_X−1,   h_Y→2h_Y−3.

The two branches retain difference2 at every scale. The corresponding inverse construction decreases the target on these recognized words. This is a reproducible infinite pair of physical repairs, not a new assertion that every even gap occurs; the latter was already available from Sessions199/221/243/244. No claim of priority is made for these particular numerical values. The new ingredient is the complete all-length carry classification for this declared two-switch mechanism.

A second +2 repair cannot simply be repeated at the same two slots: the first move changes10...01 to01...10, so repeating the edge swaps returns to the first graph and subtracts2. Rearming by growth changes the absolute target according to the two affine laws. The rigid +2 language has only one word per m≥4, hence only O(log X) such values below X. This family cannot cover the odd spectrum. The target-preserving anti-palindrome language, while much larger, moves only the converse core and exterior singleton under these guards.

## 4. Why this advances the question

The concrete distinction is between (a) an exact scalar cancellation, (b) renewal of the port word, and (c) selection of an arbitrary absolute target. Here all three are separately visible. The binary carry is the missing sidecar in an apparent ±2 coincidence; requiring it at every length collapses a large word language to a single dyadic ray. A prospective five-core composition must escape this one-apex word class or vary additional response coordinates, while retaining an exact invoice like §1.

Next open test: allow a nontransitive substituted block at one of the two switched positions, derive its full path-cover response rather than replacing it by H alone, and ask whether the resulting repair language contains a gap-free target interval with renewable ports. No such uniform selector is proved here.

## 5. Nontransitive replacement pays an exact repair tariff

Replace the middle zero run in the preceding +2 family by ANY nonempty tournament C, inserted as a module: the first three transitive outer vertices beat C, C beats the last three outer vertices and the apex. Keep the outer apex words101|101 for X(C), and011|110 for Y(C). Let h=H(C), let p2(C) count unordered covers of C by two nonempty directed paths, and define

Z2(C)=sum_{S⊆V(C)} H(C[S])H(C[V(C)\S])=2h+2p2(C).

Then

H(X(C))=h+15Z2(C),
H(Y(C))=3h+15Z2(C),
H(Y(C))−H(X(C))=2h.

This is an all-C formula. It is not inferred from small C. Collapse the runs through C in each spanning path to indistinguishable slots. The exterior has seven named vertices. Its exact slot polynomials, indexed by the number of nonempty C-runs, are

X: 15+31t+15t²,
Y: 17+33t+15t².

A k-slot pattern contributes k!p_k(C), since unordered cover paths are assigned bijectively to ordered slots. For C nonempty the constant coefficient contributes zero. Thus X=31h+30p2 and Y=33h+30p2. The displayed finite polynomials are independently checkable by enumerating only7! exterior orders and their eight gaps, with no bound on C. A legal gap may hold C iff its predecessor is absent or one of the first three outer vertices and its successor is absent or a last-three vertex or the apex. Every absent direct exterior arc requires such a gap nonempty. This certifies the two finite polynomials and the all-order substitution theorem.

If C is transitive of order l, Z2(C)=2^l, recovering the dyadic pair. For a nontransitive C, h≥3, so the difference is at least6. Exactly+2 occurs iff C is transitive. Replacing a binary zero clock by a cyclic core therefore cannot preserve the unit repair in this scaffold; its additional odd-cycle mass pays the exact multiplier h.

The two coordinates (h,Z2) are physical and distinct. Z2 is multiplicative under ordinal sum, by partitioning each component independently. Adjoining a singleton in an ordinal gap fixes h and doubles Z2. Hence the scaffold is repeatedly rearmed at every size, but its absolute target changes and its repair tariff remains2h. This is an explicit return/scale separation, and a precise requirement for a richer scaffold that cancels the h-dependence rather than pretending it is absent.
