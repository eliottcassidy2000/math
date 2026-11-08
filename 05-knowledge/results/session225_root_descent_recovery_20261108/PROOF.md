# Affine chamber returns and finite-line boundaries

**Status: PROVED, with freshly reproduced FINITE-EXACT controls and fresh independent audit PASS.**

This is a new reconstruction of unavailable Session225 work, not restoration of its original files. Chan–Pak completeness and every eventual cutoff remain OPEN. No literature-priority claim is made.

## 1. Physical conventions and source normalization

Use the verified physical definitions and operations in [the graft recovery, sections 1 and 5](https://github.com/eliottcassidy2000/math/blob/6872236362fbb7077314e25f3a66a7d0070c5402/05-knowledge/results/session225_graft_recovery_20261108/PROOF.md).

For a nonempty rooted tournament (A,a), D=H(A), E counts Hamilton paths ending at a, d=H(A−a), and U counts ordered two-path spanning covers (P,R) with R nonempty and starting at a and P possibly empty. Set H(empty)=1. Each empty path bin has one realization.

The initial graph is the directed edge 0→1 rooted at 0, with outgoing rows [2,0] and state t=(1,0,2,1). Operations are:

- O adds a global sink and retains the old root.
- Z adds a vertex beaten by every old nonroot and beating the old root; retains the old root.
- J makes the same addition as Z but marks the new vertex.

They act on (D,E,U,d) by

    O:(D,0,2U,d),
    Z:(D−E+U,d,2U−d,d),
    J:(D−E+U,D−E,D+U,D).

The cited graft proof establishes these exact formulas combinatorially. The local verifier additionally constructs literal graphs and checks the whole four-coordinate profile by subset DP.

Words act left to right. Define F=D−E+U. An EXTRA closing J gives Hamilton-path count H=F, by the first coordinate of J. F is not generally the unclosed graph's H. In particular t.OJ has H=5, F=9, while t.OJJ has H=9. We make no cycle-count or minimum-order claim.

## 2. Chamber return formula

An OJ ending sends every state to (D,E,U,d)=(h,e,h,e), with h=D+2U, e=D. Put g=h−e and N=F=h+g. A chamber basis embedding is (h,g)↦(h,h−g,h,h−g).

Let B=b_1...b_n be a nonempty O/Z block ending O, identifying Z with 1 and O with 0. Set

    L=2^n,
    q=Σ_i b_i 2^(n−i),
    a=Σ_i b_i 2^(i−1),
    c=Σ_(j<i) b_i b_j 2^(i−1−j) + Σ_(i=2)^n b_i b_(i−1),
    i=b_1.

The exact block action is

    B(D,E,U,d)=(D−iE+aU−cd,0,LU−qd,d).

Proof: U_j=2U_(j−1)−b_jd, and E_j=b_jd after the first letter. At a Z, D increases by U_(j−1)−E_(j−1). Substitute the U recurrence and sum. The initial E contributes −b_1E; the remaining E terms give the adjacent-digit summand in c. The final E is 0 because B ends O. This is the inherited Session223 digital-block identity, rederived here rather than claimed new.

Substituting the chamber embedding and applying J yields the exact matrix on column (h,g):

    T_B = [[L+a+1−q−i−c, q+i+c],
           [L−q,             q]].

An independent construction applies the literal four-coordinate moves to basis vectors (1,1,1,1) and (0,−1,0,−1). Signed basis vectors here are linear algebra, not physical graph claims.

### One-return Perron boundary

For every such block, 2^n<rho(T_B)<2^(n+1). In particular its Perron factor is not a power of two.

The matrix is nonnegative. Its second column is visibly nonnegative; its first is obtained by applying B J to (1,1,1,1), the rooted singleton state, and reading its physical chamber h,g. Its row sums are L+a+1 and L. Because B ends O, 0<=a<=L/2−1. If q>0, both off-diagonal entries are positive, so strict Perron row-sum bounds give L<rho<L+a+1<2L. If q=0, B=O^n and T_B=[[L+1,0],[L,0]], with rho=L+1 strictly inside (L,2L).

This concerns ONE one-J return block, not arbitrary products or nonbinary constructions.

## 3. Three exact infinite strong affine families

Direct matrix multiplication gives:

    P82=OOZOZOJ,       T82=[[73,12],[54,10]],
    P292=OZOOZZOOJ,    T292=[[217,90],[180,76]],
    P580=OOOZOOZOOJ,   T580=[[545,40],[476,36]].

Their invariant rows are respectively (3,−4), (5,−6), (7,−8). If I is the corresponding invariant, the readout N=h+g obeys respectively

    N'=82N+15I,
    N'=292N+21I,
    N'=580N+63I.

Proof: for each displayed T and invariant row l, verify lT=l and (1,1)T=lambda(1,1)+kappa l. This proves the identities for every chamber state, not only sampled iterates.

The common physical seed t.OJ has (h,e,g)=(5,1,4), N=9. Its invariants are −1, 1, 3. Consequently, sources t.OJ.P^m, each followed by an EXTRA closing J, realize for every integer m>=0:

    N82(m)=(238*82^m+5)/27,
    N292(m)=(880*292^m−7)/97,
    N580(m)=(1800*580^m−63)/193.

The 82 family starts 9, 723, 59271; the 292 family starts 9, 2649; the 580 family starts 9, 5409. These formulas follow by solving the affine scalar recurrences, and integrality also follows directly from physical iteration. Root changes grow with m, so this does not contradict bounded-J selector obstructions.

Every closed member is strong, at all orders. The initial root reaches every vertex of the directed edge. This property survives O because the new vertex is a global sink. It survives Z: the old root reaches some old nonroot (there are always at least two old vertices), which beats the new vertex. It survives J because the new root beats the old root, which reaches every old vertex. Finally, in any closing J, every old nonroot beats the new vertex, and the old root reaches an old nonroot and hence the new vertex. The new vertex reaches the old root and all old vertices. Thus all ordered pairs of vertices are mutually reachable. This is an all-order strongness proof, not an extrapolation from finite cycle counts.

## 4. A completely specified finite line graph

This section is FINITE-EXACT, with an all-order consequence for its fixed menu only. Enumerate every B ending O with 1<=|B|<=11 (2,047 blocks), compute T_B, deduplicate matrices, and retain detT>1. Consider every primitive row l=(a,−b), 1<=a<b<=80. Retain the directed edge target l_t→source l_s exactly when l_t T=l_s. The actual forward state operation goes from the source line to the target line. The edge direction is reversed intentionally for inverse branching counts.

There are 1,981 retained matrices. The strongly connected component containing (3,−4) has 21 vertices and 67 edges. Every target has at most 5 retained incoming forward operations. Every retained determinant within this component is at least 6. CHECK.json includes the complete ordered edge list, matrices and representative words; verify.py regenerates it with no external input. Its 21 vertices, written (a,b), are

    (3,4),(5,6),(5,7),(7,8),(7,9),(7,10),
    (11,12),(11,13),(11,14),(11,15),(13,16),(13,19),
    (15,16),(15,17),(15,19),(19,22),(19,25),
    (23,32),(25,32),(31,38),(31,45).

### Fixed-menu counting theorem

Starting from any fixed finite physical chamber seed bank with positive h,g on these lines with intercept −1, and using only the 67 retained operations, the number of resulting scalar counts at most X is O(X^(log5/log6)). The same conclusion holds with intercept +1, or a finite union of these two cases.

On a source line a_s h−b_s g=−1, write s=H_s(b_s,a_s)+(0,1/b_s), H_s=h/b_s>0. The relation l_t T=l_s implies

    T(b_s,a_s)=det(T)(b_t,a_t).

For example this follows by multiplying the covector equation by the two-dimensional alternating matrix, or by direct 2×2 expansion. Reading the first coordinate gives

    H_t=det(T)H_s + T_12/(b_s b_t)>=6H_s.

The matrices are nonnegative. Also N=(a+b)H+1/b on each line, so H and N are comparable uniformly over this finite line list. A word producing N<=X has length at most log_6(CX), with C depending on the finite seed bank. Reversing from an endpoint leaves at most 5 choices per step; multiply by the fixed endpoint and seed counts and sum the geometric bound. Collisions only decrease the number of outputs. For intercept +1 use s=H(b,a)+(1/a,0), H=g/a>0, and read the second coordinate to obtain the same lower expansion bound. t-generated physical chambers have g>0; positivity is explicitly required of the seed bank.

This does not exclude other lines, longer blocks, larger root-change budgets per return, intercept changes, or other constructions.

## 5. Recurrent finite-intercept restriction

Let T be an integer 2×2 matrix sending the entire affine line l_s x=c_s into l_t x=c_t, where the integer covectors l_s,l_t are primitive and c_s,c_t are nonzero. Then l_t T=k l_s and c_t=k c_s for an integer k. Constancy along the source line proves proportionality; Bezout for the coordinates of primitive l_s proves k is an integer.

In a finite directed line graph, the product of k around any cycle equals 1. Every nonzero integer factor on that cycle must therefore have absolute value 1. Thus every edge in a recurrent strongly connected component has |k|=1. Nonrecurrent edges may have |k|>1. Different nonunit intercepts may supply additional residue classes; this theorem does not say they are useless.

For l_s=(a_s,−b_s), l_t=(a_t,−b_t), associated positive directions v_s=(b_s,a_s), v_t=(b_t,a_t), and invertible T, the determinant identity gives Tv_s=(detT/k)v_t. If T is nonnegative this factor is positive, so k and detT have the same sign. Zero-intercept lines and singular direction collapse are separate cases, not covered by this sign conclusion.

## 6. Finite unit-line obstruction for t-generated returns

Consider finitely many lines a h−b g=±1 with coprime positive integers a<b, optionally together with the boundary h−g=1. Look only at endpoints t.A OJ for arbitrary words A over O,Z,J that lie on one of these lines. Adding any finite set of exceptional scalar counts does not make this family contain every sufficiently large odd count.

This statement is about t-generated endpoints, NOT arbitrary seed closures.

Set x=D−E, W=U−d. The initial t has d>0, W>0, 0<=E<=d, 0<x<=W. The operations become

    O:(x+E,0,2W+d,d),
    Z:(x+W,d,2W,d),
    J:(W+d,x,W+d,x+E).

These preserve the inequalities by direct substitution. In particular D=x+E<=W+d=U. After a final OJ, h=D+2U and g=2U, so g/h>=2/3.

For a line a h−b g=c with c=±1 and a/b<2/3, substitution gives (2b−3a)h<=−3c<=3. It thus has only boundedly many endpoint counts. All remaining unbounded lines have a/b>=2/3 and therefore a>=2. Let M be the least common multiple of 2(a+b) over these remaining lines, or M=2 if none remain. If N≡1 modM, integrality of

    h=(bN+c)/(a+b)

would require b+c≡0 mod(a+b). But 0<b−1<=b+1<a+b for a>=2, b>a. Every such N therefore misses all strict-slope unbounded lines.

The boundary e=h−g=1 cannot be discarded: it arithmetically admits every odd N. For t.A OJ its e equals D(t.A). Each Z or J increases D by U−E>=U−d>0, while O preserves D. Since t has D=1, e=1 forces A=O^k. Then h=1+2^(k+2), hence N=2h−1=2^(k+3)+1. Only O(logX) such outputs occur through X. The progression N≡1 modM has positive density and contains only odd numbers. Removing the bounded strict-slope exceptions, these boundary outputs, and any finite additional exceptional counts leaves infinitely many uncovered odds.

No conclusion about finitely many NONUNIT lines, arbitrary seed closures, non-OJ endpoints or larger predicates is claimed.

## 7. Reproduction and provenance

Run python3 verify.py > replay.json and compare with CHECK.json; repeat with python3 -O and compare with CHECK.optimized.json. No third-party package or input file is required. Every guard uses an explicit exception, not a removable assertion. Controls include all 2,047 block formulas; the entire graph; invariant and readout identities; 11 iterates of each family; 121 literal four-coordinate profiles through word length 4; four literal closed witnesses; and 3,280 t-generated cone/boundary controls through word length 7. Finite checks do not prove the all-order claims, whose proofs are above.

Original Session225 artifacts are unavailable. This packet was reconstructed from retained mathematical descriptions and checked anew. The fresh independent audit approves Sections 1–6, including both narrowly scoped continuation theorems and the all-order strongness bridge. No larger 65-gadget bank, 211-line or 354-line graph, old output hash, or original audit file is asserted recovered here.
