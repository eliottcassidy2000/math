# Beyond the mirror: the exact low spectrum of cyclic three-module substitutions

Session 175. **PROVED**, with finite independent controls. No global spectrum-completeness or historical-priority claim.

## Result

Let A,B,C be arbitrary nonempty tournaments, and let C3[A,B,C] orient every intermodule arc A→B→C→A. Then

    { H(C3[A,B,C]) : H(C3[A,B,C]) ≤ 45 }
      = {3,5,9,11,15,17,23,27,33,37,45}.

In particular 7 and 21 never arise in this constructor, even after entirely abandoning the mirror relation B=A*. This closes the mirror-breaking option proposed at the end of Session 174's synthesis: escaping that *architectural* obstruction requires leaving the uniform cyclic three-module architecture, not just changing the two seeds independently.

The global exclusion of 7 and 21 is inherited from `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. Completeness of the unrestricted positive odd spectrum remains OPEN. The present result gives an explicit construction-level formula and a sharper low-value classification, without a census of unrestricted tournaments.

## 1. Positive universal cyclic skeleton kernel

Let q_X(k) count ordered covers of X by k nonempty directed paths, so q_X(k)=k!p_X(k); set it to zero for k>|X|. Let a_k=q_A(k), b_k=q_B(k), c_k=q_C(k). Then

    H(C3[A,B,C]) = Σ_{r≥1} [
       3a_r b_r c_r
       + a_{r+1}b_r c_r + a_r b_{r+1}c_r + a_r b_r c_{r+1}
       + a_{r+1}b_{r+1}c_r + a_{r+1}b_r c_{r+1}
       + a_r b_{r+1}c_{r+1} ].                       (1)

**Proof.** Collapse maximal runs within modules in a Hamilton path. Consecutive run labels must follow A,B,C,A,B,C,…, since these are the only forward intermodule arcs. Every module must occur. A skeleton with 3r runs has three possible starting labels, all with counts (r,r,r). Skeletons with 3r+1 or 3r+2 runs each have three starts, and their counts are respectively a permutation of (r+1,r,r) or (r+1,r+1,r). Each skeleton is filled uniquely by an ordered path cover in each module, ordered by the runs in which its paths appear. Summing gives (1). This is a bijection; no cancellation or assumed realizability of abstract profiles is involved. ∎

The formula specializes for C a singleton, writing a=H(A), b=H(B), u=q_A(2), v=q_B(2), to

    H(C3[A,B,1])=3ab+ub+av+uv.                    (2)

For B=A* this becomes 3a²+2au+u², the Session 174 norm formula. The norm is the diagonal of a more general positive bilinear kernel; the mirror arithmetic restriction disappears off diagonal, but the size obstruction persists.

## 2. Realizable profile lower bounds

For an m-vertex tournament X,

    q_X(k) = Σ_(ordered nonempty k-block partitions) Π_i H(X[S_i])
           ≥ k! S(m,k),                           (3)

where S(m,k) is a Stirling number of the second kind. Each induced subtournament has at least one Hamilton path, by elementary insertion. Equality in (3) holds for transitive X. In particular H(X)≥1 and q_X(2)≥2^m−2.

Define F_X(2)=q_X(2)+2H(X), allowing two labelled path slots with either slot empty. If v is an extra vertex, insertion into whichever slot is selected gives

    F_(X+v)(2) ≥ 2 F_X(2).                         (4)

To see injectivity, choose one deterministic feasible insertion into each selected slot. Deleting v recovers both old slot paths, and its occupied slot identifies which of the two choices was made. Thus nontransitive X, which contains a cyclic triangle, has

    H(X)≥3,  F_X(2)≥12·2^(m−3)=3·2^(m−1).        (5)

The triangle's F(2)=12 follows from H=3 and q(2)=6. Bounds (3)–(5) apply to actual tournament profiles, not arbitrary nonnegative vectors.

## 2a. Exact fixed-size minimization, at every order

**Theorem.** For fixed positive module sizes (a,b,c), the unique minimizing module types are transitive: H(C3[A,B,C]) is minimized if and only if each of A,B,C is transitive. The minimum is given explicitly by (1) with q_X(k)=k!S(|X|,k).

**Proof.** Bound (3) is coefficientwise and (1) has only nonnegative coefficients. If any module is nontransitive, it contains a cyclic triangle, and insertion gives its q_1=H≥3, strictly above the transitive value 1. The term 3a_1b_1c_1 then increases strictly, while no other term decreases. Conversely transitive modules attain every lower bound. ∎

This removes interior orientation from the fixed-size minimization problem completely, rather than checking particular mirrored or transitive seed examples. It is a statement about this constructor, not the unrestricted minimum over strong tournaments.

## 3. Complete case reduction below 46

Cyclic rotation lets us choose which module is singled out. There are three cases.

**All three modules have at least two vertices.** The profile of a transitive two-vertex tournament is q_1=1, q_2=2, others zero. Coefficientwise positivity in (1) and (3) give

    H≥3+3·2+3·4+3·8=45.

If any module has at least three vertices, its q_2≥6 instead of 2. Already the positive singleton-increment term in (1) raises the lower bound by at least 4; hence H≥49>45. Consequently 45 occurs exactly when all three modules have two vertices.

**Exactly one singleton, with other sizes m,n≥2.** Formula (2) implies

    H≥3+(2^m−2)+(2^n−2)+(2^m−2)(2^n−2)
      =(2^m−1)(2^n−1)+2.                         (6)

If min(m,n)≥3 the bound is 51. If min(m,n)=2 and max(m,n)≥4, the bound is 47. Thus only sizes (2,2,1) and (3,2,1), up to cyclic arrangement, remain. The first has H=11. A three-vertex module is transitive (a,u)=(1,6) or cyclic (3,6); the two-vertex module has (b,v)=(1,2). Formula (2) gives respectively H=23 or 33. The formula is symmetric in A,B, so both cyclic arrangements have these values.

**At least two singletons.** Let X be the remaining m-vertex module. Formula (2) gives

    H(C3[X,1,1])=3H(X)+q_X(2)=H(X)+F_X(2).        (7)

The lower bound is 2^m+1, so m≥6 is excluded. Transitive X, for m=1,…,5, yields 3,5,9,17,33. Nontransitive X has m≥3. At m=5, (5) gives H≥3+48=51. At m=3 the cyclic triangle gives H=15. At m=4 choose a cyclic triangle and the remaining vertex. If that vertex beats all or none of the triangle, H(X)=3 and there is one cyclic triangle. If it beats one or two, H(X)=5 and there are two cyclic triangles. These statements follow directly by inserting the vertex into the triangle's three Hamilton paths, or by the eight possible neighborhoods. For c cyclic triangles, the ordered 1+3 and 3+1 partitions contribute 2(4+2c), and the ordered 2+2 partitions contribute 6, so q_X(2)=14+4c. Equation (7) gives 27 or 37.

Combining the three cases proves that the only possible values ≤45 are the stated eleven. Every value has been given an explicit construction in this case analysis, proving equality of sets. ∎

## 4. What this does and does not establish

- New local deliverable: a positive all-profile formula and exact ≤45 spectrum of this entire arbitrary-module constructor, strengthening Session 174's mirrored-family obstruction.
- Inherited tools: tournament insertion, path-cover profiles, and uniform-module skeleton factorization. The cyclic skeleton enumeration here makes that factorization explicit for every number of runs.
- No new proof of global exclusion is claimed, and no assertion that all other odds occur follows.
- Counts 13,19,25,29,31,35,39,41,43, although globally relevant, also lie outside this architecture. For a constructor pursuing unrestricted completeness, this is a concrete diagnostic rather than another symmetry restriction.
- The map sends three full path-cover profiles to one Hamilton count. It preserves the target count by a literal bijection, discards graph isomorphism and most profile information, and needs no self-duality sidecar.

## Reproduction

Run `python3 check.py > certificate.json` and `python3 -O check.py > certificate.optimized.json`, then compare. The dependency-free checker exhausts all labelled A,B of sizes 1–4 against literal tournament subset-DP for (2), checks all labelled module interiors at sizes (2,2,2),(3,2,2),(3,3,2), checks the full kernel on all 1,331 triples of labelled tournaments of sizes 1–3 using independently enumerated ordered-cover profiles, and checks selected wrappers by independent permutation enumeration. The universal conclusion is established by the proof above; the finite tests are controls, not substitutes for it. No software was installed, and no prohibited ledger was read.
