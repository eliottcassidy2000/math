# Recovering weak start-one sources and their erased deletion response

Session 197. Exact all-order weak-source membership, with a finite set of possible deletion counts. This is a separate extension of the reviewed strong-profile inverse, not a change to its scope.

## 1. Definitions and inherited strong oracle

For a nonempty rooted tournament (T,r), write D=H(T), E for Hamilton paths ending at r, U for ordered covers (P,Q) with P possibly empty and Q starting at r, and d=H(T-r), with H(empty)=1. We require that exactly one Hamilton path starts at r. The frozen `BACKWARD_PROFILE_INVERSE.md` recognizes every strong start-one triple (D,E,U), plus the singleton and a two-vertex transient, and returns a physical gadget recipe. For a strong member with D>=3 it also proves

    D <= U <= 2(D-E) < 2D,  E>=1,                    (1)

and recovers the deletion count as

    d=U-D+E if 3E>D, or if 3E=D and U>D;
    d=E otherwise.                                  (2)

The singleton has (D,E,U,d)=(1,1,1,1). The present result explains exactly what happens when singleton SCC padding makes E vanish.

## 2. Canonical weak structure

**Theorem 1.** Every non-strong rooted tournament with start count one has the unique form

    C -> Tr_k,  k>=1,                                (3)

where C is the initial strong component, r belongs to C, and (C,r) has start count one. The component C may be a singleton. Its responses satisfy

    D=H(C),  E=0,  U=2^k U_C,  d=H(C-r).             (4)

Conversely every graph in (3) has start count one and responses (4).

**Proof.** Tournament SCCs are totally ordered, and a Hamilton path traverses each SCC in that order. Its first vertex must lie in the initial component. The number of paths starting at r is the start count inside that component times the ordinary H of every later component. All these factors are positive integers. Every nontrivial strong tournament has H>=3: a Hamilton cycle gives at least its order many paths, and its order is at least3. Thus the product can be one only when every later component is a singleton and the initial rooted count is one. This proves (3), its converse, and D and E in (4).

Appending a sink preserves H and d. For every old ordered cover, the new sink can be placed at the end of either of its two paths; if the first path was empty, putting the sink there gives the singleton first path. These are disjoint and exhaustive possibilities, so U doubles. Repeating k times proves the remaining formula. The SCC decomposition makes C and k intrinsic. QED.

The use of a Hamilton cycle in a strong tournament is the standard Camion input already explicit in the inherited low-cycle proofs. No new extremal estimate is needed here.

## 3. Exact numerical oracle, including the erased endpoint

Suppose D>=1 and E=0. A realization exists exactly as follows.

- If D=1, U must be a power2, U=2^k with k>=1. The initial component is a singleton and the only possible deletion count is d=1.
- If D>=3, there must be an integer k>=1 such that U is divisible by2^k and

      D <= U_C=U/2^k < 2D.                           (5)

  This k is unique. For every integer E_C=1,...,D, apply the strong-profile inverse to (D,E_C,U_C). Keep precisely the accepted triples, and recover d from (2). The resulting set of d values is exactly the set possible for the weak triple (D,0,U).
- Other inputs, including even positive D, are rejected by these tests and the strong oracle.

For E>0, no nonempty singleton tail is possible. Test the input directly with the strong oracle, accepting only a strong graph or the singleton, and use (2) or d=1. In particular the count-one two-vertex transient is not a strong-core candidate in this branch.

**Proof.** The dyadic intervals [2^k D,2^(k+1)D) are disjoint. Inequality(1), and the analogous singleton inequality1<=1<2, prove both existence of the normalization for every actual weak graph and uniqueness of k. Exact divisibility is necessary; merely rounding U/2^k is not allowed. The missing E_C is an ordinary positive endpoint count and is bounded by D. The strong oracle is sound and complete, so testing every E_C produces every possible initial component response and nothing else. Formula(2) then gives its exact deletion count. Conversely append k physical sinks to every accepted recipe. Equation(4) gives precisely the requested weak triple and recorded d. QED.

This is a finite existential recovery of an erased coordinate. It does not assert that d is a function of the weak triple. It returns one physical witness for each distinct possible d; it need not enumerate all isomorphism classes with that data. The strong-oracle portion is polynomial in the numerical D, using the inherited memoized inverse. Normalization takes O(log U) exact halvings, with bit cost polynomial in log U; this is an arithmetic-operation statement, not an assertion of linear time in the input bit length. The entire algorithm is not claimed polynomial in the bit length of D.

## 4. A literal collision and a consequential physical branch

Consider these two strong rooted graphs, with vertices numbered from zero:

- A: rows[2,28,25,17,1], root0, profile(9,1,9), deletion count1.
- B: rows[18,28,25,1,8], root3, profile(9,3,9), deletion count3.

Each root has start count one. Appending a sink gives

- A+: rows[34,60,57,49,33,0], root0.
- B+: rows[50,60,57,33,40,0], root3.

Both have the same triple (9,0,18), but their deletion counts remain1 and3. The oracle returns exactly {1,3}; here the erased strong endpoint was respectively1 or3. Therefore extending the strong readout (2) unchanged to weak sources would be false.

For a direct physical illustration, from each weak source add a vertex a whose sole out-neighbor is its root b, then replace a by a transitive pair, retaining the first clone as root. Literal path/cover counting gives respectively

    A+ branch: (H,e,R)=(61,17,89), start count2;
    B+ branch: (H,e,R)=(57,15,87), start count2.

These actual different outputs show that the erased deletion coordinate can matter to a later construction, even though the two weak inputs have identical (D,E,U). This example is directly counted in the verifier; no new universal pair-transfer formula is assumed here.

## 5. Controls, dependencies and scope

`verify_weak_start_one.py` embeds the frozen strong inverse, binds its source/proof hashes, constructs literal recipes, and adds the normalization/oracle above. It checks every weak numerical candidate in a stated bounded rectangle against a complete bounded forward core bank, counts the core deletions directly, and checks all start-one rooted instances among labelled tournaments through order5. Every oracle witness in these physical controls is reconstructed and counted. The two collision graphs and their two adjacent branches are also counted directly. Normal and optimized runs use explicit exception checks.

The all-order argument is the SCC factorization, the dyadic normalization, and the already proved strong inverse. The finite controls do not establish completeness by themselves. This result supplies a weak-source interface for later arithmetic constructions; it does not prove cofinite C1/C2 coverage or collapse weak triples to a universal three-coordinate transfer.

Frozen dependencies:
- `BACKWARD_PROFILE_INVERSE.md`, SHA256 8a07f1bb59f7ee3af9c172a1cfec8551c87125747eeb95dddef7335df33df1fb.
- `verify_backward_profiles.py`, SHA256 3f11b79a5470736e01ce1596a5d7b5fe9b75be9a78489cc285a9cfc2e22cf4b4.
- `ROOTED.md`, SHA256 ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e.
