# The elementary A,K,D grammar is exactly the endpoint-one spectrum

Session 195. **PROVED elementary graph equivalence and scalar-support equality; independent review passed.** This does not prove the unrestricted `{7,21}` spectrum.

## 1. Definitions

A rooted tournament starts at a singleton. The elementary grammar uses:

- A: add x with x->r and u->x for every old u other than the root r, and mark x;
- K: make the same vertex birth and keep marking r;
- D: add a universal sink and keep marking r.

Let G be its scalar H support, allowing arbitrary words and weak outputs. Let E be the scalar H support of all tournaments with some vertex whose Hamilton-path terminal count equals one. A tournament and its converse have the same H. An endpoint-one vertex becomes a start-one vertex in the converse.

The endpoint-recursion packet classifies every strong endpoint-one rooted tournament (T,v). Aside from a singleton, it is obtained from a smaller strong endpoint-one core (C,u) by taking a transitive chain Z=(z_1->...->z_m) with Z->C, adding v with u->v->C-u, and specifying a binary word w where b_i=1 iff v->z_i. If m>0, b_1=1; if m=0, |C|>1. Every endpoint-one tournament has the same H as its marked strong sink component, or as a singleton, because every preceding component is a neutral transitive singleton.

## 2. Exact compilation of an endpoint macro

**Lemma 1.** If C^op is rooted at u, the converse of the endpoint macro with word w=b_1...b_m is constructed by

    A, followed by b_m,b_(m-1),...,b_1,

where a bit 1 is implemented by K and a bit 0 by D. The empty word is implemented by A alone.

**Proof.** The first A creates the new root v, with v->u and C-u->v, exactly the converse core incidences. Keep this new mark while reading the reversed bits. The vertex created for b_i is z_i. It loses to all C vertices. If b_i=1 it beats only v among the old vertices, by K; if b_i=0 it also loses to v, by D. Hence z_i->v iff b_i=1, as required by taking the converse. Every later-created z_i loses to every earlier-created z_j with j>i, giving z_j->z_i. This is precisely the converse order on the chain. No old edge is altered. QED.

Thus the complete endpoint-one inverse decomposition compiles recursively into an actual elementary grammar word. This is a graph-level statement, not just an equality of three numerical coordinates. It proves E subset G, including weak endpoint-one witnesses after deleting their neutral components.

## 3. Start-one invariance and the initial-K exception

Let s_r(T) count Hamilton paths starting at the root.

**Lemma 2.** A and D preserve s_r exactly. K also preserves s_r whenever |T|>=2. K on the initial singleton instead gives s_r=0.

**Proof.** For A, every Hamilton path beginning at the new x must next visit r, since x has no other old out-neighbor. Removing x gives exactly an old root-starting Hamilton path. For D, every Hamilton path is an old Hamilton path followed by the new sink.

For K, a Hamilton path starting at old r cannot contain x internally: x's only outgoing old edge goes back to the already used r. Thus x must be its final vertex. When |T|>=2, every old Hamilton path starting at r ends at another old vertex, which beats x, so appending x is always legal and is a bijection. When T is a singleton, this append would require the absent r->x edge, and the rooted start count becomes zero. QED.

**Lemma 3.** Every initial consecutive K operation can be deleted from a word without changing its scalar H output. The corresponding created vertices form a permanent transitive source prefix, never containing the root.

**Proof.** A word beginning with K^t creates vertices x_1,...,x_t with

    x_1 -> x_2 -> ... -> x_t -> r.

The root remains the initial r. Each x_i beats every vertex present after this prefix except earlier x_j, which beat it. In every later A or K birth, the new vertex's only old out-neighbor is the current root, which is never one of the x_i; therefore every x_i beats that new vertex. The same holds for D. Root changes only mark a fresh vertex, so the prefix is never marked later. Hence all x_i beat the entire subsequent graph, and remain a transitive source prefix. Deleting them leaves exactly the graph constructed by the word with K^t erased, with labels appropriately shifted and the same root. Their ordered factor in H equals one. QED.

After erasing that prefix, the word is empty or begins with A or D. It begins from start count one, and the first step creates at least two vertices with start count one. Lemma 2 then preserves that count forever. Therefore every scalar grammar value has a start-one rooted representative; taking its converse gives an endpoint-one witness. This proves G subset E.

## 4. Exact conclusion and the limitations it transports

**Theorem.** G=E: the arbitrary-depth elementary A,K,D grammar has exactly the scalar support of all endpoint-one tournaments, at every order and count. Every strong endpoint-one witness has a literal converse construction in the elementary grammar.

No equality of marked state sets is claimed. The four-coordinate elementary grammar permits weak graphs, and its d coordinate has no counterpart in the three-coordinate endpoint macro interface. Numerical state counts can therefore differ even though scalar supports agree.

The dimensional difference has a concrete explanation. Every compiled endpoint macro begins with A, which resets the deletion coordinate to the old h. The subsequent K and D operations retain that deletion coordinate. Hence a completed macro depends only on the old (h,e,R), corresponding under converse to the endpoint state (h,s,L), even though unrestricted primitive continuations require four linear coordinates. There is no conflict with the four-dimensional minimality theorem for the primitive grammar.

In particular the common absence of H=13 is structural, not a finite-census coincidence. All 29 endpoint-class omissions at most 10,000 are excluded from A,K,D at arbitrary depth. Conversely every positive witness produced by either mechanism transfers to the other scalar class.

This theorem does not make the elementary grammar complete for unrestricted tournament H. It proves why the independently computed omission sets agree, and why leaving this exact class requires a genuinely larger gadget family or a different operation.

## 5. Verification

`verify.py` checks the macro compiler at graph level, including reversed bit order and the root-preserving permutation of labels. It checks all elementary words through length seven, verifies the removal of their complete initial K prefix, checks the permanent source property, compares original/reduced scalar counts by exact subset DP, and verifies start-one for every reduced output. It uses no source-packet counting code, and every check is an explicit exception rather than an assertion. Both infinite directions are proved above; finite controls are not extrapolated.
