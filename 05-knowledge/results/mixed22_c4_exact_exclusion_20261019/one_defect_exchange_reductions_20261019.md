# One-defect exchange reductions for the seven-thin mixed-(1,2), c=4 frontier

Date: 2026-10-19. Status: proved sufficient balanced-pair reductions. No c=4 exclusion theorem is claimed. This note does not modify the frozen proof or active scan.

## 1. General one-defect exchange lemma

Write D(u), U(u) for strict predecessor and successor sets of u in a finite poset. For incomparable u,v define

    δ(u,v) = |D(v) \ D(u)| + |U(u) \ U(v)|.

**Lemma.** If δ(u,v)=0, then P(u<v)≤1/2. If δ(u,v)≤1, then P(u<v)≤2/3. Consequently δ(u,v)≤1 and δ(v,u)≤1 imply that u,v is balanced.

All probabilities are under the uniform distribution on the actual ambient linear extensions.

**Proof for zero defects.** Swap the labels u and v in an extension with u before v. Moving v earlier is valid because D(v)⊆D(u); moving u later is valid because U(u)⊆U(v). This is an injection into extensions with v before u.

**Proof for one extra successor.** Suppose D(v)⊆D(u) and U(u)\U(v)={t}. In particular u<t, v∥t, and U(t)⊆U(v): any strict successor of t is a successor of u distinct from t and hence a successor of v.

Partition the extensions with u before v into two classes:

1. u before v before t. Swap u,v. The result is valid: the sole potentially obstructing successor t remains after the new u, and all other successors of u follow the old v. This is an injection into extensions with v before u.
2. u before t before v. At the three occupied positions replace the labels (u,t,v) by (v,u,t), leaving every other label in place. The new v is valid in the old u position because D(v)⊆D(u), and its successors stay after it. The new u is valid in the old t position because its predecessors preceded the old u; its successors other than t follow the old v, and t is now in that old v position. The new t moves later, so its predecessors remain before it, and its successors follow the old v because U(t)⊆U(v). This is an injection into extensions with v before u. The three positions recover the input uniquely.

Each class therefore has size at most the number of extensions with v before u. Their sum is at most twice that number, yielding P(u<v)≤2/3. The case D(v)\D(u)={t}, U(u)⊆U(v) is the order-dual statement (applied to the reversed ordered pair). These are all possibilities for δ=1. ∎

No equality of the other downsets or upsets is required. In particular, large differences favorable to v do not invalidate the upper bound.

## 2. Ambient setup and bookkeeping

Use the nine-point core Q={x,y,b,a1,a2}∪Z, |Z|=4, with root relations b<y and x<a1,a2, the root critical-pair inclusions, and the prescribed slave/common incomparabilities. Let A and B be the lower and upper exterior. Then

    D(x)=A,       U(x)=B∪{a1,a2},
    D(y)=A∪{b},   U(y)=B.

A core element's only possible exterior comparison directions are A<core<B. A common z with k core comparable vertices has core incomparability degree 8−k, so the seven-thin hypothesis permits at most k−1 exterior incomparabilities.

Thus k=1 forces every A below z and every B above z. For k=2, at most one exterior relation may be missing. These facts concern full strict downsets/upsets, including transitive relations, rather than a selected Hasse-edge list.

## 3. Unconditional reductions

### 3.1 Sole comparable core vertex b

If the only core comparability involving z is b<z, then

    D(z)=A∪{b}=D(y),  U(z)=B=U(y).

Swapping y,z is a bijection between the two orientations. Hence P(y<z)=1/2.

### 3.2 Sole comparable core vertex an upper slave

If the only core comparability involving z is z<a_i, then

    D(z)=A=D(x),  U(z)=B∪{a_i},
    U(x)=U(z)∪{a_j},  j≠i.

Thus δ(x,z)=1 and δ(z,x)=0. The general lemma gives

    1/2 ≤ P(x<z) ≤ 2/3.

### 3.3 Sole comparable core vertex a common predecessor

If the only core comparability involving z is w<z for w∈Z, then

    D(z)=A∪{w},  U(z)=B,
    D(y)=A∪{b},  U(y)=B.

Both directed defects equal one. Consequently

    1/3 ≤ P(y<z) ≤ 2/3.

This is a third unconditional saturated-vertex reduction. It does not cover a sole common successor z<w.

### 3.4 Exactly two comparable core vertices, both upper slaves

If z's only core comparabilities are z<a1 and z<a2, its core predecessor/successor sets agree with x's. Every exterior element is comparable with x, whereas at most one exterior element is incomparable with z. Thus their full predecessor/successor sets differ by at most one membership in total, and both directed defects are at most one. Therefore (x,z) is balanced for every actual exterior.

This conclusion includes exterior incomparability on either side and does not require testing which exterior vertices exist.

## 4. Further slack-one reductions with a stated exterior condition

These conclusions are genuine ambient-poset reductions. The exterior condition must be proved in a particular application; it cannot be dropped to remove a core from a core-only census.

### 4.1 Two predecessors, one of them b

Suppose z's only core comparabilities are b<z and w<z, with w∈Z. If every B element is above z, then

    D(z)=(A\E)∪{b,w}, U(z)=B, |E|≤1.

Relative to y, δ(y,z)=1 and δ(z,y)=|E|≤1, so (y,z) is balanced. This allows the unique exterior incomparability to be in A, or allows none. It does not establish the conclusion when an element of B is incomparable with z.

### 4.2 One predecessor b and one successor

Suppose z's only core comparabilities are b<z<w, where w is a common vertex or an upper slave. If every A element is below z, then

    D(z)=A∪{b}, U(z)=(B\E)∪{w}, |E|≤1.

Relative to y, δ(z,y)=1 and δ(y,z)=|E|≤1, so (y,z) is balanced. This allows the unique exterior incomparability to be in B, or allows none. It does not establish the conclusion when an element of A is incomparable with z.

## 5. Scope

These proofs need neither the maximum-span majority-order hypothesis nor enumeration, and they hold for arbitrary finite actual exterior size and relations. They therefore give unconditional exclusions of the stated morphologies from any balanced-pair-free seven-thin poset. They do not classify every two-core-neighbor common vertex, and they do not prove the full mixed-(1,2), c=4 branch.

## 6. Uniform core-only sufficient criterion

For either root r∈{x,y} and a core vertex v incomparable with r, let δ_Q use strict predecessor/successor sets restricted to Q. Write q_A and q_B for the number of elements of A and B, respectively, incomparable with v. Since the root has all exterior relations, and no core element has the reverse exterior directions,

    δ_P(r,v)=δ_Q(r,v)+q_B,
    δ_P(v,r)=δ_Q(v,r)+q_A,
    q_A+q_B≤s(v):=7−d_Q(v).

In particular the entirely core-checkable condition

    max{δ_Q(r,v),δ_Q(v,r)} + s(v) ≤ 1

certifies that (r,v) is balanced for every permitted actual exterior. This is a sufficient condition, not a necessary one; it can discard information about which side can actually contain an exterior incomparability. It applies equally to common vertices and slaves. The conditional reductions in Section 4 use the sharper separate q_A,q_B formulas.

## 7. Any-pair core-only sufficient criterion

For arbitrary incomparable core u,v, set ε(r)=0 for each root and ε(t)=7−d_Q(t) for a nonroot core vertex t. A root has no exterior incomparability by the definition of Q. Write A_t={a∈A:a<t} and B_t={b∈B:t<b}. The exterior contribution to δ_P(u,v) is

    |A_v\A_u|+|B_u\B_v|
      ≤ |A\A_u|+|B\B_v| ≤ ε(u)+ε(v).

The same bound holds in the reverse direction. Consequently

    max{δ_Q(u,v),δ_Q(v,u)}+ε(u)+ε(v)≤1

certifies that (u,v) is balanced for every actual permitted exterior. This sufficient criterion can detect nonroot witnesses missed by the root-only criterion in Section 6.
