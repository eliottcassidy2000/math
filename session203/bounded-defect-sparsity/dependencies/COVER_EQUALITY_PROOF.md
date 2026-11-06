# Exact structure of a rooted endpoint-cover equality

Status: RESERVED / proof candidate awaiting independent audit.

## Scope

Let U be a nonempty tournament rooted at v. Write h=H(U), and let L count ordered spanning two-path covers (Q,P) with Q nonempty and ending at v; P may be empty. A path is its ordered vertex sequence. Let I be the set of inneighbors of v and A=U-v. Empty Hamilton counts equal one.

**Theorem.** L=h if and only if either I is empty, or I={x} and start_x(A)=0. More quantitatively,

    L-h >= sum_{x in I} start_x(A) + binom(|I|,2).    (1)

When |I|<=1 the first summand is the exact value of L-h.

This sharpens the general necessary condition L=h => start_v(U)=H(U-v). It applies to every tournament, not only strong tournaments. It does not prove staggered profile occupancy or an eventual Hamilton-spectrum threshold.

## 1. A two-terminal cover lemma

We use standard tournament endpoint facts: the strong components are totally ordered, a Hamilton path can start at any vertex of the first strong component, and can end at any vertex of the last one. For a strong component this follows by opening a Hamilton cycle at the requested vertex (the singleton case is immediate); concatenate the component paths. Conversely a Hamilton path cannot leave a later component for an earlier one. These are the usual Redei/Camion facts, not new literature claims.

First suppose S is strong and x,y are distinct vertices. At least one of the following holds:

- S-x has a Hamilton path ending at y;
- S-y has a Hamilton path ending at x.              (2)

If both failed, let K_x be the last strong component of S-x, which excludes y, and K_y the last component of S-y, which excludes x. Because the excluded vertex y is outside K_x, K_x remains a strong sink component of S-{x,y}; similarly K_y does. Uniqueness of the last strong component gives K_x=K_y=K. Since K is a sink in S-x, any edge leaving K in S must go to x. Strong connectivity of S requires at least one such edge K->x. But K=K_y is a sink in S-y and x lies outside it, so every edge between x and K is x->K. This is a contradiction, proving (2).

Now let A be an arbitrary tournament and x,y distinct vertices. There exists a spanning cover by two nonempty directed paths, one starting at one of x,y and the other ending at the other vertex.

To see this, order the strong components D1,...,Dk. If x and y lie in different components, assume x in Di and y in Dj, i<j. Let P start at x and cover Di and all components after Di except Dj. Let Q cover the components before Di and then Dj, ending at y. All cross-component concatenations are legal. Both paths are nonempty and cover A disjointly.

If x,y belong to the same component Di, apply (2) there. For example, if Di-x has a Hamilton path ending y, put that path after all components before Di to form Q. Put x before all components after Di to form P. The reverse alternative is identical with x,y swapped. This proves the two-terminal cover lemma.

## 2. Exact bad-cover decomposition

Cutting a Hamilton path immediately after v gives exactly those L covers with P empty or v->first(P). Thus L-h counts the remaining covers with P nonempty and first(P) in I.

Delete the terminal v from Q. If Q was a singleton, the result is a Hamilton path P of A beginning in I, contributing exactly sum_{x in I} start_x(A). Otherwise the result is a cover (Q',P) of A by two nonempty paths, with last(Q') and first(P) both in I. These vertices are distinct because the paths are disjoint.

For each unordered pair {x,y} subset I, the two-terminal cover lemma supplies at least one such cover whose two specified endpoints are x,y. Covers associated with different endpoint pairs are distinct. Therefore the nonsingleton-Q contribution is at least binom(|I|,2). This proves (1).

If I is empty, both contributions vanish. If I={x}, the two nonempty paths cannot have the same vertex as their two specified endpoints, so the second contribution vanishes and L-h=start_x(A). If |I|>=2, (1) is strictly positive. The theorem follows.

## 3. Complete physical normal form at the sharp Z-defect floor

This section uses the separate frozen `dependencies/SHARP_DEFECT_PROOF.md` only for its exact equality test. There, D=s for a nontrivial positive-start rooted U was shown equivalent to: v has a unique outneighbor a, L_a(U-v)=H(U-v), and end_a(U-v)=0.

**Corollary.** Every such equality source, and only such a source, has the following form. Take a nonempty tournament T, add v,a, orient v->a and T->v, and choose one of two cases:

1. a dominates all of T;
2. choose x in T with start_x(T)=end_x(T)=0; let x->a and let a dominate T-x.

In both cases s=H(T), D=s, and U is strong. The condition on x means precisely that x lies in a strong component of T which is neither its first nor its last component.

**Proof.** The sharp-defect equality test gives a unique root outneighbor a. Put A=U-v and T=A-a. The present theorem says that a has either indegree zero in A, or a unique inneighbor x with start_x(T)=0. In the latter case, deleting the terminal a gives exactly

    end_a(A)=end_x(T).

Thus the remaining zero-end requirement is precisely end_x(T)=0. The earlier equality test also guarantees T is nonempty. This proves necessity.

Conversely, in the first case L_a(A)=H(A) and end_a(A)=0 because a is a source and T is nonempty. In the second, the present theorem and the displayed endpoint identity give the same two conclusions. The unique-root-outneighbor equality test therefore gives D=s. Deleting the initial v and then a shows s=H(T)-start_x(T)=H(T), with no subtracted term in the first case.

Every vertex of T beats v. In the first case a reaches all of T directly. In the second case x is outside the first strong component of T, so a dominates that entire component and reaches every vertex of T through it. Together with v->a this proves U is strong in either case. The endpoint support facts identify the allowed locations of x. QED.

This is a complete graph normal form, stronger than merely matching the numerical boundary. It does not preserve B when deleting v,a and does not imply all-target staggered occupancy. The first case is the previously known simple cyclic construction. The second case is essential: for T transitive of order three, its middle vertex has both endpoint counts zero and supplies a genuine additional family.

Sections 1-2 stand independently of the sharp-defect packet. The computation checks all 202,013 rooted labeled inputs through order 6 and all 502,170 unordered terminal pairs. It also verifies the normal form for every equality input and physically constructs both cases from every labeled T through order 6. These controls are not the all-order proof.
