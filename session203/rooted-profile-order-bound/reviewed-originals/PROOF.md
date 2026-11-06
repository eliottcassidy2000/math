# A logarithmic order bound at bounded rooted start count

**PROOF CANDIDATE, pending independent audit.** This is a finite-search completeness bound for actual rooted profiles. It proves neither that a requested profile exists nor F(8,137), and it gives no eventual odd-spectrum threshold.

## 1. Statement and boundary conventions

For a nonempty rooted tournament (U,v), write h=H(U), s=start_v(U), and let L count ordered spanning two-path covers (P,Q), where P is nonempty and ends at v and Q may be empty. Put

    B(U,v)=h+L-s.

Empty paths have one realization on the empty vertex set, so H(empty)=1. These are the physical B and s of the Session203 staggered-profile criterion.

**Theorem.** If s>0, then there is an induced rooted subtournament (U',v) with the SAME pair (B,s) and

    |U'| <= s+1+floor(log_2 B).                    (1)

U' is obtained simply by deleting every singleton strong component occurring after the root's strong component. No internal arc is changed. In particular, if U is strong then no deletion is needed, and the stronger underlying inequality is

    B >= L >= 2^(|U|-1-s).                       (2)

Negative exponents in (2) have their ordinary rational meaning. Equivalently, L*2^s>=2^(|U|-1).

The singleton has (h,s,L,B)=(1,1,1,1), satisfies (1)-(2), and is handled separately below. A two-vertex tournament is not strong. If its root has positive start count, the root is its source and deleting the other singleton component gives the rooted singleton with the same (B,s)=(1,1).

## 2. Elementary path facts

Every tournament has a directed Hamilton path by successive vertex insertion. Every vertex of a strong tournament starts a Hamilton path and ends a Hamilton path. One direct rooted proof starts with the chosen vertex and maintains a Hamilton path on a visited set K. Unless K is all vertices, strong connectivity supplies an arc x->w from K to its complement. The new vertex w can be inserted without changing the first vertex: append it if the old last vertex beats w; otherwise, since some visited vertex beats w, the indicators [p_i->w] contain a 1 followed eventually by 0, giving an insertion gap. Repeating produces a Hamilton path beginning at the chosen vertex. Converse gives the ending version. In particular, a strong tournament D has H(D)>=|D|. These statements include a singleton.

The strong components of a tournament are totally ordered by dominance: all arcs from an earlier component point to every later component. Hamilton paths visit these components in that order. Thus if the components of U-v are D_1,...,D_k in dominance order, then

    s_v(U) = q * sum_(a in D_1, v->a) start_a(D_1),
    q = product_(i=2)^k H(D_i).                  (3)

This follows by choosing the first vertex after v in D_1 and concatenating the component paths. Empty products are one.

## 3. Strong sources: two physical path-cover injections

Let U be strong of order n>=3. Since h>=s,

    B=h+L-s>=L.                                  (4)

### Case A: U-v is strong

Each summand in (3) is at least one, so the number of out-neighbors of v is at most s. Hence v has at least n-1-s in-neighbors.

For every subset S of its in-neighbor set, choose a Hamilton path of U[S] and append v. Take a Hamilton path on the complementary vertices as Q. All edges into the appended v are legal. When S is empty the designated path is simply (v). These are valid ordered covers counted by L, and different S give different designated vertex sets. Therefore

    L >= 2^(indegree_U(v)) >= 2^(n-1-s).          (5)

### Case B: U-v has at least two strong components

Let D be its terminal component, let m=|D|, and put E=V(U) minus (D union {v}). There is a vertex w in D with w->v. Otherwise D has no outgoing arc to the rest of U: it loses to every earlier component and, under that supposition, to v. This contradicts strong connectivity of U.

Fix a Hamilton path of D ending at w, which exists by Section2. For every subset S of E, choose a Hamilton path of U[S], follow it by that fixed D-path, and then by v. This is the designated path P. Take a Hamilton path on E minus S as Q. Every joining arc from S to D is forward because all earlier components dominate D, and w->v is fixed. Empty S and empty Q are both allowed. The two paths partition every vertex of U. Different S give different first-path vertex sets, so this is an injection of all 2^|E| subsets into the covers counted by L. Hence

    L >= 2^(n-1-m).                              (6)

The factor q in (3) contains H(D). Since s>0, the sum in (3) is a positive integer. Consequently

    s >= q >= H(D) >= m.                         (7)

Combining (6)-(7) proves L>=2^(n-1-s). This completes (2) for nontrivial strong sources and gives n<=s+1+floor(log_2 B). The singleton satisfies the same inequality directly.

## 4. Weak sources: preserving both coordinates

Now suppose s>0 but U need not be strong. Its root lies in the initial strong component C_0, since a Hamilton path cannot return to an earlier component. Write the later components as C_1,...,C_k, and put

    q = product_(i=1)^k H(C_i).

Let h_0,s_0,L_0,B_0 denote the rooted quantities inside C_0. Component concatenation gives h=qh_0 and s=qs_0. It also gives

    L=qL_0,     B=qB_0.                          (8)

For the L identity, a path ending at v cannot contain any vertex of a later component, so its vertex set is contained in C_0. The other path contains the remaining C_0 vertices, possibly none, followed by all later components in order. For each internal cover, that tail has exactly q choices. This proves both directions of the cover bijection, including empty ordinary paths inside C_0.

Every singleton later component has H=1. Delete all of them. Formula(8) shows that h,s,L, and therefore BOTH B and s are unchanged. In the resulting U', retain the notation q for the product over its remaining later components.

If there is no later component, U'=C_0 and Section3 applies. Otherwise q>=3. Each remaining later component has order at most its Hamilton count, and each such count is at least three. For integers a_i>=3, sum a_i<=product a_i; the one-factor case is equality and the inductive step uses ab>=a+b for a,b>=3. Thus

    sum_(i>=1) |C_i| <= q.                       (9)

If C_0 is a singleton, s_0=B_0=1, so s=q and |U'|<=1+q=s+1, proving (1). Otherwise Section3 and (9) give

    |U'| <= s_0+1+floor(log_2 B_0)+q.

Since s=qs_0 and s_0,q>=1,

    s_0+q <= s+1

by (s_0-1)(q-1)>=0. Also q>=3 and B=qB_0 imply

    floor(log_2 B_0) <= floor(log_2 B)-1.

Substitution proves (1). No assumption about the internal arcs, order, or path count of the original singleton tails was needed.

## 5. Consequence for the staggered-profile test

For a specified odd base X>=137, a proposed F(8,137) certificate has

    B(U,v)=X-2j,    s_v(U)=r+3j,
    r in {1,2},    0<=j<=8.

If that profile exists at any order, the theorem supplies an actual certificate of order at most

    r+3j+1+floor(log_2(X-2j))
      <=27+floor(log_2 X).                       (10)

Thus each specified base has a finite exhaustive tournament-search universe with a proved order cap. This is a completeness guarantee, not a claim that such enumeration is cheap. It is not evidence that every base has a certificate: the all-base occupancy statement F(8,137) remains OPEN. Finite success at many bases still does not prove it.

## 6. Controls and attribution

The proof uses only tournament component order, elementary Hamilton-path insertion, and the stated physical cover definition. The component product and rooted-start identities are inherited standard mechanisms; no literature-priority claim is made. The contribution is their combination into the B,s-preserving order bound and its application to the precise staggered-profile certificate.

`verify.py` independently counts Hamilton paths and the designated two-path covers, finds strong components, and performs the literal deletion. It checks every labeled root-zero tournament through order5, the strong-source subset-cover injections, the singleton and rooted-edge boundaries, and extra weak-tail examples with multiple nontrivial components and removable singleton components. Normal and optimized receipts must agree. The finite controls do not replace the all-order proof above.
