# The binary R register cannot be closed by one fixed one-port context

Status: PROVED, independently audited; see audit/REVIEW.md. This is a limited obstruction, not a disproof of F(8,137).

## The tempting arithmetic register

Start at the root-first edge and apply arbitrary O/Z words. In the inherited nine-coordinate state p=(h,s,e,L,R,d,C,b,T), s=d=1 throughout, and

    O: R -> 2R,       Z: R -> 2R-1.

At length m, R runs exactly once through every integer in [2^m+1,2^(m+1)]. This follows by induction: the two maps give the even and odd members of the next interval. Also b=2^(m+1), and e is the parity of R, including the empty word. Thus R is already an exact binary register for all integers at least two.

If a fixed physical closing context converted this register to B=2R+an odd constant while retaining bounded positive start, it would give an immediate odd tail. A broader family of potentially useful readouts is

    alpha R + beta b + gamma e + kappa.                 (1)

The following theorem excludes every such exact identity for one fixed one-port context on the whole O/Z orbit.

## Context definition and transfer inheritance

A fixed context is a finite nonempty tournament G with anchor v and a marked vertex r in G (possibly r=v). Glue v to the old root of U, and make every old vertex in U-v dominate every vertex in G-v. Mark r in the resulting tournament, denoted G[U]. The context is fixed independently of the word used to construct U.

For r distinct from v, the published all-context transfer theorem in session203/start-two-gates/PROOF.md gives a linear action at every cover level k. In particular:

- H(G[U]) is a fixed linear combination of (h,e,R,d);
- L_r(G[U]) is a fixed linear combination of (C,L,T,b).

For r=v, the same needed L conclusion has a direct proof. In a designated root-ending cover, if the designated context skeleton consists only of v, its optional bulk prefix contributes E. The other bin has an optional arbitrary bulk prefix f, giving [A]Ef=L. If a nonroot context vertex precedes v in the designated skeleton, its initial bulk prefix is arbitrary and gives f instead of E, so the total is [A]f^2=b. All remaining context incidences are fixed. Precisely, L_v(G[U])=H(G-v)L+(L_v(G)-H(G-v))b: the first coefficient counts context skeletons with designated singleton v, and the second counts all other designated root-ending skeletons. Both coefficients are nonnegative. This includes singleton G, where the map is simply L_out=L. The H transfer is independent of the marked r and is proved directly below.

For every nontrivial U, start_r(G[U])=c_G s(U), for a fixed nonnegative integer c_G. To see this, a Hamilton path starting at context vertex r cannot enter the old bulk from any context vertex other than v. Its old bulk vertices therefore form one run immediately after v. Choices of that run are exactly the old root-start Hamilton paths, with v removed. Its preceding and following context pieces give a fixed number c_G of choices. Nontriviality excludes the empty-bulk exceptional case. No assumption that c_G is positive is needed below.

## The unavoidable h coefficient

**Lemma.** In the linear expression for H(G[U]) in coordinates (h,e,R,d), the coefficient of h is exactly H(G-v), with H(empty)=1. It is therefore strictly positive for every tournament context G.

Write A=U-v. In the squarefree vertex algebra of A, let f enumerate an optional arbitrary path, E an optional path ending at an inneighbor of v, and S an optional path starting at an outneighbor of v. Then

    h=[A]ES, e=[A]E, R=[A]fS, d=[A]f.

An output Hamilton path can have at most two old-bulk runs: an initial run before its first context vertex, and a run immediately after v. No other entry from a context vertex into bulk is possible.

First consider paths in which v is the first context vertex. The remaining context vertices form a Hamilton word Q of G-v. The optional bulk prefix before v contributes E. If Q is nonempty, put z=1 when v beats the first vertex of Q, and z=0 otherwise. A nonempty bulk run after v always bridges to Q, while an empty run is legal precisely when z=1. Thus the contribution for Q is

    [A] E((S-1)+z) = h+(z-1)e.

If Q is empty, the contribution is h. Summing over all Q yields h coefficient H(G-v).

When v is not the first context vertex, the initial bulk run is arbitrary, contributing f. The optional run after v contributes either S or S-1, according to whether the next context edge permits an empty run; if v is last, it contributes S. These terms are combinations of R and d only, and cannot alter the h coefficient. This proves the lemma and provides a direct path-counting interpretation of that coefficient.

Consequently the physical scalar

    B(G[U],r)=H(G[U])+L_r(G[U])-start_r(G[U])

is a linear functional of the old nine-state whose h coefficient is H(G-v)>0. The other two summands use only (C,L,T,b) and s, respectively.

## The O/Z orbit has no further linear relation

Every state in the rooted-edge O/Z orbit satisfies s=d. Conversely, the following eight actual word states span that eight-dimensional hyperplane. In coordinate order (h,s,e,L,R,d,C,b,T), use the words

    empty, O, Z, OO, OZ, ZO, ZZ, OOZ.

After omitting coordinate d, their row determinant is -1152, as recorded with exact states in SPAN_WITNESS.json. Their span therefore has dimension eight. Since all have s=d, their span is exactly that hyperplane. This finite determinant proves an exact span statement; no finite-to-infinite inference about occupancy is involved.

**Theorem.** There is no fixed context G,r and no constants alpha,beta,gamma,kappa for which (1) equals B(G[U],r) for every U in the rooted-edge O/Z orbit.

If such an equality held, replace its constant kappa by kappa*s, since s=1 on the orbit. The difference of the two linear functionals would vanish on all eight displayed states, so it would be a multiple of s-d. In particular its h coefficient would be zero. But the physical context scalar has h coefficient H(G-v)>0, while (1) has h coefficient zero. This contradiction proves the theorem.

## Scope and surviving route

This rules out one fixed one-port closing context that reads the exact binary R register, even allowing the length register b and the last-bit register e in the desired scalar formula. It does not exclude different contexts on different word classes, simultaneous multi-port responses, target-dependent rebuilding, or an arithmetic invariant using h and L rather than eliminating them. It does not exclude general finite-register constructions and does not disprove F(8,137).

The exploratory enumeration of contexts of orders at most five suggested the lemma, but no bounded context enumeration is used in its proof. Only the universal physical transfer theorem, direct path decomposition, and explicit eight-state determinant are dependencies.
