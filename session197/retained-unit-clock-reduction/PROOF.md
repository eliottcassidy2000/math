# Retained-root unit binary clocks are hidden outdegree-one clone clocks

Session 197. An exact physical inverse reduction, with a more general sink-module transfer. This does not classify all unit binary amplitudes or prove scalar coverage.

## 1. The retained-root operation and its initial boundary

For a nonempty rooted tournament (T,r), write

- h=H(T), e=the number of Hamilton paths ending at r;
- d=H(T-r), with H(empty)=1;
- R=the number of ordered path covers (P,Q) of T, where P may be empty and Q is nonempty and starts at r.

The operation Z appends a vertex beating r and losing to every other old vertex, while retaining the old root r. Its physical transfer is

    Z(h,e,R,d)=(h-e+R,d,2R-d,d).                    (1)

Thus for t>=1 successive births have

    h_t=(R-d)2^t+h-e-R+2d,
    e_t=d,  R_t=(R-d)2^t+d,  d_t=d.                (2)

In particular the eventual affine offset is R-h+e-2d. The first offset h_1-2h is R-h-e. Consequently the exact recurrence

    h_(t+1)=2h_t+beta for every t>=0                (3)

holds if and only if e=d; in that case beta=R-h-d. This criterion includes the first step. An arbitrary first transient is not an exact returning clock at the original count.

The singleton has e=d and beta=-1 and gives the constant count 1. The structural reduction below assumes |T|>=2; every nonconstant unit-offset instance automatically avoids this singleton exception.

## 2. An arbitrary physical sink-module transfer

Let (U,r) be any nonempty rooted tournament with profile (h,e,R,d). Let S be any nonempty tournament on disjoint vertices. Form T by keeping U and S intact and setting

    U-r -> S -> r.                                 (4)

This specifies every new arc, without changing the root. Put G=H(S), and let

    B_S(2)=sum_(A subset V(S)) H(S[A])H(S[V(S)-A])

count ordered two-path covers with either slot allowed empty. Then the exact profile is

    H(T)=(R-d)B_S(2)+(h-e-R+2d)G,
    e_r(T)=d_r(T)=dG,
    R_r(T)=(R-d)B_S(2)+dG.                         (5)

In particular

    R_r(T)-H(T)-d_r(T)=(R-h+e-2d)H(S).             (6)

### Proof of the module identities

Every Hamilton word has at most two S fragments: leaving an S fragment for the exterior requires its unique possible exterior destination r, which can occur only once. For an ordered cover reserved at r, each S fragment must terminate its path slot, since r is already the first vertex of the reserved slot and cannot be revisited. There are therefore at most two S fragments in this setting also.

For either statistic, a one-fragment skeleton lifts in H(S) ways. A two-fragment skeleton lifts in B_S(2)-2H(S) ways, using two ordered nonempty paths covering S. All exterior incidences are uniform on S, so the skeleton count is independent of the internal tournament S. Hence each statistic is a linear combination of these two quantities.

Determine the two coefficients using S=Tr_1 and S=Tr_2. These are respectively one and two actual Z births, because the later newborn loses to the earlier one. Their H values from (1) are h-e+R and h-e+3R-2d. Their reserved-cover values are 2R-d and 4R-3d. Also H(Tr_1)=H(Tr_2)=1, while B_(Tr_1)(2)=2 and B_(Tr_2)(2)=4. Solving gives the first and third rows of (5). This is an exact fragment decomposition for arbitrary S, not interpolation asserting that arbitrary tournaments are transitive.

After removing r, the graph is (U-r)->S, so d_r(T)=dG. Every Hamilton path of this deletion ends in S, which points to r; appending r gives e_r(T)=dG. This proves all rows and (6).

## 3. The terminal component is forced by a unit offset

Suppose now that a physical rooted T with at least two vertices satisfies e_r=d_r. Let S be the terminal strongly connected component of D=T-r in its transitive SCC order.

Every Hamilton path of D ends in S. Moreover every vertex of S occurs as an endpoint of such a path. For completeness, in a strong tournament choose a longest path ending at a specified vertex. If a vertex outside that path points to any path vertex, it can be prepended or inserted immediately before the first path vertex it beats. If no outside vertex can be inserted, every path vertex points to every outside vertex, contradicting strongness. Thus the path is Hamilton. Concatenating arbitrary Hamilton paths of earlier components proves the endpoint assertion for D.

The equality e_r=d_r says that every Hamilton path of D can be extended by r. Therefore every vertex of S points to r. All earlier components point to S, so S is exactly the module in (4), with U=T-S. Formula (6) applies and proves

    H(S) divides beta=R_r(T)-H(T)-d_r(T).            (7)

If beta is +1 or -1, then H(S)=1. A nontrivial strong tournament has a Hamilton path ending at each of its at least three vertices, so it has at least three Hamilton paths. Therefore S is a singleton v. Its orientation is

    T-{r,v} -> v -> r,

so v has exactly one out-neighbor, namely r. This v is unique because it is the terminal component of T-r.

**Unit-clock reduction.** Every exact, nonconstant retained-root clock h_(t+1)=2h_t+beta, beta in {+1,-1}, already has a unique old vertex v with sole out-neighbor r. Repeating Z at r literally replaces this v by a growing transitive module. No rebuilding or change of scalar representative is required.

The complete literal clone identity, indexed by positive clone multiplicity m, is

    H(T[v -> Tr_m])=((H(T)+beta)/2) 2^m-beta.        (8)

The original graph is m=1, and each further Z birth appends the next clone after all previous clones. The same old root r remains exposed. The amplitude is positive in the nonconstant case.

Conversely, if v has sole out-neighbor r and its literal transitive-clone sequence has the form A2^m-beta with beta in {+1,-1}, then v is the universal sink of T-r, so e_r=d_r. Adding clones of v is exactly Z at r, and its profile offset is beta. Thus this is an equivalence of physical marked interfaces, not merely a numerical implication.

## 4. Scope, inverse meaning, and a hostile control

The reduction supplies a graph-level inverse: compute the terminal component of T-r. If an exact unit Z clock is claimed, that component must be the unique singleton v described above. Unit-amplitude obstructions for literal outdegree-one clone sources therefore also obstruct this retained-root implementation. They cannot be evaded merely by changing which of the two vertices is called the root.

This does not rule out a first transient, changing the physical graph at the same scalar count, adaptive port operations, or a nonunit offset. Without the unit condition the terminal component need not be a singleton. For example, take U to be the directed triangle with rows [2,4,1], root 0, and S the four-vertex tournament [8,1,3,6]. Construction (4) gives a physical profile

    (H,e,R,d)=(69,5,69,5), beta=-5.

It has a nontrivial terminal SCC after deleting the root and an exact 2H-5 return. This is not a counterexample to the unit theorem.

An initial transient really can matter. The graph [8,1,3,6], marked at r=2, has (h,e,R,d)=(5,1,9,3). Its Z counts are 5,13,25,49,... . The eventual offset is -1, but the first offset is +3, and the original deletion terminal component need not be the existing singleton of the unit theorem.

There is no claim here that every attainable odd scalar admits a unit Z marking, that unit clocks alone cover the target spectrum, or that the original count can be recovered from an eventual clock after ignoring its initial transient.

## 5. Verification and dependencies

The standalone verifier constructs (4) for arbitrary physical U and S controls and checks every row of (5). It separately checks terminal SCCs, the divisibility (7), the forced singleton for unit states, three literal returning births, and the nonunit H69 example. The control file is an inherited exact representative bank through order seven; no completeness of that bank is needed for the all-order proof.

The only inherited transfer is the physical retained-root Z row from the independently reviewed pair-entry bridge, whose core.py is included unchanged. Its proof snapshot is SHA256 fcff7ca192c4055243aeab2344582686214aa7aa25504d1233f3b6b33465a856. Equations (5)-(8) and the terminal-component reduction are proved here. All computations use standard-library Python and exact integers. Finite checks corroborate the structural proof; they are not a finite-to-infinite premise.

## 6. Two inherited all-order scalar-branch exclusions

Theorem 3 of `session197-prime-atoms/PROOF.md`, frozen SHA256 eb3fd77753f216e36e1ab8ce26c9cbd9185cfe54ddd96acc47503d2c44863f31 and independently reviewed in `session197-prime-atoms-audit/REPORT.md`, SHA256 69cf399e545474befc9fc321de5704da9d5337caaee451fc13524e787cf2e9f7, excludes literal transitive-clone families 7*2^m-1 and 7*2^m+1 at every order.

Applying the physical equivalence (8) therefore proves:

- No rooted H=13 tournament can use repeated Z to realize the exact returning branch H->2H+1 from its original count.
- No rooted H=15 tournament can use repeated Z to realize the exact returning branch H->2H-1 from its original count.

An exact affine return from the original count already forces e=d by Section 1, so no additional unspoken profile hypothesis is needed. The known H13 graph [2,12,9,17,7], marked at r=1, does have the other branch: its profile is (13,3,15,3), and it gives 13,25,49,97,... . Its hidden clone vertex is 0, of sole out-neighbor 1, with amplitude 6 and exponential constant +1.

These are exclusions for the fixed physical Z operation. They do not exclude an isolated first step, a different operation, a graph rebuild, or either scalar count itself. The cited amplitude-seven theorem is used only for these two corollaries; the module identities and equivalence above are independent of it.
