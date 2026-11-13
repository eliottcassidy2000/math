# Physical neutral returns in the indegree-one constructor language

Status: **PROVED** algebraic constructions; **VERIFIED** independent literal path counts. Scope: the two explicitly defined constructors below, not arbitrary tournament extensions. Chan-Pak remains **OPEN**.

## Inheritance and target

The closest mechanism is the indegree-one two-birth law in `../../berggren_target_return_20261111/reset/TARGET_RETURN.md`: a rooted state of total count h and defect b advances to count 3h+2b. The hostile boundary is that an equal Hamiltonian count does not determine future counts, and that repairing a conditioned cycle count does not control the total. The corrected near miss is a one-shot retarget described below: it can preserve the full count but cannot repeat immediately. The least-used sidecar here is deletion count d; it closes an exact three-coordinate physical constructor language.

This note supplies (i) an all-source, same-full-count four-arc macro rewrite retaining the literal final root and its predecessor; (ii) arbitrary-depth neutral choices with growing defect; and (iii) diamonds that rejoin at a complete continuation state. Its output counts have a factor 51. It does **not** supply general target-preserving 3h+2/3h+4 returns or a proof of Chan-Pak.

## State, physical constructors, and proof

Rows encode out-neighbor bitsets. Let r have indegree exactly one. Define h=H(V), d=H(V-r), s=the number of Hamilton paths starting at r, and F=the number of ordered two-path covers whose first nonempty path ends at r, allowing the second path to be empty. The existing rooted-cover identity is F=h+d-s. Put b=F-h=d-s and x=h-d. Thus h=x+d.

The integer state lies in x>=0, 0<=b<=d, d>=1, x even, d odd, and x+b=h-s>=1. These are necessary invariant cone conditions; sufficiency for tournament realizability is not asserted. For the last inequality, if a is the sole vertex beating r, then a,r followed by any Hamilton path on V minus {a,r} is a Hamilton path not starting at r.

Define **L** by adding a vertex z beaten only by r, then marking z. Define **T** by adding a vertex a beaten only by r, then a vertex z also beaten only by the *same old r*; z beats a, and z is marked. T costs two vertices and L one. All old names and arcs survive. Every output root has indegree one. A strong input remains strong, because each new vertex lies on a directed cycle through the old strong component. An arbitrary weak input need not become strong: the transitive triangle (6,4,0), rooted at its middle vertex 1, retains vertex 2 as a global sink after every L/T word. This is the smallest hostile control, since the sole indegree-one root of the transitive pair becomes strong after L.

For an arbitrary birth beaten only by r, splitting Hamilton paths at the newborn gives h'=h+F-s. The new root has deletion count h and start count h-s. For an indegree-one input this yields

    L(x,b,d) = (x+2b, d-b, x+d).

For T, after its first birth the old r still has indegree one. Every Hamilton path starting at old r must next use the newborn a, since a has no other predecessor. Its start count is therefore d. Deleting r leaves a global source a over V-r, so its deletion count is also d and its defect is zero. The second birth consequently gives

    T(x,b,d) = (2x+2b, d, 2x+2b+d).

These arguments prove the formulas for every indegree-one source, not merely the tested orders. The two-step advance S=L^2 is

    S(x,b,d) = (x+2d, x+b, 2x+2b+d),
    H(S)=3h+2b.

The invariant cone is preserved by both matrices. In projective coordinates t=x/d,u=b/d, their maps are

    L:(t,u) -> ((t+2u)/(t+1), (1-u)/(t+1)),
    T:(t,u) -> ((2t+2u)/(2t+2u+1), 1/(2t+2u+1)).

T sends every input to t+u=1, and T^2 sends every input to the ray (x,b,d)=(2,1,3)q, where q=2x+2b+d=H(L(V,r)). In particular T^3=3T^2, and the projective ray (2/3,1/3) is fixed. This is an exact rank collapse of a physically realized language, not a claim that all abstract cone points have realizers.

## Universal full-count neutral rewrite

Read words chronologically. Set A=TTTLL and B=TTLTL. Each adds eight vertices. Direct matrix multiplication gives

    A(x,b,d)=(24,9,27)q,
    B(x,b,d)=(22,12,29)q,
    q=2x+2b+d.

Both total counts equal 51q, while their defects are 9q and 12q. The initial TT is common. From its last root r0, write the last four new vertices a,b,c,d. Their sole predecessors at birth are:

| suffix | a | b | c | d |
|---|---|---|---|---|
| TLL | r0 | r0 | b | c |
| LTL | r0 | a | a | c |

Exactly four arcs differ: r0-b, a-b, a-c, b-c. Toggling those four arcs simultaneously is an involution between the two outputs. Every older vertex and arc, the four newest names, the final root d, and its sole predecessor c are retained literally. The root used by the *construction* changes at births; the root used by this *rewrite* is the same physical final vertex on both sides. These are different assertions and neither is silently substituted for the other.

The simultaneous four-arc rewrite preserves the full Hamiltonian count. No claim is made that four successive single-arc flips each preserve it. An arc-by-arc metagraph must retain the intermediate count excursions or treat this operation explicitly as a macro.

Before the TT preparation, the suffix TLL has H-row (12,12,5), whereas LTL has H-row (9,6,9). Thus the same rewrite is neutral exactly on

    3x+6b=4d.

Within the nonnegative projective cone this is the segment t+2u=4/3, 0<=u<=2/3. TT guarantees a realizable point on this segment for every source. Realizability of every other lattice point on the segment is not asserted.

## Infinite composition and complete-state diamonds

Write vA=(24,9,27), vB=(22,12,29), and rho=(2,2,1). Then rho*vA=93 and rho*vB=97. A prepared state q*vA or q*vB therefore obeys the exact transition table

| incoming type | next A | next B |
|---|---|---|
| A | 93q*vA | 93q*vB |
| B | 97q*vA | 97q*vB |

At fixed q the neutral rewrite chooses either outgoing type without changing H=51q. It can therefore choose which multiplier, 93 or 97, the next eight-vertex macro uses. This rule repeats indefinitely, and defects grow with q. The four-arc switch itself adds no vertices.

For a fixed original source q0 and a word of m>=1 macros, the state is

    q0 * 93^a * 97^b * vE,
    a+b=m-1, E in {A,B}.

The integers a,b count A and B among the first m-1 macros; E is the last macro. Since 93=3*31 and 97 is coprime to 93, these give exactly 2m distinct states and m distinct Hamiltonian counts at that construction depth. Each state has binomial(m-1,a) formal macro words. This is a census of words and continuation states, **not** a count of nonisomorphic tournaments; a tournament could admit more than one construction address.

The first nontrivial diamonds are the exact matrix identities

    ABA = BAA,
    ABB = BAB.

For example ABA and BAA both give q0*93*97*vA. Their total count is 460071q0. They retain the original input as an induced labeled subtournament and use the same new-label budget and final root. Every later word in L,T gives identical full Hamiltonian counts from either endpoint. This is a stronger reentry than the single A/B neutral switch, whose next values of q differ.

The coordinate state is precisely observable in this language. Knowing H before extension, H after L, and H after LL gives

    b=(H_LL-3H)/2,
    d=2H+2b-H_L,
    x=H-d.

Thus two states have all future L/T counts equal if and only if their (x,b,d) coordinates agree. This does not imply equality under arbitrary external tournament attachments: such contexts may require a larger boundary tensor.

The factor 51 is permanent in this prepared macro family. These neutral diamonds consequently prove a coherent infinite metagraph family, not odd-count coverage or the desired affine target return. The family shows how growing defect can support exact repeatable neutral choices; it does not lower defect to one.

## Why the simplest retarget does not iterate

Compare S=L^2 with T. Their total-count difference is x+2b-2d=F-s-2d. Retargeting the latest vertex from its immediate newborn predecessor to the old root changes exactly two incident arcs and is full-count neutral precisely when x+2b=2d. Under that guard, the output T has defect d instead of S's defect x+b.

The next application of the same guarded repair would need x'+2b'-2d'=0. Substitution of T gives -2(x+b), which is strictly negative because h-s=x+b>=1. Hence this particular neutral retarget can never immediately repeat on its own. This is a theorem about one exact physical rule, not a prohibition on other sparse repairs or growing-state constructions.

## Verification and boundaries

Reproduce:

    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/constructor_language.py

`constructor_language.json` records 100 independent literal state checks (L,T,LL,TT on all 25 indegree-one interfaces in inherited ordinary tournament classes through order five), eight direct full-path-count A/B pairs through order thirteen, exact changed-arc lists, and symbolic diamond matrices. Strongness is checked under the correct strong-input premise. The source bank is the previously certified `../../metagraph_aspect_sequences_20261111/records.jsonl.gz`; no new isomorphism classification is inferred from this check.

An exhaustive matrix-word enumeration by vertex cost finds no equal-H/different-b pair at costs 1 through 7, and the first at cost 8 is A/B. At cost 7, TLTLL and TLLLT already have equal H and b but different d, a direct hostile control for dropping the deletion sidecar. This minimum statement is finite and restricted to the declared positive L/T word language.

Independent audit by the full-count lane rederived L, T, the rank collapse, both macros, both diamonds, and the non-repetition obstruction from physical path decompositions. Computation corroborates these proofs; it does not replace their quantifiers.
