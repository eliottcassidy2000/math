# What failed, which structure caused it, and what survives

**PROVED scoped / FINITE-EXACT.** The fixed-root obstruction below is
compatible with the successful changed-port grammar: its hypothesis is
that the distinguished vertex remains the shared vertex.

## Seven coordinates for fixed-root left and right contexts

Set Z=(h,s,e,q,F,G,P), with e the full paths ending at r, G=F in the
converse tournament, and P=P_(T−r)(2). For A⋆B with both operands of
order at least2, let δ_A=G_A−h_A−q_A+e_A. Then

    h' = q_A h_B +(e_A−q_A)s_B −δ_A q_B +(G_A−q_A)F_B,
    s' = s_A q_B,                 e' = q_A e_B,
    q' = q_A q_B,
    F' = P_A F_B +(F_A−P_A)q_B,
    G' = q_A G_B +(G_A−q_A)P_B,
    P' = P_A P_B.

The first four laws are inherited from `surplus_selector_20261113`.
The e and G laws follow by reversing the construction; P follows because
deleting r leaves the ordinal sum(A−r)→(B−r), and the two-bin polynomial
is multiplicative for ordinal sums. Singleton operands are identities,
handled separately. Thus this is a closed seven-coordinate response
algebra for fixed-root contexts on either side. Associativity is also
literal: nonroot blocks are totally ordered from left to right.

Independent controls:7,952 physical small gluings,112 permutation checks,
216 associative checks. These controls supplement the formulas; they do
not establish arbitrary-order completeness by extrapolation.

## A sign cone forbids the hoped-for endpoint-type switch

Call a difference E-neutral if Δh=0 and ΔF=Δs, and S-neutral if Δh=0
and ΔG=Δe. Orient a pair so Δq>0. Consider the cone

    ΔP>0,  ΔG>Δq,  ΔG>Δe.

All eleven inherited literal pairs lie in this cone after choosing their
order. For every rooted tournament q>0, P>0, G>=q>=e: G contains the
first-singleton-root contribution q; deleting a final root injects the
e full paths into the q deletion paths.

Left gluing preserves the cone because

    Δq'=q_A Δq,   ΔP'=P_A ΔP,
    ΔG'−Δe'=q_A(ΔG−Δe)+(G_A−q_A)ΔP,
    ΔG'−Δq'=q_A(ΔG−Δq)+(G_A−q_A)ΔP.

Right gluing preserves it because

    ΔG'−Δe'=Δq(G_B−e_B)+(ΔG−Δq)P_B,
    ΔG'−Δq'=Δq(G_B−q_B)+(ΔG−Δq)P_B,

and both other coordinates are again multiplied by positive factors.
Every displayed strict inequality persists. Therefore NO finite word of
fixed-root left/right gluings makes these pairs S-neutral, even if
intermediate steps are allowed to lose E-neutrality. This is an all-order
obstruction, not the absence of hits through8. A global converse changes
the meaning of the endpoint observer; it does not refute the fixed-gauge
statement. Changing the selected vertex does escape the hypothesis.

Furthermore, right gluing an E-neutral pair with ΔP≠0 gives

    Δ(F'−s') = ΔP(F_B−q_B).

Here F_B=q_B iff r is a source: the singleton first path gives q_B,
and any incoming u→r gives at least one additional cover starting with
the path(u,r). Hence right E-neutrality forces a source-rooted context.
In that case the slope Δs/Δq is unchanged. This explains the unproductive
right-context search without relying on a size cutoff.

## Finite heads and quotient loss

The number J_n of rooted tournament isomorphism classes whose fixed-root
left response admits a nonzero E-neutral direction is

    n:    2  3  4   5   6    7    8
    J_n:  2  3  5  13  33  131  994.

The inherited lower bound J_n>=2a_(n−1)+a_(n−2), n>=3, where a_n=A000568,
gives excesses0,0,3,5,7,26 for n=3..8. The prime-looking short tail does
not persist to the next term. This is a useful refinement sequence, not
a formula for A000568.

There are22 fixed-root input/output slope types cumulatively through8.
None accepts1/3,4/9,6/11. The separate two-port census DOES accept all
three. Its occurrence totals use all roots of chosen ordinary type
representatives and are not isomorphism-class counts. Do not compare
those totals directly with J_n.

The initial numerical direction search discarded output-port incidence.
Adding a port restores composition with four profile coordinates, via
PG1; retaining one root but adding only endpoint scalars cannot do so.
All nontrivial fixed-root accepted slopes in this finite head decrease;
two-port rules can increase them, e.g.1/3→1/2, or become negative.
Thus a monotone-slope principle must retain its fixed-root hypothesis.

## Boundaries that remain after the local repair

- **Unit readiness:** the inherited53,084,231 exclusion is for indegree
  one with defect one. A renewable half-slope pair need not have that
  property. PG4 shows why two nonidentical repaired states cannot both
  stay unit-ready at the same h and E.
- **Induced containment:** every gluing word retains its input as an
  induced subtournament. Incoming
  [fill-tree, induced boundary](../fill_tree_matrix_20261113/INDUCED_BOUNDARY.md)
  proves that no H13 source embeds in any H27 target, at any order. A
  changed output port cannot make13→27 into an induced-extension rule.
  A global target tree needs different parents or genuine surgery.
- **Density:** the nine small half-return matrices have the exact
  O(X^(3/4)) upper bound in [LANGUAGE](LANGUAGE.md). Repeating successful
  repairs is insufficient to populate almost every target.

These are different failures: endpoint-sign obstruction, missing port
incidence, forbidden containment, and insufficient numerical branching.
They require different remedies and must not be merged into one scalar
“bad count” list.
