# Odd-prime Sylow-action bijections for invariant reversal quotients

**Status: PROVED.** This is the same construction on the full generality of the inherited THM-4531 twisted Mallows–Sloane principle. It is a theorem on quotient-object isomorphism classes; “automorphism group” on that side means the stabilizer of the entire quotient object.

Let V=F_2^{E(K_n)} and let U<=V be any S_n-invariant linear subspace. Tournaments form an affine V-space under edge reversal. Let Q=T/U be their quotient by reversals from U, and let

    U^perp={X in V: |X intersection u|=0 mod 2 for every u in U}.

For each odd prime p there is an explicit bijection between

- S_n-orbits of quotient objects q in Q; and
- S_n-orbits of Royle-sign-even graphs X in U^perp,

which preserves the S_n-conjugacy class of a Sylow p-subgroup of the corresponding stabilizer groups. Consequently it preserves, simultaneously, the presence of every embedded p-group symmetry H<=S_n.

## Proof with all changed objects identified

The characters phi_X are constant on a U-reversal class exactly when X is in U^perp. They give a Walsh basis of functions on Q over R=Z[1/2], with orthogonality factor

    c=|Q|=2^(binomial(n,2)-dim U).

The induced S_n-action permutes these characters with exactly the same sign cocycle as in the full tournament case. Taking signed orbit invariants therefore gives an integral matrix M_U in GL_N(R), with columns indexed by sign-even graph orbits inside U^perp and rows indexed by quotient-object orbits. The argument does not divide by n!; only the Boolean-cube Walsh inverse introduces denominators, all powers of two.

Grouping character orthogonality by S_n-orbits gives

    M_U diag(|Aut X|) M_U^t = c diag(|Stab_Sn(q)|).

The general weighted-pairing theorem produces a p-unit matching with equal p-valuations of the two stabilizer orders. Now choose a Sylow p-subgroup P of Stab(q). For every X in U^perp, phi_X is genuinely a function of q, so the signed term epsilon(Y)phi_Y(q) is constant on P-orbits of graph copies exactly as in `SYLOW_ACTION_ADDENDUM.md`. Its nonzero sum modulo p forces a P-fixed graph copy Y. Equal p-valuations make P Sylow in Aut(Y). Thus the matched Sylow actions are conjugate in S_n. The cofactor recursion over F_p specifies the bijection and its inverse. □

## Important special cases

1. U=0 gives tournaments versus all sign-even graphs, the main result.
2. U=cut space gives tournament **switching classes** versus sign-even graphs with all degrees even. The target has both predicates: Royle-sign-even and Eulerian. On the source, the group is the permutation stabilizer of the switching class, which can be larger than the automorphism group of any representative tournament.
3. U=span(K_n) gives tournament types modulo global converse versus sign-even graphs with an even number of edges. The group on the left includes permutations that carry a representative tournament to its converse. For odd p, such a p-subgroup necessarily lies in the orientation-preserving kernel of the possible quotient of order two, consistently with the separate converse-block construction.
4. Sums or intersections of invariant subspaces are allowed without a new argument, provided the target is the actual annihilator U^perp.

The cocycle and Fourier quotient construction are inherited from THM-4531, `01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`, Theorem G. The extra conclusion here is the action-preserving p-unit matching. This does not imply an S_n-equivariant map of labelled quotient objects, or a single bijection preserving every odd prime at once.
