# A one-arc conjugate pair that keeps a whole external reader

The [shared-apex bouquet](BOUQUET.md) also gives a local positive
construction. In a five-vertex tournament take the two triangles
0→1→2→0 and 0→3→4→0. Complete the cross arcs:
1→3, 1→4, 3→2, 4→2. Call this T−; reverse just 1→4 for T+.
Now inflate vertex 2 (not an endpoint of the switched arc) by any
nonempty tournament D, all external arcs uniform. Then exactly

    H(T−[2←D])=6 Q(D)−H(D),
    H(T+[2←D])=6 Q(D)+H(D),

where Q is the ordered optional two-path-slot count. In particular
one arc reversal changes H by exactly 2 H(D) for every D. Transitive
D of arbitrary order m preserve an exact +2 flip between
6·2^m−1 and 6·2^m+1. This is not a novelty claim for those numerical
families: previous signed forms and literal attachments appear in
[the seam session](../chanpak_target_address_20261115/seams/SINGLE_SEAM.md).
The receipt here is a single fixed arc outside an arbitrary module.

## A portable way to certify the entire blow-up response

Order the remaining vertices by a permutation. An internal gap
where the direct arc is missing is mandatory: it must receive a
nonempty D-run, and the original two arcs via vertex 2 must exist,
otherwise reject the permutation. An internal gap with both a
direct arc and a two-arc route via 2 is optional. The beginning or
end can be optional when the corresponding one-arc attachment
exists. For b mandatory and a optional gaps that permutation has
signature z^b(1+z)^a. Sum those signatures over the external
permutations. The coefficient at z^k multiplies P_k(D), ordered
covers by k nonempty paths, by concatenating the k runs in their
chosen order. Ignore the constant term for nonempty D. This is
a specialized form of the prior run-word substitution kernel.

For the two specified five-vertex graphs the (b,a) multiplicities
are, by the 24 external permutations:

    signature        (0,2)  (1,1)  (0,1)  (1,0)
    T−                  4      2      1      0
    T+                  4      2      1      2

The resulting polynomials are 5+11z+6z² and 5+13z+6z².
This proves the claim since Q=P₂+2H. Notice what survives:
both optional-cover coefficient and the external constant are
unchanged, and only the one-run coefficient varies. Each inflated tournament is strong:
the five-vertex skeleton is strong and a nonempty module with all
external arcs uniform preserves paths between its vertices through
the skeleton.

For an exact target pair (N,N+G) in this reader, the required interface
is uniquely H(D)=G/2, Q(D)=(N+G/2)/6. Reject nonintegral values. The
[previous fixed-H dyadic catalogue](../chanpak_graceful_ternary_20261115/BERGGREN.md)
gives a finite exact inverse in this lane. For G=2 the module must
be transitive and the lane is sparse; it is a robust arithmetic
move but not an absolute anchor.
