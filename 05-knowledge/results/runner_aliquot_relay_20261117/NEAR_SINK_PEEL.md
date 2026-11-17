# Peel the near-sink before searching for an aliquot inverse

Follow-up to [the aliquot relay](README.md) and [the induced inverse obstruction](NO_UNIVERSAL_INDUCED_INVERSE.md). This is an all-order reduction of the **near-sink** branch, not a theorem about arbitrary contact cuts or Chan–Pak coverage. The endpoint identity and the suffix bound were found in a follow-up research pass; the ordinal-tail reduction below packages the survivors as two strong-core tests.

## Endpoint cancellation

For a tournament with marked root r, let H be its Hamiltonian path count, g the count after deleting r, and E and O the numbers of ordered covers by a nonempty first path starting at r and an optional second path, with respectively even or odd first length. Let B count Hamiltonian paths ending at r. Then

    O − E = B ≤ g.

To see the equality, concatenate the first path with the reverse of the second and keep the cut position k≥1. Along this permutation the arcs before the cut must point forward and those after it backward. The seam itself can point either way. Thus a permutation allowing a cut has all its forward steps before all its backward steps. If the first a steps point forward, the two possible cut positions are a and a+1, retaining only legal positions in [1,n]. They cancel unless a=0; that leaves k=1, the reversal of a Hamiltonian path ending at r. Removing the terminal r gives the injection proving the bound. Also g≤H: each Hamiltonian path of the deletion admits insertion of r at an endpoint or at a forward-to-backward transition in the root's incident orientation.

Use the [existing exact extension laws](../fourteen_axis_20261117/README.md#3-a-common-continuation-keeping-the-full-14). An ordinary sink S sends (H,E,O,g) to (H,E+O,E+O,g); a near-sink N beating only r sends it to (H+2E,E+O−g,E+O,g).

## The last old near-sink cannot do the job

For prime p with 7∤p, the aliquot predecessor count of 7p is H=p+8, in the sense that s(7p)=H. Inverting that arithmetic edge by appending N requires E=3H−28, which is odd. Peel ordinary sinks and near-sinks at the same root r without deleting r. If there is any old N in the resulting chronological suffix, let h be the count just before the last one. At its birth the new count is H, since all later letters are S. Its coordinate obeys

    E_N = 2E_old + B_old − g ≤ 2E_old = H−h.

With no S afterward this is already too small for H≥29. With exactly one S it gives E≤2(H−h)+g≤2H−g<3H−28. With two or more S, E is even. Thus for p≥23 the peeled suffix is either empty or a single S. The argument also works at p=19: the only possible equality in the one-S bound would need h=g=1 and E_old=13. A count-one tournament is transitive (any directed triangle gives at least three paths by inserting the remaining vertices into its three paths); its marked E is zero or a power of two, since the first path can choose any subset of the root's out-neighbors, with each part then in its unique transitive order. Hence 13 is impossible.

## The same filter can be imposed on a strong core

A tournament of count 7p, p≥19 prime, has exactly one nontrivial strong component: otherwise the ordered component counts multiply to 7p via a factor 7, which is forbidden. All old components strictly before the component containing r remain before the unique new component containing r and the appended vertex, so must be singleton components. Discard that transitive prefix. It affects neither H nor (E,O): it cannot lie in the first rooted path, and it attaches to the beginning of the second path in exactly one way, including when the second path was empty.

Let U now be the whole ordinal tail after the old strong component of r. For nonempty U, the new E and O before the final N are both (E_core+O_core)P₂(U)/2. Indeed the old two paths extend via an arbitrary two-slot cover of U, and among those covers the first tail length is equally often even or odd. The same cut cancellation proves this on U, now with k=0 permitted, so both cuts always survive. If |U|≥2, P₂(U)/2 is even: pair each subset S with its complement in the sum of H(U[S])H(U[U−S]); this gives an even number of odd summands by Rédei parity, with the empty count set to one. The inverse needs odd E, so |U|≤1.

Consequently any near-sink aliquot inverse at p≥19 reduces, up to an inert transitive prefix, to **one of two equations on a strong tournament of count H=p+8**:

1. Direct attachment: E=3H−28. The marked core has no old near-sink beating only r.
2. One sink then attachment: E+O=3H−28. Again the marked core has no old near-sink beating only r.

Conversely either strong-core equation gives the required count by the exact laws. This reduces the interface search but does not show either equation is solvable for all remaining primes.

## Literal witnesses for the second branch

Outgoing masks, with bit i denoting vertex i:

| p | Strong core masks | r | (H,E,O,g) | after S | after S,N |
|---|---|---:|---|---|---|
| 31 | (88,73,107,16,6,91,24) | 0 | (39,42,47,19) | (39,89,89,19) | (217,159,178,19) |
| 47 | (104,121,99,100,45,64,16) | 2 | (55,67,70,9) | (55,137,137,9) | (329,265,274,9) |

These profiles, including the appended graphs, were checked using subset path DP. The sevenfold arithmetic lane itself cannot host three successive prime quotients >3: if s(7p)=7q, then p=7q−8≡q+1 (mod 3). The only three-prime chain ending at 3 is 83→13→3, giving 581→91→21, and the last count is forbidden for tournaments. Hence even a successful selector for the near-sink equations would have to leave that lane to recurse indefinitely.
