# Exact count of Hamilton-count-neutral vertex arrivals

**PROVED elementary corollary of the audited SCC equality criterion.** The
finite replay is separate. This counts labelled extensions of one fixed
induced parent, not unlabelled types or automorphism orbits.

Let T be a fixed tournament on n named vertices with k strongly connected
components, ordered C1->...->Ck. Add m new, distinct named vertices; old
arcs remain unchanged. There are 2^(mn+binom(m,2)) possible extensions W.
Then precisely

    m! binom(m+k,k)                                         (1)

of these satisfy H(W)=H(T). Consequently a uniformly random such extension
is count-neutral with probability

    m! binom(m+k,k) / 2^(mn+binom(m,2)).                    (2)

For m=1, the neutral cuts number k+1: the new vertex lies before all the
components, between two consecutive components, or after all of them.
In particular a strong parent has exactly two neutral cuts, while a
transitive n-vertex parent has n+1.

**Proof.** The OCF equality criterion in PROOF.md says that an
induced deletion preserves H exactly when every removed vertex is a
singleton strong component of the larger tournament. Therefore H(W)=H(T)
if and only if all m added vertices are singleton components of W.
The original strong components stay intact and in their original order.
The new singleton components can be placed in the k+1 gaps between and
outside them. Choose a weak composition of m into k+1 gap sizes, in
binom(m+k,k) ways, and a linear order of the m named new vertices, in m!
ways. Both determine every new arc, and they are uniquely recovered from
W's ordered strong-component decomposition. There is no overcount. Divide
by the number of all independent new-edge orientations for (2). QED.

The convention k=0 for the empty parent gives m! transitive labelled
extensions, consistent with H(empty)=1. For m=0 the count is1.

**Preserved and lost information.** The neutral **number** depends only on
k and m, but the actual neutral cuts require the ordered SCC partition.
The value H(T) alone does not recover that partition, and the formula says
nothing about the distribution of positive increments. For example, source
padding leaves H unchanged but increases k, and hence increases the number
of neutral one-vertex cuts at the new order.

No new OCF, SCC factorization or general random-moment formula is claimed.
This is an explicit counting consequence relevant to order growth.
