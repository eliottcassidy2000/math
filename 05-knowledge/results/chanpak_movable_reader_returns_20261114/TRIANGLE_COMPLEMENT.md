# A nontransitive center covers targets rejected by every transitive-spine apex

The two-wing formula includes the chart with no left wing,
\[
    N=K+h(h_{\rm bin}(y)-1)+Py+S(y\bmod2).
\]
A background research pass supplied the following small witnesses;
we independently replayed the literal graphs and the all-order
[one-apex inverse](https://github.com/eliottcassidy2000/math/blob/c598c9ba1c079e42d82c61347185b995b03c84cc/05-knowledge/results/chanpak_long_cycle_and_transfer_20261114/HALF_INVERSE.md).

Use the directed triangle with outgoing masks (2,4,1), and a transitive
right wing with contact word 1110000, so y=120. The apex beats either
{0} or {0,1} in the triangle. The shared baseline K=5 and the two exposed
(P,S) pairs are (3,1) and (6,2), with h=3. As h_bin(120)=107,
the two counts are 683 and 1043. A single arc flip changes one
into the other. Full outgoing masks are
\[
\begin{split}
 &(1018,2044,2041,1008,992,960,1920,1792,1536,1024,57);\\
 &(1018,1020,2041,1008,992,960,1920,1792,1536,1024,59).
\end{split}
\]
Both graphs are strong: the apex reaches the whole triangle via
the marked contact, the triangle beats the spine, and every spine
vertex can follow its spine to the last vertex, which beats the apex.
The transitive-spine inverse returned no witness at either target;
both are genuine additions to the old spectrum when a nontransitive
center is allowed. This verifies strict improvement, not a guarantee
that any fixed center bank covers the whole complement.

The displayed one-sided formula is a specialization of
[the two-wing proof](bilateral.md), not an additional transfer theorem.
