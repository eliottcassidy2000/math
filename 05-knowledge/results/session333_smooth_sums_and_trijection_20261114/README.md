# Session 333: exact smooth sums, coupled frame weights, and disconnected carriers

This session resolves two versions of the [Wave 6](../chanpak_waves_20261113/README.md) smooth-sum proposal in opposite directions: the arithmetic antichain lemma is true, and the naive literal compiler into a list of two-hub frames is false. There are nonetheless new absolute counts, not merely modular coverage.

## What can be built at the prescribed count?

- [Absolute smooth rail](ABSOLUTE_SMOOTH_RAIL.md): **every \(1+2^A3^B\) with \(A\ge3,B\ge0\) is an actual unmarked strong Hamiltonian-path count.** In the two-hub sector \(B\ge3\), the spine word \(11\,0^k\,2^\ell\,311\) gives exactly \(1+36\cdot2^k3^\ell\). Block occupancy proves the finite-dimensional specialization, so nine direct Hamiltonian counts certify the whole identity. Smaller \(B\) use literal near-transitive or two-visit inflations. The result is a whole unbounded two-parameter sector, though still a sparse one.
- [Exact endpoint decoder and new slopes](ENDPOINT_FRAMES_AND_SLOPES.md): every positive \(N\) can be realized exactly as the number of Hamiltonian paths **between specified endpoints** in a strong one-apex tournament on \(\lfloor\log_2 N\rfloor+5\) vertices. There is also genuine *unmarked* progress: the earlier offset-one predecessor chain \(15\to45\) of slope parameters continues by literal strong frames through \(135,405,1215\), with the dyadic exponent decreasing by one per step. At parameter \(a=9\), the absolute counts are \(961,1441,2161,3241,4861\). No extrapolation to all further tripled slopes is made.
- [Absolute marked targets](ABSOLUTE_MARKED_TARGETS.md): an arbitrary positive target is also an exact three-run marked path count and the exact number of perfect matchings of a simple disconnected bipartite graph with one fewer edge than vertex, all using \(O(\log N)\) vertices. At a tight uniform tournament cut, that matching count appears as a normalized highest-layer coefficient, but lower layers and their factorial constraints cannot be erased.

## Why the simple additive converse fails

[Two-hub sums](TWO_HUB_SUMS.md) proves the positive antichain decomposition by a decreasing arithmetic inverse, and the forward physical frame decomposition with careful endpoint conventions. Every nontransitive two-hub graph has **at least two weight-one frames**: insert the hubs into the fixed increasing spine path. Thus the physical frame multiset is not a divisibility antichain, even after excluding a designated baseline. Triply covered free sites are precisely an open interval common to all three runs, and the opposite hub order creates the same weighted frame. Moreover the companion frames cannot be dropped: [sharp ternary cost](TWO_HUB_FRAME_DEBT.md) proves that the minimum count admitting \(\beta\) ternary sites is
\[
2\cdot3^\beta+\Theta((\sqrt6)^\beta).
\]
An explicit strong family attains that correction scale. [Many-hub debt](MANY_HUB_DEBT.md) extends the unavoidable multiplicity to \(k!\) for \(k\) hubs and gives the sharp next scale for fixed \(k\). On the other hand, a bounded number of *spine colour runs* is provably too sparse for cutoff, and almost every target requires linearly many such runs in its log-size budget if the two-hub model is to realize it.

## What the trijection and the small numbers say

[Functional trijection and root returns](FUNCTIONAL_TRIJECTION_AND_ROOT_RETURNS.md) extends the chart to directed functional pseudoforests with one tree and possibly several unicyclic components. Odd stabilizers give an exact fixed-point Euler test for a chart; arbitrary disconnected \(n-1\)-edge graphs fail it already at seven vertices. [Smooth endpoint fibres](SMOOTH_ENDPOINT_FIBRES.md) gives a literal path-to-two-cycles-plus-path bijection on a *single anchored frame*. It preserves the count and all actual ternary choices, but only when the required closers are present.

[34,45,69](NUMBERS_34_45_69.md) tracks the older finite seam lower edge, the order-eight minimum, and the new links: one marked-source change sends 45 to 69 with the same Fermat block and size; a different 69 graph literally has two-visit coefficient 34; and the totient pair \(45+\varphi(45)=69\), \(1+\varphi(69)=45\) is rigid in its natural 3/5 factor family. The complete low-count ternary census includes 67 and 69 but not 45. A second ternary site costs at least 105 even in an arbitrarily large host.

The next Chan–Pak step is not another residue proof. It is a *positive regrouping or an unmarked decoder* that assembles multiple prescribed terms and retains only the correct other frames. This session proves that one large smooth term can be put at an exact target, that arbitrary marked targets are real, and that the literal antichain and fixed-run shortcuts cannot provide the missing assembly. The full conjecture remains open.
