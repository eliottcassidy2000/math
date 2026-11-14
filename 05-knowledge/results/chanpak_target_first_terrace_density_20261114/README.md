# First triangle terrace: exact factorization, not exact coverage

**Proved, all orders.** Among strong tournaments with the minimum (n-2)
directed triangles, the distinct ordinary Hamiltonian path counts at most
(X) number (O(X^{49/50})). Yet every odd residue modulo every power
of two has a strong witness on this same terrace, on (O(k^2)) vertices
for modulus (2^k). The obstruction is global, not a missing binary
digit. In particular 7, 21, and 63 all fail *on this terrace*; 63 has
a literal two-hub tournament witness outside it. This does not add a
global exception to Chan–Pak.

Proof in [TERRACE_DENSITY.md](TERRACE_DENSITY.md);
run `python replay.py` for the exact inequalities and small structural
cross-checks. The first triangle terrace classification and interval
packing were established in [Session 337](https://github.com/eliottcassidy2000/math/blob/26519368143bbd88e824d8aba2c86ef2c0adffb0/05-knowledge/results/session337_rooted_line_surfaces_20261114/TERRACE.md),
using the [Grinberg–Stanley odd-cycle expansion](https://arxiv.org/abs/2307.05569).
The entire terrace, not just a fixed-radius region of the line metagraph,
is covered by the argument below.

**Proved, structural.** A strong tournament at triangle excess (d) is
an ordinal sum of minimum-triangle strong tournaments and singletons
with at most (6d) arcs reversed, in descending score order.
This is not yet a spectrum bound at positive excess. The reversals may
connect far-apart blocks and the scalar path counts may then no longer
suffice. See [TERRACE_DENSITY.md](TERRACE_DENSITY.md#5-above-the-terrace).

[One long reversal](SINGLE_WRAP.md) derives the exact new count
when the outer edge of a reservoir is flipped: (H(W(B))=H(B)+L_B(1)).
For ordered sums both coordinates multiply, so choosing the pieces for a
prescribed absolute target requires a sum of two products built from
the same physical pairs. Some residual parity constraints are exact.
