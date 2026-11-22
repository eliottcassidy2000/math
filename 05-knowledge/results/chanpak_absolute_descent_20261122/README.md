# Absolute descent in a fixed directed-edge construction

**2026-11-22, euler.** This session asks for exact integers, not local phases. For a successful descent each edge must describe a real tournament operation, and the register passed to the next edge must carry enough information to certify the target.

| Result | Scope |
|---|---|
| [The backedge chain](BACKCHAIN.md) | **PROVED + FINITE-EXACT.** Fixed backwards arcs with variable transitive gap widths give a nonnegative three-register construction. It compiles strong products, a coupled \(2ab-1\) move, and an exact two-register inverse branch. The full grammar gives 6,437 absolute counts below 100,002 absent from the one-apex transitive chart, with an explicit 21-vertex witness for 99,991. |
| [Strict summand gluing fails at one root](SUMMAND-GLUING.md) | **PROVED.** Retaining disjoint induced operands costs at least their product; sharing one root still overshoots the desired reduced sum by at least four for nontrivial operands. The lower toll is sharp already for two triangles. |
| [Independent contact-switch classification](CONTACT-CARRIES.md) | **PROVED + FINITE-EXACT.** For disjoint one-apex contact blocks, the mixed term is an exact product plus one adjacent correction. The only two-block additive exception is a carry at a touching boundary. No three disjoint ordered switches form an independent additive digit bank. |
| [Exact finite carry descent](AFFINE-CARRY.md) | **PROVED.** A physically certified common-base affine type graph admits exact integer inversion by a finite carry automaton. Regular-language emptiness, combined with finite verified factors, is a sufficient certificate for eventual coverage. Ordinal products of a bank of exponent \(\delta<1\) still have exponent at most \(\delta\), so they cannot repair a thin direct bank. |
| [Dyadic absolute bands and fixed-backedge compiler](DYADIC-BANDS.md) | **PROVED + FINITE-EXACT.** All \(2^m+3\leftrightarrow2^m+5\) for \(m\ge6\) are strong literal flips on \(m+1\) vertices. A 16-factor band gives \(\Theta(m^{16})\) exact offsets on at most \(m+2\) vertices. A bounded total number of reversed edges cannot supply a cofinite spectrum. |
| [Root-wrap decoder and a thin binary tape](ROOT-WRAP.md) | **PROVED.** A signed rooted state yields a bounded complete exact inverse within the grammar. Its two shortest physical moves have density-zero output. |
| [Single-backarc port depth](PORT-CHAIN-DEPTH.md) | **PROVED.** Complete search within a fixed port alphabet has logarithmic block depth and a sharper prime bound. No fixed word returns with slope two even on a prepared unbounded affine state line. |
| [Paley source repair threshold](PALEY-REPAIR.md) | **PROVED + FINITE-EXACT.** Induced binary and quaternary digits fail for Paley order seven regardless of how many vertices are added. To reach 379 and 377 requires exactly one and three old reorientations, respectively. |

## The remaining finite target

The actual [Chan–Pak / bof–Royle question](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf) remains open. There is no claim that the tested edge bank covers all large odds or that matching residues establishes count equality.

A useful next search is to catalogue physical type-to-type edges, keeping the auxiliary registers that make a backwards step valid. For a common-base scalar controller with these guards truly captured by its types, test up front whether it directly covers all but finitely many odd inputs left after dividing out a fixed verified factor bank. The resulting complement is an exact finite-state object, not a statistical sample. In the backchain specifically the new inverse uses an auxiliary odd count: finding a terminal register for every sufficiently large odd height is still the missing step. The summand obstruction says not to look for plain reduced sums by gluing two nontrivial induced operands with only one shared root.


**Replays.** Run the `.py` checkers with Python 3; `backedge-absolute-check.py` uses `backedge-chart.py` from the same directory. Compile `binary-old-two.cpp` with C++20. `cofinite-decode.py` takes an ordinary height (for example `29`) and emits a literal rooted witness if found.
