# Independent audit: fundamental cycles and orientations

**Result:** Pass. No correction to `FUNDAMENTAL_CYCLES.md` or `fundamental_cycles.py` is needed.

For a specified Hamilton path P on n>1 vertices, P is a spanning tree and the number of chords is β=m−n+1. Each chord gives a simple fundamental cycle with exactly that chord. Projecting its incidence vector to the chord coordinates gives a standard basis vector, proving independence. For any even-degree edge subset, symmetric-differencing with those fundamental cycles indexed by its chords removes all chords. What remains is an even-degree subgraph of a forest, which must be empty. Thus these vectors also span the F₂ cycle space.

A directed version of P fixes n−1 directions; the β chords independently choose directions, giving 2^β extensions. Reversing P gives a disjoint set of extensions, because n>1. One unoriented path therefore contributes 2^(β+1) orientation–directed-path incidences. Different paths can share an orientation, so the sum is not necessarily the number of successful orientations. As a hostile example, K3 has 3 undirected paths and β=1: 12 incidences but only 8 successful orientations, distributed as 6 with one directed Hamilton path and 2 with three.

## Independent checks

- Constructed square-sum graph edges, enumerated paths with independent bitmask DFS, constructed basis vectors and enumerated every path's 2^β orientation extensions. Counts agree for G15, G16, G17, G23, G25, G26 and G27.
- G15–G17 each have one unoriented path, one chord {1,3}, fundamental cycle 1–3–6–10–15–1, and exactly four successful edge orientations.
- G23 has m=28, β=6, U=3, and 384 successful orientations, each with one directed Hamilton path.
- G25 has m=32, β=8, U=10, with 4,608 orientations containing one directed Hamilton path and 256 containing two. Hence 4,864 successful orientations and 5,120 incidences. Orientations with zero directed paths are counted by subtracting the enumerated successful orientations from 2^m, not by enumerating all 2^m.
- For all 75 simple graphs on n=2,...,5 containing the fixed path 0–1–...–(n−1), directly enumerated all edge subsets, selected the even-degree subsets, and verified equality with the generated fundamental-cycle span.
- Ran the original `fundamental_cycles.py` in normal and optimized modes, successfully and with identical output. Ran independent `audit.py` in normal and optimized modes, successfully and with identical output; all failures use explicit exceptions.

`normal.txt`, `optimized.txt`, `original.normal.txt` and `original.optimized.txt` record the runs. These finite checks corroborate the general argument rather than replacing it. No original source files were modified.
