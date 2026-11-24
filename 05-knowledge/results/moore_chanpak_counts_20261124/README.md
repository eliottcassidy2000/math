# Moore tournament path counts and the limits of fixed-template Chan–Pak grammars

**Scope.** The exact count and rail constructions below do not require a degree-57 graph. Results explicitly about that graph are conditional. Chan–Pak remains open. Here \(H\) is the number of directed Hamilton paths of a tournament, and \(C\) is the number of its directed Hamilton cycles up to rotation. By contrast a spanning undirected graph path counted in both word directions has an even count.

## The new bridge: physical +2 seams, and what they cannot do alone

A pentagon already supplies **six unbounded exact +2 families**. In each row the same physical exterior arc is reversed at every chain length \(k\ge1\), and both tournaments are strong:

| Lower \(H\) | Upper \(H\) | Method |
|---:|---:|:---|
| \(59\cdot2^k-9\) | lower \(+2\) | compressed seven-site pentagon, one prepared arc |
| \(32\cdot3^k+30\cdot2^k-23\) | lower \(+2\) | same compressed host, marked port 0 |
| \(68\cdot3^k-18\cdot2^k-35\) | lower \(+2\) | same host and exterior arcs, opposite marked port 6 |
| \(245121\cdot2^k-74075\) | lower \(+2\) | doubled pentagon, a first apex and an outer apex |
| \(16\cdot9^k-20\cdot6^k+33\cdot4^k-29\cdot2^k+11\) | lower \(+2\) | compressed pentagon, two coupled equal widths |
| \(66\cdot9^k-180\cdot6^k+187\cdot4^k+64\cdot3^k-196\cdot2^k+71\) | lower \(+2\) | same host, a second two-port choice |

For example the first row starts at **109 ↔ 111**; the fourth starts at **416,167 ↔ 416,169**. These are constructive pairs, not phase hits. The middle two have **equal unedited profiles** at opposite blue ports, yet the same exterior edits choose different correction rails. See [compressed ports](moore-compact-count-rails.md), [the two-apex construction](moore-chanpak-fork.md) and [the coupled-width pair](moore-coupled-rails.md). Small pairs in the previously certified range are new adjacency witnesses, not claims of new scalar representability. The repo's separate [rooted-bank run through \(2^{35}\)](../../../01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md) is not replaced by these rails.

**General obstruction.** A [Dirichlet-series proof](bounded-templates.md) shows that an arbitrary finite bank of fixed tournament templates, even with substitutions nested without limit and with ordered joins, produces only \(O_{\epsilon}(X^\epsilon)\) distinct heights up to \(X\), for every \(\epsilon>0\). The proof uses ordered path covers as a profile: strongly connected substitution shifts the profile parameter up by one via a full cyclic lap plus an extra path. A convergent Euler product then controls joins. Thus a single Moore57 construction, or even every fixed-size Moore57 chart as a seed, **cannot be the missing ingredient if the rest of the proof only uses fixed-template substitutions**. The positive seams and this impossibility coexist: a unit difference is not a count selector.

**A two-port equality gate.** In the first genuinely coupled row, \(|\Delta H|=2\) holds **exactly when the two chain widths are equal**; the sign of \(\Delta H-2\) detects which width is larger. For arbitrary fillers \(X,Y\) the exact signed response is
\[
2h_Xh_Y+5(Q_Xh_Y-h_XQ_Y)+4(R_XQ_Y-Q_XR_Y),\quad
(h,Q,R)=(Z_1,Z_2,Z_3).
\]
The antisymmetric terms give a precise correction mechanism absent from single-port rails. For transitive chains \(h=1\), \(Q=2^k\), \(R=3^k\), both ratios \(Q/h\) and \(R/Q\) increase strictly. The first antisymmetric term has absolute value at least 10 for unequal lengths, and the second has the same sign; equality therefore gives \(+2\) exclusively. There are 195 distinct triples through order seven and [none of their ordered pairs with a nontransitive filler gives \(\pm2\)](coupled-fillers.py), despite many pairs with opposite signs in the two correction terms. Also see [the exact two-parameter profile inverse](profile-transfer.md), which tests an integer by Mersenne divisors rather than by phase.

## Compare the actual readouts

These exact values are for the specified charts, not graph invariants. In particular the graph itself has 10 (pentagon) or 240 (Petersen) undirected Hamilton path words including both directions.

| Graph | Tournament chart | Sites | \(H\) | \(C\) |
|---|---|---:|---:|---:|
| Pentagon | graph-spine sparse / dense | 5 | 9 / 11 | 1 / 1 |
| Pentagon | cyclic Gram decoder | 5 | 15 | 2 |
| Pentagon | compressed hidden-star pair | 7 | 135 / 147 | 12 / 12 |
| Pentagon | doubled, order 01432 / 01243 | 10 | 6,461 / 6,821 | 298 / 340 |
| Petersen | graph-spine sparse / dense | 10 | 5,213 / 6,961 | 278 / 376 |
| Petersen | compact affine | 10 | 8,265 | 464 |
| Petersen | compressed hidden-star chart | 16 | 834,157,149 | 27,494,558 |

The full 50-vertex Hoffman–Singleton counts are **not computed** here. Retaining only paths that visit every displayed 5-cohort contiguously gives exact contributions and hence lower bounds: **887,480,820,590,625** in the compact affine chart, **2,496,220,189,141** in the sparse graph-spine chart and **18,811,695,423,864** in its dense mate. The endpoint/subset DP in [counts.py](counts.py) takes arbitrary orders of the ten blocks and exact internal path endpoint tables. It does not assume the Hamilton path uses the displayed global order.

The two direct Petersen sparse graph-spine classes have the same \(H=5213\) and give the same **59 distinct heights** in their six-chord subcubes; neither cube has a literal +2 single-chord move. Yet their first endpoint has split signatures \((5213,8769,3726)\) vs \((5213,7513,2898)\) (they interchange at the last endpoint). Inflating that endpoint with a transitive \(k\)-chain produces respectively
\[
 3726\cdot3^k-2409\cdot2^k-1147,\qquad
 2898\cdot3^k-1181\cdot2^k-1119 .
\]
In the sparse native graph-spine model the first endpoint has indegree \(d-1\), so it can occupy at most \(d\) separated visits after inflation. In particular these port signatures could have up to seven slots at degree 7, or 57 at degree 57 **if** a graph Hamilton spine exists. A complement spine is available at degree 57 and permits all graph edges as backward arcs with at most 58 visits at that endpoint. These are upper bounds, not existence of full support. See [forks.py](forks.py) for the exact visit counter.

**Parity pressure.** For every even self-converse tournament, \(H\equiv1\pmod4\): an involutive arrow reversal has no fixed site and pairs directed odd cycles, so the odd-cycle packing formula modulo 4 gives the congruence. In particular a +2 seam cannot have two even self-converse endpoints. Pair-preserving mate switches are limited to multiples of four in \(H\), though the parity of \(C\) changes by an exact affine rule on the mate bits. See [diagonal fugacity](moore-diagonal-fugacity.md). This is why the useful +2 geometry must be allowed to leave the blue stratum when it starts on an even blue chart.

## How to make degree 57 relevant

The [star-contact and cap calculation](degree57-ports.md) gives a conditional degree-57 test: a 60-port cap would force the exact equitable partition of sizes \((60,1650,1540)\). A root-uniform construction is not justified either: degree 57 would have at least nine vertex orbits. Both suggest looking for distinct, compatible root responses rather than multiplying a single universal count.

The concrete next conjecture is **not** that the larger graph's count fills the gaps. It is that a Moore-like unique-center tiling can coordinate a growing number of *nonuniform* contacts so that signed count transfer closes under an exact target inverse. For a useful first test, find two or more interacting contacts whose joint change is not determined by fixed-module path profiles, while preserving the split signatures needed by a later step. The smaller forks calibrate what exact unit control should look like; the obstruction above says it must be combined with genuinely growing quotient complexity, even if degree 57 is eventually constructed.

## Replay

Run `python counts.py`, `python forks.py`, `python moore-compact-rail-final.py`, `python moore-chanpak-lift-check.py`, `python moore-diagonal-fugacity-check.py`, and `python moore-star-cap-history.py`. Also run the [coupled rail checker](moore-coupled-rails.py), [asymmetric filler census](coupled-fillers.py), [profile replay](moore-chanpak-profile-check.py), and [strong-template transfer check](bounded-templates-check.py). These are finite certificates. The fixed-template obstruction is a separate all-orders proof.

Primary context: [Chan–Pak §6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf); [Grinberg–Stanley](https://arxiv.org/abs/2307.05569). Inheritance: [direct graph-spine results](../moore57_graph_spine_tiles_20261124/), [compact/native routes](../moore57_compact_blue_routes_20261124/), [projective nonlinear pencils](../moore57_nonlinear_tournament_pencils_20261124/).
