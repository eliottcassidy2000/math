# Counting lenses are not the native line operator

2026-10-24, Session 52. **PROVED, inherited mechanisms:** the signed-pair character formula, ordinary and self-converse partition counts, and the elementary operator distinctions below. **FINITE-EXACT:** independent replay through the declared bounds. **CORRECTION:** the unguarded historical `sc_partition_formula.md` has a false criterion and must defer to the already corrected THM-283. **NOT CLAIMED:** a new A000568 algorithm, a new self-converse enumeration theorem, a new large-order class census, or recovery of native metagraph wiring from counts.

## 1. Inheritance and the precise question

The useful prior attempts solve different problems:

1. `01-canon/theorems/THM-283-sc-burnside-formula.md` already states the correct self-converse partition filter. Its endpoint-monodromy explanation was repaired and independently audited in `05-knowledge/results/converse_tournament_count_channels_20261021/`. That packet also gives the existing reduction of self-converse counts to two ordinary odd-partition channels. This note re-proves and checks those mechanisms; it does not originate them.
2. `01-canon/theorems/THM-587-per-level-signed-cycle-index-antipodal-euler-number.md` gives a signed character polynomial whose coefficients determine the **single-arc-flip weighted class quotient**. `THM-588-metagraph-laplacian-algebraic-connectivity-four.md` explicitly defines that operator and its spectrum. Its operation is not simultaneous complementation of all non-path arcs.
3. `01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md` already gives class BFS, Hamiltonian-path fibre enumeration, and the checksum `sum H(T)/|Aut(T)| = 2^binom(n-1,2)`. These are inherited infrastructure. Its comparison of the **single-tile** quotient with simple arc adjacency is not a comparison of full tile-complement lines with arc flips.
4. `05-knowledge/results/metagraph_growth_20261020/enumeration.md` already separates the two graph sorts, supplies the minimal order-four non-descent witness, and gives order-three marked-growth and order-seven outgoing-partner losses. The source bytes recovered for this lane are the published copy in `session50-publication/clean-repo/`; the current sparse `math/` checkout does not contain that directory.
5. `05-knowledge/results/tournament_enumeration_20261021a/README.md` already computes A000568 through 200 with the inherited response DP, odd-denominator integer arithmetic and packed keys. It reports a single-host 7.94-second run, not a portable guarantee. Re-grouping the signed trace by partitions here is a small independent checker, not a competitive new counting algorithm.

The closest mechanism is signed pair-orbit monodromy. The hostile is the transposition versus the four-cycle: pair-orbit length without endpoint swapping gives the wrong criterion. The corrected near miss is the old “arbitrary even cycles and fixed points” rule. The least-used sidecar is the literal Hamiltonian-path presentation, which determines where a full tile complement lands.

The live concepts are: ordinary orbit counts / converse-fixed counts / signed Fourier levels / native line incidence / Hamiltonian-path fibres. The anchor is precise separation of these outputs; the niche is a portable coefficient-level replay; the wildcard is the smallest explicit operator mismatch. No arbitrary visual analogy is needed.

## 2. Three operators with different domains

Let `d=binom(n,2)` and let `X_n={0,1}^d` be labelled tournaments. Relabelling acts by signed coordinate permutations. The full arc-reversal involution `R` complements all `d` coordinates and commutes with relabelling.

The labelled hypercube adjacency flips one coordinate. Its quotient on ordinary tournament classes is the weighted row operator

    A_n[i,j] = number of arcs of a representative of class i
               whose reversal has class j.

This is well defined under relabelling and has row sums `d`. It is generally not symmetric in this class basis, although it is reversible with labelled-orbit-size weights. Its spectrum is the restriction of hypercube adjacency to invariant functions.

In contrast, fix the directed Hamiltonian path `0->1->...->n-1`. A tiling has `m=binom(n-1,2)` free non-path arcs. The native involution `kappa` complements **all m free bits simultaneously**, preserving the path. Its unordered pairs are the user's blue/black lines. Projecting endpoints to tournament classes, or to converse-merged classes, produces a line multigraph. The fibre over an ordinary class is `H/|Aut|`, not one tiling in general. The projected lines give a weighted correspondence, not a function on classes.

A third operation flips just one of the `m` tile bits. THM-1430's single-tile adjacency comparison concerns this third operation. Equality of two simple edge supports, when established, would not identify their weights, loops, spectra, or the full-complement pairing.

### Minimal operator mismatch at n=3

Order the two classes as transitive, cyclic. Both are self-converse. Both tilings are reflection-fixed, so the native graph has exactly one blue edge. The two operators are

    native line incidence C_3 = [[0,1],[1,0]], spectrum {1,-1};
    single-arc quotient A_3   = [[2,1],[3,0]], spectrum {3,-1}.

The transitive triangle has two arc flips staying transitive and one becoming cyclic. Every arc flip of the directed cycle becomes transitive. Thus even at the first nontrivial order the native operator is not the arc quotient, including after converse merging. Meanwhile `P_3(t)=1+t^2` correctly gives the arc eigenvalues `3-2k` at `k=0,2`. It does not give the native eigenvalue `1`.

### Minimal failure of kappa descending to classes at n=4

Use free pairs `(0,2),(0,3),(1,3)` and bit 1 for the smaller vertex pointing to the larger. Tilings 0 and 1 have rows

    (2,4,9,3), (6,4,8,3),

and are isomorphic by `0->1->2->3->0`. Their complements are tilings 7 and 6, with score multisets `(0,1,2,3)` and `(0,2,2,2)`. One is transitive; the other and its converse are not. The targets remain different after converse merging. Exhaustion through order three verifies minimality. This is the published October 20 witness, independently replayed here.

**Loss ledger.** Source: fixed-path tilings. Target: (possibly converse-merged) tournament classes. Map: forget the presentation. Preserved: the isomorphism class and its converse orbit. Lost: location inside the Hamiltonian-path fibre. Missing sidecar: a Hamiltonian path modulo automorphisms, equivalently the actual tiling. Cheapest decisive test: the masks 0 and 1 above. The strongest survivor is a well-defined incidence matrix obtained by aggregating all presentations; there is no contradiction with counting its nodes.

## 3. The signed character polynomial and what its values count

For a permutation `sigma`, let its induced cycles on unordered pairs have lengths `ell_c` and orientation-product signs `s_c`. Choosing the smaller-to-larger orientation as coordinate gauge, one relabelling step contributes minus exactly when the order of the two endpoints reverses. Define

    P_n(t) = (1/n!) sum_sigma product_c (1+s_c t^ell_c).

A degree-k Boolean-Fourier character is indexed by a k-element subset of pairs. Its trace under `sigma` vanishes unless that subset is a union of pair cycles; including a whole cycle contributes `s_c`. Therefore the coefficient of `t^k` in the product is the trace on the kth Fourier level. Averaging characters projects onto relabelling-invariant functions, so `[t^k]P_n` is a nonnegative integer multiplicity. The hypercube adjacency is `d-2k` on that level, proving the arc-quotient spectral statement.

At `t=1`, the product counts labelled tournaments fixed by `sigma`; the average is the ordinary class count `a_n=A000568(n)`. At `t=-1`, it counts tournaments fixed by `R sigma`; averaging yields the trace of converse on ordinary class functions. Each self-converse class contributes one, so `P_n(-1)=s_n`, the self-converse count. Equivalently, a self-converse class has an anti-isomorphism coset of size `|Aut|` at each of its `n!/|Aut|` labelled representatives; its total contribution to the twisted sum is `n!`.

Consequently the number of converse-merged nodes is `(a_n+s_n)/2`, and the number of merged nodes representing two ordinary classes is `(a_n-s_n)/2`. These are exact **node counts**. They do not specify native line adjacency or its spectrum. On the separate converse quotient of the arc operator, only even Fourier levels remain, with polynomial `(P_n(t)+P_n(-t))/2`.

The safe meaning of `P_n(-1)` here is an alternating graded trace and a converse fixed-class count. No differential or topological complex is constructed in this packet; no Betti-number conclusion is needed or claimed.

### Grouping by vertex-cycle type

For a partition `lambda=(L_1,...,L_r)` and `z_lambda=product_L L^m_L m_L!`, the same polynomial is the sum of `Q_lambda(t)/z_lambda`, where

    Q_lambda(t) = product_i I_(L_i)(t)
                  product_(i<j) (1+t^lcm(L_i,L_j))^gcd(L_i,L_j),

    I_L(t) = (1+t^L)^((L-1)/2)                          if L is odd,
    I_L(t) = (1+t^L)^(L/2-1) (1-t^(L/2))              if L is even.

Within a vertex cycle, each non-antipodal distance class is one length-L pair orbit of positive sign. An even cycle also has an antipodal orbit of length `L/2` whose accumulated endpoint swap gives sign minus. Between two distinct vertex cycles, there are `gcd` orbits of length `lcm`, all with positive sign. This proves the formula, including repeated parts and singletons.

This grouping uses `p(n)` vertex-cycle types rather than all `n!` permutations for the trace computation. It is a standard conjugacy-class regrouping, not a novelty or complexity breakthrough claim; coefficient arithmetic is not unit cost. The existing response DP is the relevant faster scalar-count baseline.

## 4. Correct anti-automorphism types, necessity and sufficiency

**Inherited theorem, proved here for clarity.** A vertex permutation can reverse every arc of some tournament if and only if it has at most one fixed point and every other cycle length is `2 mod 4`.

**Necessity.** Every tournament automorphism has odd order: if an automorphism had even order, a power would be a nontrivial involution, and one of its transposed vertex pairs would have its arc reversed. If `sigma` is an anti-automorphism, `sigma^2` is an automorphism, hence has odd order. A vertex cycle divisible by four would give `sigma^2` an even-order cycle, so is impossible. If `sigma` had an odd cycle of length `L>=3`, then the odd power `sigma^L` would still be an anti-automorphism but would fix all L vertices on that cycle. An anti-automorphism cannot fix two vertices, because it must reverse their arc. This also excludes two singleton cycles. The remaining non-singleton cycle lengths are exactly twice odd integers.

**Sufficiency and exponent.** On an unordered-pair orbit of length `ell` and product sign `s`, the anti-invariance constraint is consistent exactly when `s=(-1)^ell`; when consistent, it has two orientation assignments. In a cycle of length `L=2u` with u odd, the `L/2-1` non-antipodal orbits have even length L and positive sign, while the antipodal orbit has odd length u and negative sign. Both pass. Between two allowed even cycles, pair orbits have even length and positive sign. Between the optional fixed point and an even cycle they again have even length and positive sign. There is no pair of fixed points. Thus every orbit passes, independently.

Writing the non-singleton lengths as `L_1,...,L_k` and the fixed-point count as `f in {0,1}`, the number of pair orbits is

    e(lambda) = sum_i L_i/2 + sum_(i<j) gcd(L_i,L_j) + f*k.

The number of labelled anti-fixed tournaments is `2^e(lambda)`, and twisted Burnside yields

    s_n = sum_(lambda: f<=1, all other L_i=2 mod 4) 2^e(lambda)/z_lambda.

The empty partition gives `s_0=1`; the singleton gives `s_1=1`. These conventions are consistent with the trace formula.

The ordinary all-odd formula, used only as the known baseline, is

    a_n = sum_(lambda all odd) 2^(sum_i (L_i-1)/2 + sum_(i<j) gcd(L_i,L_j))/z_lambda.

A permutation with an even cycle cannot preserve a tournament because its antipodal pair has negative sign. Every all-odd cycle type passes all pair-orbit constraints. This independently proves the ordinary specialization.

**Smallest hostile controls.** The old rule wrongly admits the identity type `(1,1)` at n=2, but two fixed points give zero anti-fixed tournaments. It also wrongly admits `(4)` at n=4, whose antipodal pair orbit has length two and negative sign, again giving zero. Positive controls are `(2)` with 2 anti-fixed tournaments, `(2,2)` with 16, `(6)` with 8, and `(1,6)` with 16. The correct criterion and these mechanisms were already in THM-283 and the Session22 correction lineage; the newly required action is a warning on the stale historical file.

## 5. Reproduction, independence, and limits

Run from any directory, with only Python's standard library:

    python3 -B /path/to/this/packet/check.py
    python3 -B -O /path/to/this/packet/check.py --receipt /tmp/counting-lenses-optimized.json

Both normal and optimized runs passed **399,751 explicit checks**, and their JSON receipts are byte-identical. Bounds:

- Every vertex-cycle partition for `0<=n<=20`: literal signed pair-cycle traversal versus the closed inventory; ordinary and anti-fixed counts versus the separate direct partition exponent; full coefficient positivity/integrality; both specializations.
- Every permutation in `S_n` for `0<=n<=8`: full signed coefficient averaging versus the cycle-type average, not merely the two endpoint evaluations.
- Every labelled tournament bitmask and every permutation for `0<=n<=5`: literal relabelling versus equality and bitwise-complement equality. This control does not use pair-cycle signs or either partition filter.
- Inherited independently audited Session22 inputs for both scalar vectors through n=20, frozen in `inputs/counts-0-20.json` with their parent receipt hash. This is not a fresh external OEIS retrieval or a newly computed class census.
- Literal n=3 operator matrices; explicit n=4 isomorphism and unequal merged targets; exhaustive minimality through n=3.

`receipt.json` contains all coefficients of `P_n` and the full ordinary/self-converse/merged table through n=20. Example rows are:

| n | ordinary | self-converse | merged |
|---|---:|---:|---:|
| 3 | 2 | 2 | 2 |
| 4 | 4 | 2 | 3 |
| 5 | 12 | 8 | 10 |
| 6 | 56 | 12 | 34 |
| 7 | 456 | 88 | 272 |
| 8 | 6880 | 176 | 3528 |
| 9 | 191536 | 2752 | 97144 |
| 10 | 9733056 | 8784 | 4870920 |

The checker's two implementation paths are independently coded, but share this author. External peer review is not implied. Frozen source hashes and full scope appear in `provenance.json`; `manifest.json` authenticates the deliverables. No source or canonical file was changed by this lane.

## 6. What actually remains useful for the user's native metagraph

The counting problem has already-successful all-order formulas and fast finite-bound implementations. A potentially valuable new route would therefore have to do one of the following, not merely reproduce those counts:

- retain the Hamiltonian-path presentation long enough to compute actual blue/black outgoing targets;
- prove a compressed family of marked node types closes under both native line incidence and cross-order extension/deletion;
- give a genuinely better scalar-count algorithm with a fair comparison against the existing packed response DP;
- recover specific native graph information, such as multiplicity-weighted walks or connectivity, with an explicit map and retained boundary data.

The n=4 witness refutes deterministic complement on unmarked classes, not global counting methods. The n=3 matrix mismatch refutes substitution of the arc spectrum for the native line spectrum, not the signed character formula. Counts, class generation, and wiring should remain three separately typed deliverables.
