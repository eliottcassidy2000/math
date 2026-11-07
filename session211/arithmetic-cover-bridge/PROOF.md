# Comparison defects, arithmetic operation graphs, and a common path-cover algebra

**Status: PROVED elementary interfaces and obstructions; FINITE-EXACT controls as specified below.** No literature-priority claim. The tournament spectrum conjecture remains OPEN. The cover reciprocity mechanism is classical in character; the purpose here is a typed bridge and explicit hostile controls, not a claim to have invented the cover polynomial.

## 1. Inheritance and the actual objects

The pinned source snapshot is `214586af6c8803d8b21959951e8f151e64bff36f` of `eliottcassidy2000/math`. Source manifest: [SOURCES.json](SOURCES.json).

The closest proved mechanism is Session210's [gate and ordered-cover correspondence](https://github.com/eliottcassidy2000/math/blob/214586af6c8803d8b21959951e8f151e64bff36f/session210/graded-spectrum-gates/PROOF.md). It already proves the tournament gate operator and finite observation theorem. We extend the operator to loopless digraphs, then determine which comparison information that extension cannot restore.

Two earlier `+2`/doubling objects must be separated:

- The literal graph in [the 15/16 seam](https://github.com/eliottcassidy2000/math/blob/214586af6c8803d8b21959951e8f151e64bff36f/05-knowledge/results/seam_threshold_20260925.md) has edges `{x,x+2}` and `{x,2x}` on `1,...,N`. Its undirected planarity threshold is inherited, not re-proved here.
- The [odd-core/dyadic-height address](https://github.com/eliottcassidy2000/math/blob/214586af6c8803d8b21959951e8f151e64bff36f/07-reflections/lrc-add-mult-gauge-stack-s513.md) writes `x=2^h o` with `o` odd. Horizontal `o -> o+2` means `x -> x+2^(h+1)`, not literal `x+2` except on the odd row. This historical reflection supplies object provenance, not a new theorem dependency.

The source atlas identifies `{3,10}` as the odd indices `(H-1)/2` and cardinalities of signed odd-cycle-packing orbits; it does not establish an independently classified number-theoretic family with those two gaps. The square-sum labels `15,25` are graph-order milestones, not those indices. We preserve these types throughout.

The canonical hostile is the loss of missing comparisons when a sparse arithmetic graph is completed to a tournament. The corrected near miss is treating a gate's local degrees as sufficient for Camion's theorem on an arbitrary digraph. The least-used sidecar is the count of unordered pairs with zero or two directed arcs.

Concept board: operation addresses; missing/double comparisons; reversal; ordered path covers; local gates; parallel-coherent permutohedra. The anchor is a map that actually preserves path counts; the niche is the arithmetic graph's directed-versus-undirected boundary; the wildcard is whether every nested gate observation can detect that boundary.

## 2. Ordered covers for arbitrary loopless digraphs

A loopless simple digraph `D` may have neither, one, or both arcs on an unordered pair. Let `H(D)` count directed Hamiltonian vertex words, with the algebraic convention `H(empty)=1`. Let `P_D(q)` count ordered `q`-tuples of possibly empty directed paths that partition the vertex set. In particular `P_D(1)=H(D)`.

For a permutation word `pi`, let `f_D(pi)` be the number of adjacent pairs that are arcs of `D`. Let `D^c` be the arc complement: for distinct `u,v`, the arc `u->v` is in `D^c` exactly when it is absent from `D`. This is not the underlying undirected complement of a general digraph.

**Proposition 1 (common algebra).** For `n=|V(D)|`,

    P_D(q) = sum_pi binom(q+f_D(pi), n),
    P_D(-q) = (-1)^n P_(D^c)(q).                         (1)

For nonempty `D`, this is a monic rational polynomial of degree `n` with zero constant coefficient. If `pc(D)` is the minimum number of nonempty vertex-disjoint directed paths covering all vertices, then, for every integer `q>=0`,

    P_D(q)>0  if and only if  q>=pc(D).                 (2)

**Proof.** Concatenate the nonempty bins of a cover into `pi`. Every absent adjacency must be cut. With `k` nonempty bins there are `binom(f_D(pi),n-k)` ways to choose the retained adjacencies and `binom(q,k)` increasing bin labels. Vandermonde summation gives (1)'s first identity. Each of the `n!` words contributes leading coefficient `1/n!`; at `q=0` every term vanishes for `n>0`. For the second identity use

    binom(-q+f,n)=(-1)^n binom(q+n-1-f,n),
    f_(D^c)(pi)=n-1-f_D(pi).

The empty digraph has `P=1` and satisfies both identities separately. Equation (2) follows by putting a minimum cover in distinct bins, and conversely deleting empty bins from any counted cover. QED.

This gives a common invariant for sparse arithmetic graphs, symmetric square-sum digraphs, and tournaments. It does not make their path-count spectra equal.

For a tournament `T`, `T^c=T^op`, and reversing every path gives `P_(T^op)=P_T`. Hence (1) specializes to tournament parity `P_T(-q)=(-1)^n P_T(q)`. For a symmetric digraph, complement is also symmetric, but it is generally a different graph. Reversing a Hamiltonian word acts freely for `n>=2`, so its `H` is **even**, rather than the tournament's inherited odd `H`.

In the symmetric square-sum case, `H=2p`, where `p` is the number of unoriented Hamilton paths. The independently checked values `p(Q_23)=3` and `p(Q_25)=10` are supplied by the companion square-sum packet, not by this verifier. Its port-refined birth rule retains which two old vertices become endpoints in distinct components of a spanning two-path forest; aggregate cover values alone do not retain those specified ports. Together with the inherited tournament omissions `H=7,21`, they prohibit a universal construction with `(H(T)-1)/2=p(Q)`. The actual common algebra uses `2p`, not `1+2p`.

## 3. The gate works beyond tournaments, but preserves the defect

For `D` define `Gamma(D)` by adjoining `a,b` and all arcs

    a -> b -> D -> a,

leaving all original arcs unchanged. For a nonempty core this is strong: a core vertex reaches `b` through `a`, `b` reaches every core vertex, and any two core vertices are connected through `a,b`. Every Hamiltonian cycle must traverse `a->b`; deleting those two vertices bijects cycles modulo cyclic rotation with Hamiltonian words of `D`. Therefore

    HC(Gamma(D))=H(D)  for nonempty D.                 (3)

For core order at least two the degree-defined gate is unique. Indeed `b` has outdegree `n`, `a` has indegree `n`, and an interior vertex with outdegree one must point only to `a`, which cannot be a gate head. The empty core gives a transitive pair, not a strong digraph, and its cycle count is zero despite the convention `H(empty)=1`.

For `i=0,1,2`, let `d_i(D)` count unordered pairs carrying exactly `i` arcs. Set `delta(D)=d_0(D)+d_2(D)`, so `delta(D)=0` precisely when `D` is a tournament.

**Proposition 2 (comparison-defect conservation).**

    d_0(Gamma(D))=d_0(D),
    d_2(Gamma(D))=d_2(D),
    d_1(Gamma(D))=d_1(D)+2n+1.                         (4)

Consequently no iterate of the gate turns a nontournament into a tournament. Complement swaps `d_0,d_2`, and `Gamma(D)^c` is isomorphic to `Gamma(D^c)` by exchanging the two gate vertices.

**Proof.** Every newly added unordered pair carries exactly one arc, and no old pair changes. There are `2n+1` new pairs. The complement claim follows by writing the three cross directions. QED.

A local gate therefore transports a count without repairing the comparison alphabet. In particular, strongness of a general gated digraph neither forces a Hamiltonian cycle nor the tournament's oddness or forbidden-value laws.

**Proposition 3 (universal gate operator).** For every loopless core, including the empty one,

    P_Gamma(D)(q) = q[P_D(q+1)+q P_D(q)-P_D(q-1)].       (5)

**Proof.** Work in the commutative square-zero ring with variables `x_v`, and write

    F = sum_(S subset V) H(D[S]) x_S.

Then `P_D(q)=[x_V]F^q`. A path using `a` but not `b` ends at `a`; one using `b` but not `a` starts at `b`. Thus both corresponding subset series are `F`. For a nonempty core subset, paths using both gate vertices have count `H(D[S])+P_(D[S])(2)`: either `b` and `a` are the two endpoints, or `a,b` are consecutive, with arbitrary left/right core segments. At the empty subset the count is one, not two. Therefore the whole subset series of `Gamma(D)` is

    F_Gamma = F + x_a F + x_b F + x_a x_b(F+F^2-1).

Taking the `x_a x_b` coefficient of its `q`th power gives

    q F^(q-1)(F+F^2-1)+q(q-1)F^q
      = q[F^(q+1)+q F^q-F^(q-1)].

Finally extract `x_V`. This proves (5) at all positive integers, hence as a polynomial identity. For `q=0` both sides vanish. In particular `L(1)=q^2`, correctly giving the empty-core transitive pair. QED.

The same triangular observation argument as in Session210 now shows that `n-1` observations `H(D),...,H(Gamma^(n-2)(D))` suffice to recover `P_D` for a known nonempty order `n`; no observations are needed when `n=1`. After `r` gate operators, the coefficient of `P_D(r+1)` is `r!`, since reaching `r+1` from one requires all upward shifts. The observations recover values at `1,...,n-1`; together with zero constant coefficient and known leading term, interpolation is unique. This is a sufficient invariant-recovery bound, not a minimum-query claim. Tournament parity lowers the sufficient bound as in the inherited packet.

## 4. A literal arithmetic graph gives three different macrocosms

Let `D_N` have vertices `1,...,N` and arcs `x->x+2`, `x->2x` when their targets are at most `N`. Repeated arcs are counted once. Let `U_N` replace every arc by both directions, and let `T_N` be the numeric transitive tournament, with every arc `x->y` for `x<y`.

**Theorem 4 (all-order typed comparison).** For every `N>=3`,

    H(D_N)=0,   pc(D_N)=2;
    H(U_N) is positive and even,   pc(U_N)=1;
    H(T_N)=1,   pc(T_N)=1.                              (6)

Moreover `P_(D_N)(2)=4`. All three gates are strong, but their Hamiltonian cycle counts are respectively zero, positive even, and one. The arithmetic gate `Gamma(D_N)` has exactly **four Hamiltonian paths** for every `N>=3`, despite having no Hamiltonian cycle.

**Proof.** Every `D_N` arc increases its label. A Hamiltonian word, if present, must be exactly `1,2,...,N`, but the required arc `2->3` is absent. On the other hand, the odd labels in increasing order and the even labels in increasing order form a two-path cover using only `+2` arcs. This proves the first line.

In fact there are exactly two unordered two-path covers. Vertex `1` must start a path. If `2` starts the other, that other path can only use evens, so the path starting at `1` must contain all odds, in their forced `+2` order; the even path is then forced too. If `2` does not start the other path, its only predecessor is `1`, so the first path starts `1,2` and thereafter uses only evens. The other path must start `3` and contain all remaining odds in `+2` order. In either case the odd path ends at the largest odd label, which has no doubling successor within `1,...,N`, so it cannot take any evens after the odd segment. Thus the covers are exactly the two displayed constructions. Ordering their two bins gives `P_(D_N)(2)=4`. Setting `q=1` in (5) gives `H(Gamma(D_N))=H(D_N)+P_(D_N)(2)=4`.

For the symmetric graph take the word

    largest odd, ..., 5,3,1,2,4,6,..., largest even.

Every step is a symmetric `+2` edge except `1--2`, which is a doubling edge. This is a Hamiltonian word. Reversal gives a free involution, proving positive evenness. The unique increasing word is the only Hamiltonian word of `T_N`. Equation (3) proves the gate statement. QED.

The exact number of arcs of `D_N` is

    max(0,N-2)+floor(N/2)-1_(N>=4),

since the two operation types overlap only at `2->4`. Hence its defect is `binom(N,2)` minus this count. For `U_N`, every pair is either empty or bidirected, so `delta(U_N)=binom(N,2)`. For `T_N` the defect is zero. These defects remain unchanged under every gate iteration.

Both small cores `N=1,2` have `H(D_N)=1`. This boundary is required: the obstruction begins at three. Forgetting direction creates the explicit arithmetic Hamilton path; completing missing comparisons instead creates the single transitive path. Neither change is count preserving.

The operation address itself has a second exact warning. If `A(x)=x+2` and `B(x)=2x`, then

    B A = A^2 B,

not `BA=AB`. By contrast, the row-scaled horizontal operation `R(x)=x+2^(v_2(x)+1)` satisfies `BR=RB`, because `R` changes the odd core by two and preserves dyadic height. Thus the literal affine graph and the rectangular odd-core grid have different operation laws. Their shared words `add two` and `double` are insufficient to identify them.

## 5. Every nested cover observation can still miss the defect

**Theorem 5 (minimal comparison-alphabet hostile).** On vertices `0,1,2`, let `E` have exactly the arcs

    0->1, 1->0, 1->2.

Let `T` be a transitive three-tournament. Then

    P_E(q)=P_T(q)=q^3,
    (d_0(E),d_2(E))=(1,1),
    (d_0(T),d_2(T))=(0,0).                              (7)

All nested gate polynomials and Hamiltonian-path counts agree, and all positive-depth gate Hamiltonian-cycle counts agree, although one tower consists of tournaments and the other never does. Order three is the smallest possible order for this failure.

**Proof.** The six words of `E` have forward-adjacency distribution `(1,4,1)`: the sole Hamiltonian word is `0,1,2`, the sole zero-forward word is `0,2,1`, and the other four have one forward adjacency. The transitive three-tournament has the same distribution. Thus

    P_E(q)=binom(q,3)+4binom(q+1,3)+binom(q+2,3)=q^3.

Equality of all subsequent polynomials follows from (5). Equations (3) and (4) give the cycle and defect statements. For two vertices the three comparison types have polynomials `q(q-1)`, `q^2`, and `q(q+1)` for zero, one, or two arcs. They are distinct, and smaller orders have no defect. QED.

This is stronger than equal scalar counts: even the complete cover invariant and the entire nested response sequence fail to determine whether the original object was a tournament. In particular, polynomial parity is necessary for tournament covers but is not a tournament-recognition criterion.

The face incidence expression `sum_pi(1+z)^f_D(pi)` agrees for the hostile pair too. Its interpretation as actual sinks of a parallel-coherent **orientation** is only justified for the tournament member. For an unordered pair with both digraph arcs, both adjacent orders are locally allowed; with neither arc, both are locally forbidden. A single ordinary oriented permutohedron edge cannot have respectively two local sinks or zero local sinks on that one-dimensional face. Thus extending the formula does not extend the orientation equivalence. The comparison sidecar must remain attached.

## 6. Reproduction and scope

Run with standard Python only:

    python3 verify.py > RESULT.json
    python3 -O verify.py > RESULT.optimized.json
    cmp RESULT.json RESULT.optimized.json

The checker exhausts **all labelled loopless digraphs of orders zero through four**, with independently free ordered arcs: `1,1,4,64,4096` cases. For each it compares direct endpoint subset DP and ordered-cover convolution against the word-cut expression, checks complement reciprocity at `q=0,...,4`, checks the gate operator at `q=1,2,3`, and directly counts gate cycles modulo rotation. Gate strongness and comparison-defect conservation are checked too. Empty-core behavior is explicit. This is not merely a tournament census.

Further controls cover `D_N,U_N,T_N` through `N=10`, verify the explicit arithmetic path and minimum path-cover readouts, and compute four nested path responses for the order-three hostile. Literal and row-scaled operation identities are tested for `x=1,...,100`. Every check uses an explicit exception and survives `-O`; no repository modules or external libraries are imported.

The all-order statements follow from the proofs. The square-sum `3,10` input belongs to the independently enumerated companion packet; this checker does not establish eventual square-sum Hamiltonicity, Camion's tournament theorem, Rédei's theorem, or the tournament exclusions. The construction does not prove any global spectrum completeness claim. It identifies the exact comparison information a valid cross-domain construction must preserve or explicitly replace.
