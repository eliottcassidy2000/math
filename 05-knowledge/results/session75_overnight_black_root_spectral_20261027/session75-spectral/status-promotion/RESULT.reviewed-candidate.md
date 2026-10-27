# Spectral budgets for bounded tournament substitution grammars

Session75; 2026-10-27. **PROVED self-contained in the stated interfaces;
FINITE-EXACT benchmarks; pending independent parent audit.** No historical
priority claim. This packet neither closes native all-word reachability nor
proves a ROOT positivity statement.

## Results first

Write `K(T)uv=1` if `u→v`, `−1` if `v→u`, and `0` on the diagonal; `||·||*`
is the sum of singular values, not the trace of a skew matrix. Thus K is
`Adj−Adjᵀ`, with no factor 1/2; switching the sign convention changes none
of the norms or budgets. Let a baseline
be an arbitrary-depth, disjoint substitution tree whose quotient arities are
at most `b≥2`. Transitive leaves and unbounded ordered joins are allowed,
because they compile into binary joins. Put

    α=1/(b+1), h₂(x)=−x log₂x−(1−x)log₂(1−x),
    C_b=2√(α(1−α))/h₂(α).

1. **All shapes, including combs:** `||K(B)||*≤C_b n log₂n`.
2. If `P` has `K(P)K(P)ᵀ=nI−J` (in particular a correctly constructed
   Paley tournament), every such baseline differs on at least

       m ≥ max(0, ceil([n(n−1)−C_b n^(3/2)log₂n]/4))                 (E)

   unordered named pairs. For fixed b, this is `(1/2−o(1))` of all arcs. If the grammar
   contains transitive tournaments, the upper bound `m≤n(n−1)/4` for some
   baseline makes the leading constant optimal.
3. For the join/C3 grammar, `b=3`, `C_3≈1.0674827506<13/12`.
   Every flat Paley order `n≥512` therefore requires at least
   `n(n−1)/8` arc exceptions. This is an explicit sufficient threshold,
   not a claimed least threshold.
4. **Different currencies:** a final perturbation of rank `r` must have
   `r≥n−1−C_b√n log₂n`; thus `k` additive rank-two skew corrections require
   `k≥ceil((n−1−C_b√n log₂n)/2)`. Arbitrary block reversals or arbitrary
   module-occurrence exceptions need not be rank two.
5. Fixed-size legal path complements, allowing arbitrary legal path
   relocations and free global converse, require

       t ≥ ceil( max(0,n(n−1)−C_b n^(3/2)log₂n)/(4(n−1)) ).         (C)

   For fixed b, the bound is `n/4−O(√n log n)`, stronger than the separate nuclear
   Lipschitz bound `Ω(√n)`. Size changes or clone expansion need their
   explicit propagated arc tariffs; they are not silently covered.
6. **Finite exact distinction:** distance from Paley7 to the whole
   arbitrary-depth join/C3 grammar is **3** arc flips; for Paley11 it is
   **11**. In particular Paley11 needs at least two fixed-size complements from
   this grammar, by the ten-arc tariff; attainment is not asserted. These
   are not the inherited native-complement distances
   **6 to forest frames / 15 to transitive tournaments** at order eleven.
7. **Regular-to-regular control:** the C3 tower `R_k` on `n=3^k` has trace
   norm `(2k/√3)n`. For odd `k`, Paley over the *field* `F_(3^k)` versus
   any relabeling of this tower needs at least

       n[(n−1)−2k·3^((k−1)/2)]/4                                 (T)

   arc flips. At `n=27` this is `54`, hence at least `3` fixed-size
   complements. At `n=243` it is `9234`, hence at least `39` complements.
   Both endpoints are regular, so degree/score variance cannot see this
   separation. These are lower bounds, not attained-distance claims.

## 1. Inheritance and exact interface

The byte-pinned sources are in `SOURCE_PINS.json`.

* `paid_root_cycle_budget_20261026/README.md`, T2, already defines a **named**
  rooted ordered-join/C3 tree with transitive blocks, a marked Hamilton path,
  and a named undirected exception set `D`. Arc queries use the LCA quotient
  orientation XOR membership in `D`. It already proves
  `(B,D,P)→(B^op,D△E(P),P)` under legal path complement. The source's
  `O(n log n)` is **combinatorial bits**, not a proved trace-norm bound.
* Session73 `session73-ports/REPORT.md`, G1, is the inherited exact weighted
  feedback-arc substitution identity. Its positive pair-factor vertex costs
  are not native source masses. G2 retains actual port endpoints and parallel
  component-incidence edges; repeated toggles of one physical pair combine
  by XOR. It does not make arbitrary shared-port gluing a substitution.
* `THM-1960` plus full slug in the manifest gives regular-child spectral
  splitting, with the explicit irregular-child failure. We reprove the
  special tower identity; we do **not** apply that splitting to arbitrary
  unbalanced or irregular parse nodes.
* The composition-classification packet defines disjoint modules and the
  decorated classifier. We need only the elementary module predicate and
  the supplied parse, not any classification or finite prime census.
* Session74's exact Paley11 native distances stay in their own currency.

The user-supplied trace-norm separation is recovered here as a self-contained
all-shape theorem. We have not located an earlier repository proof of that
precise asymptotic; we do not attribute this proof to the combinatorial-bit
sentence. Original finite-bank templates may be a strict subset of all
arity-at-most-b tournaments; the proof deliberately works in the larger class.

### Formal parse

Each leaf is one named vertex. An internal node has `2≤q≤b` nonempty,
pairwise disjoint child sets and a tournament quotient on its `q` positions.
Every cross pair of child vertices receives that quotient orientation.
A transitive leaf of integer size `w` expands into `w` named leaves via
binary joins. An ordered join of any arity is likewise associatively
binarized. A finite seed bank of maximum size `b` is covered after including
binary joins and converses if needed. Unary nodes are removed. Reusing one template description for several
disjoint copies (a compressed description DAG) is harmless: expand its
occurrences to this tree. Identifying their physical vertices is different
and is not justified merely by that shared description.

Every subtree is a module of the resulting baseline: an outside vertex has
one uniform direction to all its vertices. The induced tournament on any
nonempty vertex subset has a parse of maximum arity `b`: intersect every
child with the subset, discard empty children, restrict the quotient, and
suppress unary nodes. A restricted quotient need not still belong to the
original bank; it belongs to the enlarged arity-b class used for the bound.

Ports are **references to these physical vertices**, possibly repeated.
Their names do not create extra vertices or erase pair incidences. Any
operation that actually identifies vertices from different child sets must
first establish a new lawful disjoint parse; the theorem cannot be applied
to overlapping syntax as though it were a disjoint substitution.

## 2. A balanced module exists even in a maximally unbalanced parse

**Lemma.** If `n≥2`, the supplied arity-b parse has a proper module of size

    n/(b+1) < m ≤ bn/(b+1).                                     (2.1)

Start at the root and descend while some child has size greater than
`bn/(b+1)`. The final node has size greater than this threshold, and every
child has size at most the threshold. Its largest child has size greater
than `n/(b+1)`, because there are at most `b` children. This child is the
required proper module. Singleton children are allowed. The bound gives
both sides of the cut size at most `bn/(b+1)`.

For such a module, reorder vertices with it first. The cross skew matrix is

    X = [ 0       1_m sᵀ ]
        [ −s 1_mᵀ   0   ],   s∈{±1}^(n−m).

Its two nonzero singular values are `√(m(n−m))`, so

    ||K(B)||* ≤ ||K(B[M])||* + ||K(B[V\M])||* +2√(m(n−m)).        (2.2)

The complement need not be a module. Its **induced parse** is what enables
the next induction. This is the missing step in an original-depth argument.

### Entropy induction and the constant

For `0<x≤1/2`, the function `g(x)=2√(x(1−x))/h₂(x)` is decreasing.
Its derivative has the required nonpositive sign exactly when

    (1−x)ln(1−x) − x ln x ≥0.

The expression on the left is concave on `[0,1/2]`, has value zero at both
ends (using `0 ln0=0`), and hence is nonnegative. Thus for
`x∈[α,1−α]`, `2√(x(1−x))≤C_b h₂(x)`.

Induct using (2.2), the induced-parse property, and

    n log₂n − m log₂m −(n−m)log₂(n−m) = n h₂(m/n).

The singleton case is zero. This proves claim 1 at **arbitrary original
parse depth**, without a regularity or equal-leaf-weight assumption.

For `b=3`, `h₂(1/4)=2−(3/4)log₂3`. Since `3^5<2^8`,
`h₂(1/4)>4/5`; hence `C_3<5√3/8<13/12`, the last inequality following
by squaring positive quantities (`2700<2704`). These rational integer
comparisons, not the displayed decimal, certify our simple explicit bound.

### Hostile original-depth sum

In the C3 comb `S_j=C3[S_(j−1),1,1]`, `|S_j|=2j+1`, the exact cross-node
trace contribution is `2√(4j−1)`. The triangle sum over the **original**
nodes is therefore `Θ(n^(3/2))`, not `O(n log n)`. It is merely a loose
upper bound; it does not disprove the theorem. The new balanced-module
induction is essential. The exact checker includes a comb for which this
naive sum's integer lower enclosure exceeds the rigorous `13n ceil(log₂n)/12`
upper bound on the actual matrix norm.

## 3. The trace dual certificate gives the sharp leading edge density

Let `A=K(P)`, `B=K(Baseline)`, and let `m` count unordered disagreeing pairs.
With Frobenius inner product `⟨X,Y⟩=tr(XᵀY)`, every agreeing pair contributes
`2` and every disagreeing pair contributes `−2`, so exactly

    ⟨A,B⟩=n(n−1)−4m.                                            (3.1)

If `AAᵀ=nI−J`, then `||A||op=√n` and `||A||*=(n−1)√n`.
Nuclear/operator duality gives

    n(n−1)−4m ≤ √n ||B||* ≤ √n U,
    m ≥ [n(n−1)−√n U]/4.                                       (3.2)

The first inequality remains correct when its left side is negative.
Taking nonnegative integer ceilings proves (E). It holds uniformly under
all relabelings, all legal parse shapes, all bank template choices, and all
port decorations that do not change the decoded tournament.

For comparison, triangle/Frobenius alone gives

    Δ=(||A||*−U)_+,  ||A−B||F²=8m,
    m ≥ Δ/4,                  (one arc has nuclear norm 4),
    m ≥ Δ²/(8r),              r=rank(A−B)≤n−1 for odd n.         (3.3)

The second is quadratic rather than the first's `n^(3/2)` scale, but
asymptotically gives only `n²/8`. The flat-target dual certificate (3.2)
improves the leading lower bound to `n²/4`.

For any tournament, a uniform random total order has expected backward
count `n(n−1)/4`; some order does no worse. Transitive baselines in the
bank therefore supply the matching leading upper bound.

### Explicit threshold and cycle currency

For `n≥512`,

    (13/12)√n log₂n ≤ (n−1)/2.

Indeed `(n−1)/(√n log₂n)` is increasing for `n≥512` (differentiate its
logarithm), and at `512` the squared inequality is the integer comparison
`169·512·81≤36·511²`. Substitution in (3.2) gives `m≥n(n−1)/8`.
For a named exception graph on all n vertices,
`β(D)=m−n+c(D)≥m−n+1`; hence (E) also gives a cycle-debt lower bound
with `n−1` subtracted. Edge, cycle, block, and occurrence debts remain distinct.

## 4. Flat targets include finite-field Paley, proved here

Let F be a finite field of odd cardinality `q≡3 (mod4)` and χ its quadratic
character with χ(0)=0. Define `A_xy=χ(y−x)`. Since χ(−1)=−1, this is a
±1 skew tournament matrix. Its row sums vanish. For `a≠b`,

    Σ_x χ(x−a)χ(x−b)=−1.                                       (4.1)

To see this, rescale/translate to `Σ_t χ(t(t−1))`. The number of pairs
`(t,z)` satisfying `z²=t(t−1)` equals `q+Σ_tχ(t(t−1))`. Completing squares
identifies this with `(2t−1−2z)(2t−1+2z)=1`; each nonzero first factor gives
exactly one pair, so the number is `q−1`. This proves (4.1), and therefore
`AAᵀ=qI−J`.

For prime q this is the usual square-residue construction modulo q. For
nonprime prime powers it is a **finite-field** construction, not arithmetic
modulo the composite integer. The exact control at q=27 uses
`F_3[x]/(x³−x−1)`; the cubic has no root in F3 and is irreducible.
We do not transfer a prime-only inherited theorem to composite residues.

## 5. Rank budget, blocks, occurrence counts, and clones

Let `E=A−B` have rank r. On `ker E ∩ 1^⊥`, a subspace of dimension at least
`n−r−1`, `Bv=Av` and `||Av||=√n||v||`. By the singular-value min-max
principle, B has at least `n−r−1` singular values at least `√n`. Thus

    ||B||* ≥ (n−r−1)_+√n,
    r ≥ n−1−U/√n.                                             (5.1)

Alternatively restrict the domain to an orthonormal basis of that subspace.
The resulting rectangular matrix has all those singular values equal to
√n, and restriction by an isometry cannot increase the nuclear norm. This also proves the claim without needing an
interlacing formula for an arbitrary perturbation.

A final **additive** decomposition `E=Σ_(j≤k)E_j`, each skew and rank≤2,
has rank≤2k. Formula (5.1) proves the stated block-rank budget. Shared
endpoints, overlapping supports, or cancellations do not invalidate rank
subadditivity, provided those actual summands are retained.

* Reversing one individual arc has rank2, nuclear norm4, Frobenius²8.
* Reversing a **uniformly oriented** complete bipartite interface A×B
  changes the matrix by `±2(1_A1_Bᵀ−1_B1_Aᵀ)`, rank2, nuclear norm
  `4√(|A||B|)`, but it reverses `|A||B|` individual arcs.
* Reversing an arbitrary nonuniform rectangle is not automatically rank2.
  With overlapping sequential block reversals, each actual intermediate
  rectangle must be checked, rather than borrowing its first orientation.
* Replacing an arbitrary size-n module by an arbitrary tournament can make
  Paley in one occurrence. There is no unqualified occurrence-count theorem.
  A per-occurrence rank cap R instead gives the valid bound
  `kR≥n−1−U/√n`. Bounded-arity quotient changes have controlled quotient
  rank; arbitrary internal replacements do not.

For lawful substitution `Q[T_1,…,T_q]`, write `w_i=|T_i|`,
`H_vi=1/√w_i` on block i and zero elsewhere. Then exactly

    K=diag(K(T_i))+H [diag(√w_i) K(Q) diag(√w_i)] Hᵀ.            (5.2)

H is an isometry, so the cross part's nonzero singular values are those of
the displayed **weighted** quotient. Dropping the internal blocks loses
spectral mass; dropping the weights loses multiplicity. No regular-child
spectral splitting is claimed unless all child mean vectors are kernels.
A compressed singleton with weight n and zero internal matrix does not
represent an arbitrary n-vertex tournament. Clone mergers are lawful only
with their full internal parse and external uniformity, or with all lost
information charged as exceptions. The checker tests the Frobenius losses
and verifies that Paley7 has no proper nontrivial module.

## 6. Native complement clock and a regular-to-regular example

For a legal directed Hamilton path P, let its signed path matrix be L
(nonzero only on path-neighbor pairs, with the actual direction). Then

    K(tau_P T)=−K(T)+2L.                                       (6.1)

Consequently

    | ||K(tau_PT)||*−||K(T)||* | ≤2||L||*≤4(n−1).               (6.2)

Up to a unitary diagonal gauge, iL is the undirected path adjacency matrix,
so its eigenvalues are `2cos(jπ/(n+1))`. Thus, for odd n,

    ||L||*=2cot(π/(2(n+1)))−2,
    tariff=2||L||*=4cot(π/(2(n+1)))−4 ~8n/π.                   (6.3)

For even n replace cot by csc. Formula (6.2) independently yields a nuclear
clock lower bound `( (n−1)√n−U )_+ /(2||L||*)`. It is weaker than the arc
clock bound here, but applies to any starting class with a trace cap.

For the stronger bound, the inherited exact decoder update changes D by
XOR with the n−1 physical path edges, and converses its baseline. After t
complements, with arbitrary legal relocations and free global converses,
the baseline is still B or B^op and `|D|≤t(n−1)`. If free relabellings
are also allowed, transport the baseline, path and D jointly; the baseline
is B or B^op up to that accumulated relabelling. Add converses to the bank
without increasing b, then apply (3.2). This proves (C), **without needing
the intermediate tournaments to remain modular**. Module shattering is
allowed and paid by D. With a nonzero initial exception count m0 replace
t(n−1) by m0+t(n−1).

If allowed operations change size, use the exact number of final named
pairs affected by each propagated update. In particular cloning a vertex
can multiply the number of old exceptions; charging merely one old record
is invalid. No arbitrary native-word realization theorem is supplied.

The bounded-bank hypothesis is essential: if b grows to n, a bank containing
the target as one quotient represents it with zero exceptions. The finite
formulas retain their b-dependence, but their fixed-b asymptotic conclusions
are not uniform over growing banks.

### C3 towers

Let R1=C3 and R_k=C3[R_(k−1),R_(k−1),R_(k−1)]. Every child is regular.
The subspaces of blockwise zero-sum vectors and block means are orthogonal
and invariant. The former gives three copies of the old nonzero spectrum;
the latter gives `3^(k−1)` times the C3 spectrum. Thus

    N_k=3N_(k−1)+2√3·3^(k−1), N_1=2√3,
    N_k=2k·3^(k−1)√3=(2k/√3)n.                                (6.4)

This is the inherited regular-child mechanism, with no irregular extension.
For odd k the field Paley target has the same order n and flat norm.
Substituting (6.4) into (3.2) proves (T). Both score vectors are constantly
`(n−1)/2`. Formula (T) is uniform over relabelings and converse of the tower.

## 7. Optional positive vertex weights (not ROOT source-law replacement)

For positive weights a_v, put W=Σa_v, D=diag(√a_v), and K_a=DKD.
The same balanced **vertex-count** cuts satisfy

    ||K_a(B)||* ≤ W ceil(log_((b+1)/b) n).                       (7.1)

At a cut the cross norm is `2√(W_AW_B)≤W_A+W_B`; recursively each original
vertex weight is charged at most the displayed count-balanced depth. Leaf
matrices are zero. This is a safe weighted theorem, not an entropy-of-mass
bound: a two-leaf highly unequal weight split already defeats the latter
naive induction.

For m_a=Σ_(disagreeing {u,v})a_ua_v and a flat unweighted target A,

    ⟨DAD,DBD⟩=W²−Σa_v²−4m_a,
    m_a ≥ [W²−Σa_v²−||DAD||op ||DBD||*]/4,                     (7.2)
    ||DAD||op≤(max a_v)√n.

One may retain the actual weighted operator norm for a stronger certificate.
There is no claim that unequal native source masses may be replaced by
uniform vertices or fresh probabilities, nor that spectral debt pays a
signed ROOT merging current.

## 8. Finite-exact distance certificate and reproduction

`certificate.json` stores complete subset costs and one optimal labelled
parse for P3, P7 and P11. Every subtree in any candidate parse induces a
partition of its named set into two or three nonempty parts. A binary root
chooses the cheaper orientation of its interface. A cyclic root chooses
the cheaper of the two C3 orientations. Transitive ternary roots already
compile to binary roots. Thus the Bellman recurrence is exhaustive for
all finite join/C3 parses, not just balanced trees or fixed-depth examples.

The checker uses a separately implemented restricted-growth partition
enumerator and direct arc counts to verify **every Bellman inequality**
against the supplied costs. Together with the zero singleton cases, these
inequalities inductively lower-bound every possible parse; the explicit
parse attains each asserted final distance. A separate direct all-pair
comparison verifies the witness. Paley7's lower bound is additionally
checked by enumerating all 232 masks within two flips and testing exact
recursive grammar membership. No floating-point eigenvalue test bears any
proof claim. The source pins, finite-field Gram identity, arithmetic constant
and threshold comparisons, XOR complement tariff, shared-port/clone hostiles,
and unbalanced-tree cut property have explicit RuntimeError guards.

Reproduce from this directory:

    python3 check.py > receipt.normal.json
    python3 -O check.py > receipt.optimized.json
    cmp receipt.normal.json receipt.optimized.json
    sha256sum -c SHA256SUMS

The finite universes and their counts are printed in the receipt. The
all-order conclusions come from the proofs above, not an extrapolated scan.
