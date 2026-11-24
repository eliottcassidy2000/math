# Apex gaps at logarithmic scale; the permanent exponent; the gamma connection

## 1. New: a linear-in-gap description bound

Let $H$ be the ordinary count of a tournament with a marked apex, $R_+=h+P$ the [exposed response](../chanpak_f1_physical_carry_20261114/ORDINAL_DEFECT.md), and $d=R_+-H$. Restrict to **at least two** apex outcontacts. For $D\ge2$ and $H\le X$, the number of distinct heights witnessed with $0<d\le D$ is at most
$$
 2^{128D}(1+\log_2(2X))\qquad(D\le X).
$$
The constant 128 is deliberately loose. In particular **almost every odd target has no such witness with $d\le(\log_2 X)/256$**. The fixed $c<4$ barrier in the prompt follows, but this rules out much more: all $o(\log X)$ gap budgets and a positive-constant logarithmic budget. No statement is made about one contact or about suitable larger gaps.

### Inherited common budget

The [fixed-gap normalization](../chanpak_f1_physical_carry_20261114/FIXED_GAP_RAYS.md) has finite pieces $G,C,B$ with total order at most $3D+1$ and one free transitive tail. The apex contacts lie in $C$. The two ordinal wings have heights $m,q$ (empty means one). With $p=\binom{\#\text{contacts}}{2}$, the nonnegative pair-defect identity gives $mq E_C\le d$. The least costly marked contact pair contributes at most $E_C/p$; its AND and OR contractions both contribute to that cost. The cheaper has height at most $E_C/(2p)$, so if $W$ deletes the two marked contacts in $C$,
$$
 H(G)H(W)H(B)\le D/(2p).
$$
In fact the joint backward arc budget in pivot orders for $G,W,B$ is $D/(4p)$, as in the [inherited sparse encoding](../chanpak_f1_physical_carry_20261114/SPARSE_GAP_ENCODING.md). It uses the odd-cycle packing identity to charge every triangle two paths and $\sum h_i-3\le\prod h_i-1$ for the three positive heights. A pivot recursively puts its in-neighbors before it and its out-neighbors after it; every backward cross-arc then gives a distinct directed triangle through the pivot. Recursion charges no triangle twice.

The inherited encoding inserts those arcs among $O(D^2)$ positions, giving $2^{(D/(4p))\log_2 D+O(D)}$ finite templates. Even that bound covers its stated endpoint $c=4p$, by $\log D=\log\log X-\log\log\log X+O(1)$ at the endpoint. But we can remove the whole $\log D$ factor in the description cost.

### New tradeoff: a large triangle packing leaves fewer arcs to encode

Within any of the three blocks $T$ fix an inclusion-maximal family of $k$ vertex-disjoint directed triangles. Write $t$ for the total number of directed triangles and $h=H(T)$. The complement of the selected vertices is transitive. In the [Grinberg–Stanley packing expansion](https://arxiv.org/abs/2307.05569), the selected triangles alone contribute $3^k$. Every other triangle meets at most three of the selected ones. Allow it together with any subset of the remaining selected ones. Packs so produced contain exactly one nonselected triangle and are disjoint across its choices. Therefore
$$
 h\ge 3^k+2(t-k)3^{\max(k-3,0)}.
$$
In particular $t\le15h\,2^{-k}$: use $h\ge3^k$, the bound $t\le k+(27/2)h3^{-k}$, and $k(2/3)^k\le1$ for $k\ge1$. If $k=0$ the block is transitive.

This provides **two compatible encodings in the same pivot order**:
1. Mark the $3k$ packed vertices. Among the unmarked vertices record only their true transitive order relative to the pivot order. Its inversion count is exactly the number of backward arcs among them. A permutation on $u$ positions with at most $b$ inversions has an inversion vector of length $u$ summing to at most $b$; there are at most $\binom{u+b}{b}$ candidates. This costs $O(D)$ bits even for all three blocks together. It is not necessary to find a pivot order whose unmarked subsequence is already sorted.
2. Only pairs touching a marked vertex still require an arc bit. For one block there are at most $12Dk$ such pairs since its order is at most $4D$. At most $15D2^{-k}$ are backward by the tradeoff and the pivot bound. For parameters $M,b$, the number of subsets of at most $b$ out of $M$ is bounded by $\binom{M+b}{b}$ (pad with one repeated dummy) and has log cost at most $b\log_2(e(1+M/b))$. This bound increases with $b$. At $b=15D2^{-k}$, it is at most $15D2^{-k}(2k+2)\le30D$ bits, using $1+(4/5)k2^k\le2^{2k}$ for $k\ge1$. The zero-packing case has no such choices.

To account for the constant, the three marked-vertex subsets cost at most $4D$ bits; the joint transitive-complement inversion vector costs at most $4D+D/4$ bits; the three remaining arc sets cost at most $90D$ bits. Block sizes take at most $6D$ bits, since $(4D+1)^3\le 2^{6D}$. The apex contact mask, two named contacts' incident arcs within the core, and suffix flag take at most $4D+9D+D$ bits by direct overcounting; the two contacts may be given canonical names and inserted after assigning the labels of $W$ since labels do not matter. All cross-wing orientations are prescribed. Total is less than $128D$. This encodes a finite overcount, not an assertion every code is valid. Each finite template supplies at most $1+\log_2(2X)$ heights up to $X$, by the same exact dyadic tail inverse, proving the bound.

**Algorithmic consequence.** Enumerate the marked positions, inversion vectors and sparse incident arc sets instead of all orientations, with the above caps. Recover a literal template, recompute its exact profiles by subset/endpoint DP on $O(D)$ vertices, reject the wrong defects, and test the remaining exponent by integer division and a power-of-two check. This is an exact $2^{O(D)}\operatorname{poly}(\log X)$ search for a given target and gap budget, not a heuristic or an unrestricted lower bound. **replay.py** checks the packing bound, the pivot budget and actual reconstruction through all labeled tournaments of order six.

## 2. The permanent decision result and the source of $4/3$

This section uses the pinned [arithmetic border lower bound](https://hexagonmath.org/pdf/2611.00954v1) and [Bürgisser's corrected graph factor theorem](https://arxiv.org/html/1812.06828v1) as inputs, as in [PR #249](https://github.com/1337isnot1337/solved-conjectures/pull/249). It does not re-prove or claim to have rerun the formalization.

Let $g=y-\operatorname{per}_m(X)$. The ideal of its graph in complex affine space is $(g)$, as division by the monic linear polynomial in $y$ shows. Follow a generic accepting path of an exact complex algebraic computation tree: at each branch keep a nonempty Zariski-open subset of the irreducible graph, and represent every register by a numerator/denominator pair. Each arithmetic gate needs at most four division-free operations for its pair. Some equality test must have a nonzero ambient numerator vanishing on the graph. Otherwise the accepting path would impose only nonvanishing conditions in ambient space (and tautological equations); there is then an off-graph point on the same path, avoiding all poles. This yields a nonzero multiple of $g$ with straight-line cost at most $4T+3$, independent of multiplicity. Notably one extracts **one** equality witness, not all predicates of the entire tree.

Bürgisser's graph-factor estimate gives a border program for $g$ with size $O(M(m^2)(T+1))$. Specialize $y=0$; converting the operation count to all-entry circuit size costs a constant factor (count used inputs and constants). With $M(m^2)=O(m^2\log m)$ and source border bound $2^{a(\log_2 m)^{4/3}}$, the quantitative result is
$$
 T+1\ \ge\ \frac{2^{a(\log_2 m)^{4/3}}}{O(m^2\log m)}.
$$
The polynomial loss cannot change the exponent $4/3$. This is for exact complex $+,-,\times,\div$ and equality tests, arbitrary constants and defined divisions. It says nothing directly about bit-complexity on Boolean inputs.

**Where $4/3$ actually enters.** The exponent was already in the source theorem, not produced by PR #249. In its flat case (Proposition E.10 and Theorem E.11), specifying the balanced cut uses $k$ bits and the bound for the scaled end-set functional is of order $\sqrt{k}$. The penalty is therefore $k^{3/2}$ against a geometric saving of order $\ell k$; take the least useful $\ell\asymp\sqrt{k}$. The hub separation is order $k^{-1/2}$, so the saving forces $\log(\text{circuit size})\gtrsim\ell k^{3/2}\asymp k^2$. Meanwhile the constructed input has $\log m\asymp k\ell\asymp k^{3/2}$. Eliminating $k$ gives $4/3$. The contracting case supplies a stronger bound. This is a scale comparison with a preserved observable and a cost, not a Bott period.

**A literal transfer back to an existing tournament chart.** The [split-vertex weighted construction](../chanpak_tournament_sl2_20261124/WEIGHTED-COMPLEXITY.md) is an exact projection of the Hamilton-cycle polynomial into the independently weighted Hamilton-path polynomial of one fixed strong $2m$-vertex tournament. Apply the same generic-path and graph-factor argument to the Hamilton-cycle graph, whose degree is $m$, using the pinned source's Hamilton-cycle border bound. Then restrict any proposed weighted-tournament decision tree to the split-vertex projection: each arc weight is an original variable, zero, or one. Thus deciding the weighted Hamilton-path graph of that fixed host also costs $2^{\Omega((\log m)^{4/3})}$ in the same complex model. This does **not** transfer to unweighted counts or to the sparse-apex range. The two proofs do share a method: find a small retained certificate before charging cost, rather than counting everything in the original presentation.

## 3. Bott stabilization versus recursive scale

The repo's [marked Bott session](../chanpak_bott_cayley_periods_20261115/BOTT.md) distinguishes stable real period eight (and the matrix-algebra step $\operatorname{Cl}_{0,n+8}\cong\operatorname{Cl}_{0,n}\otimes M_{16}(\mathbb R)$ in its convention) from recurrence of an integer height. A concrete preservation test: for marked transfer matrices $M_i$ and endpoints $u,v$,
$$
 (u\otimes e_1)^{\mathsf T}\prod_i(M_i\otimes I_{16})(v\otimes e_1)
 =u^{\mathsf T}\prod_iM_i\,v .
$$
The analogous unmarked trace instead multiplies by 16; modulo 8 it loses every residue. Stabilization is useful only with the right mark. Likewise the apex encoder retains the two cheapest contacts and the inversion data, and the permanent extraction retains the graph-generic equality numerator. There is a useful way to compare the numbers without identifying the operations: in the permanent proof take $k=r^2$ and $\ell=100r$. Doubling $r$ multiplies the cubic input-log budget by 8 and the quartic lower-bound budget by 16, so the exponent is $\log_8 16=4/3$. The Bott statement also contains 8 and 16, but its move is additive index shift by eight and a tensor with a matrix algebra, not this doubling or a relation between the two costs. No stable eight-step identification produces the $4/3$ scale or an absolute count inverse.

## 4. Where $8\pi$ comes from, and what the ball formula really shares

[Schoenfeld (1976), §6 and Corollary 1](https://www.ams.org/mcom/1976-30-134/S0025-5718-1976-0457374-X/S0025-5718-1976-0457374-X.pdf) gives the explicit two-sided $|\pi(x)-\operatorname{li}(x)|<\sqrt x\log x/(8\pi)$ bound for $x\ge2657$ under RH. The **leading-constant mechanism** is visible before the explicit finite-range estimates. Count positive zero ordinates by the Riemann–von Mangoldt formula with leading density $\frac1{2\pi}\log(t/(2\pi))\,dt$. Under RH the weight $1/|\rho|$ is asymptotic to $1/t$. Doubling for conjugate zeros and integrating $t^{-1}\log t$ gives
$$
 \sum_{0<|\gamma|\le T}\frac1{|\rho|}
 =\frac{1}{2\pi}\log^2 T+O(\log T).
$$
Indeed the conjugate factor 2 cancels the integration factor $1/2$. At Schoenfeld's smoothing scale $\log T=\tfrac12\log x+O(\log\log x)$ (see (6.6)–(6.9)), squaring that half contributes $1/4$, hence $1/(8\pi)$. Partial summation from the Chebyshev bound to $\pi(x)$ removes one logarithm without changing this leading coefficient. The exact inequality and cutoff additionally need the tail and lower-order bounds and finite checks in the paper. The leading-constant calculation does not establish sharpness.

The [prior volume session](../session127_gsma_norm_collapse_20261031/volume/THEORY.md) supplies a **real common normalization**, not a ball-packing explanation of the eight. For $V_s=\pi^{s/2}/\Gamma(s/2+1)$, analytically continuing the unit-ball formula, one has
$$
 V_s\,\pi^{-s/2}\Gamma(s/2)=2/s,\qquad
 \xi(s)=\frac{(s-1)\zeta(s)}{V_s}.
$$
The archimedean factor of the completed zeta function is thus exactly reciprocal to the ball-volume expression up to $2/s$. Both are normalized by the Gaussian; the zero density is governed by the **phase at complex argument**, not by the volume of a real eight-ball. Indeed $V_{s+8}/V_s=(2\pi)^4/((s+2)(s+4)(s+6)(s+8))$, a shift law, not a periodic invariant. And $V_3=4\pi/3$ is a geometric volume; its $4/3$ is not the logarithmic-complexity exponent.
