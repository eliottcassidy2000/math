# An ambient rank-62 syllable bridge

**Exact finite linear algebra, with all-order consequences.** A *literal letter* of a span-three ribbon specifies the set of gaps 1, 2, 3 at which the new vertex beats an older vertex. A *syllable* is one of the 13 words listed in the [adaptive finite cone packet](../chanpak_cone_collision_placement_20261115/SUFFIX_CONES.md). These are different alphabets. The 8-dimensional cone counts and the 62-dimensional unrestricted count transfer are not interchangeable: the reduced rank can be one while the full rank is two. This note constructs a literal right cap that makes the distinction harmless on a long syllable run of an arbitrary ribbon.

## Exact bridge and the one-letter cap

Use the pinned minimal integral 62-dimensional transfer, with letter matrices $M_a$, positive physical path counts $hM_Wv$, and $U_j$ the full product for syllable $j$. Products act on column states from right to left in chronological order. Put $E=\sum_j \operatorname{im}U_j$. Take $Q$ to be the eight-row Hankel map given by the empty suffix and suffix syllables indexed $0,1,2,3,5,6,7$. These indices are zero-based in the cone packet order. Taking the same prefix basis as the cone, the eight reduced matrices $R_j$ computed by the invertible Hankel block are exactly those of that packet: the inputs to the Hankel formula are the same integer counts of literal words.

The replay verifies the following exact linear identities:

* $E$ has dimension 46, $Q(E)$ has dimension 8, and $QU_j=R_jQ$ on $E$ for every syllable. This restricted identity is not asserted on all 62 states.
* On the full space, $QU_j-R_jQ$ has image contained in $\operatorname{im}R_j$. Equivalently $QU_j=R_jQ'_j$ for a suitable $Q'_j$. Hence for every nonempty syllable word $W$, $\operatorname{rank}(QU_W)\le\operatorname{rank}(R_W)$ without any left padding: use this identity on its first unit, then intertwine on $E$ thereafter.
* The kernel of $QU_j$ is annihilated by $M_aU_j$ whenever $a$ has no distance-three backward edge. For last-syllable indices $1,4,6,7,8,9,11,12$ any of the eight letters works. Call an allowed $a$ a **cap**. Since $\ker(QU_W)\subseteq\ker(M_aU_W)$ when $a$ caps the last syllable of $W$, we obtain
  $$
  \operatorname{rank}(M_aU_W)\le\operatorname{rank}(QU_W)\le\operatorname{rank}(R_W).
  $$
  Every 13-unit syllable starts with a letter not using distance three, so one trailing syllable always suffices as well. Some padding really is necessary: $(0.3.1)(0.3.1)$ has reduced rank one but full rank two; appending literal letter $0$ gives full rank one.

The 38-dimensional invisible kernel on $E$ is especially short-lived on the syllable alphabet: the sum of its one-syllable images has dimension 12 and every next syllable kills that subspace. This does not say it dies under arbitrary continuation letters.

A separate exact check applies the same four tests (image in $E$, intertwining on $E$, first-unit residual, universal no-distance-three cap) to the 49 cyclic phases of the previously catalogued short trace-3 and trace-5 units. Precisely 23 pass: the original 13, the five omitted unrotated representatives $0.12$, $0.0.3$, $0.0.13$, $0.0.123$, $0.13.12$, and five additional phases $2.1$, $1.0.3$, $2.3.1$, $0.3.0$, $0.13.0$. They share the same state interface and padding theorem. This does **not** certify their growth cone; the existing positive certificate is for 13 only.

## Growth in arbitrary literal contexts

Fix $\rho=1.7699$. Let $B$ be a run of the original 13 syllables inside a literal span-three word $P B S$, with arbitrary literal prefix $P$ and suffix $S$. Write $|B|$ for its number of vertices.

* If $B$ has full rank at least two, use its final syllable as the cap for the preceding run. By the rank inequality the preceding run has reduced rank at least two. Then
  $$ H(P B S)\ge H(P)\,\rho^{|B|-15}/4000. $$
* If instead there is a qualifying cap immediately after $B$ and the capped factor has full rank at least two, the whole run has reduced rank at least two, and the bound improves to $H(P B S)\ge H(P)\rho^{|B|-12}/4000$.

Indeed the raw transfer from any literal prefix (not just a syllable prefix) is nonnegative, and after four retained syllables it lies in the cone of the relevant suffix tag. The comparison inequalities and expanding functionals therefore telescope from that actual prefix. Monotone Hamiltonian path insertion gives $H(P)$ as a lower bound for the readout after the first four retained syllables, cost at most 12 vertices. Further vertices, including the discarded terminal syllable in the first case, also cannot lower the count. If the run is too short to start the cone the asserted bound is trivial by insertion. This directly justifies use in arbitrary contexts without quotienting away the physical interface.

For a full-rank-at-least-two literal word, every contiguous factor has full rank at least two. With $r$ disjoint 13-syllable runs, covering all but $s$ vertices of a word of order $n$, applying the bound chronologically yields
$$ H\ge \rho^{n-s-15r}/4000^r. $$
In particular if $s+r=o(n)$ this retains growth rate at least $\rho$; the formula allows linear errors as well, at its explicit penalty. Qualifying caps save up to three powers of $\rho$ per run; capped factors may overlap at cap letters because the inequality only charges the run itself.

## When primeness certifies the full-rank condition

If a literal factor $B$ has rank one, its integer matrix is an outer product of integral column and row vectors (choose a primitive integral generator of its column space). Consequently $H(PBS)=\alpha(P)\beta(S)$ for integer functions. Choose their signs so that $\alpha(\varnothing),\beta(\varnothing)>0$; every $H(PB)$ and $H(BS)$ is positive, so both functions are positive on literal contexts.

For disjoint tournament pieces $C,D$, one can inject pairs of Hamiltonian paths into Hamiltonian paths of their union: keep the orders of both lists and at each step choose whichever foremost vertex beats the other. After removal the chosen vertex beats both remaining heads. Restriction recovers the original lists, so $H(C\cup D)\ge H(C)H(D)$ even if the cross edges are arbitrary. Applying this to $PB$ and $BS$ shows
$$ \alpha(P)\ge H(P)\alpha(\varnothing),\qquad
\beta(S)\ge H(S)\beta(\varnothing). $$
Thus a rank-one factor cannot occur between two nontransitive induced flanks of a prime-count word.

For a *strong* span-three ribbon every literal prefix of length at least seven is nontransitive: otherwise its source must be within the first four physical vertices, as every incoming-position predecessor in the prefix is beaten by it; this source then also beats every ambient later vertex by the span bound. The dual holds for suffixes. Therefore every contiguous factor with at least seven vertices on either side in a strong prime-count ribbon has full rank at least two. The segmented growth bound applies to runs in this bulk, charging the end vertices to $s$. It does not require declaring strong prime counts globally rank-two at the very ends.

### A quantitative test for hypothetical near-floor primes

Suppose a family of strong prime-count ribbons satisfies $H\le C5^{n/3}$ for a fixed $C$. In the seven-trimmed bulk, cover as much as possible by disjoint syllable runs, merging any adjacent ones. If $s$ vertices of the full word are uncovered and there are $r$ runs, then $r\le s+1$. The segmented inequality gives
$
s\log\rho+r(15\log\rho+\log4000)\ge n\left(\log\rho-\tfrac13\log5\right)-\log C .
$
Thus $s\ge 0.00197625\,n-O_C(1)$, even using the conservative constants without cap savings. Any bounded-distance-from-floor prime family must genuinely evade syllable parsing on a positive fraction of its letters, not just at a growing but sublinear number of seams. This is a necessary structural test, not an assertion that such prime families exist.

## Obstruction to a literal affine counter in this syllable grammar

For any fixed nonempty 13-syllable word $W$ and fixed literal ends $P,S$, the sequence $H(PW^kS)$ cannot be a nonconstant affine function of $k$ for all sufficiently large $k$. If every full power $U_W^k$ has rank at least two, the run-growth bound above gives an exponential lower bound. Otherwise a power has rank one: rank zero is excluded by the positive path count of every literal word. Its one-dimensional image remains invariant under $U_W$, so subsequent powers are exact scalar multiples, and the counts are eventually geometric. This is an obstruction to a literal exposed linear counter, not to normalized channels, nonperiodic programming, or an arbitrary-letter multiplexer.

There is also a cone-free obstruction on **every fixed-bandwidth literal alphabet**, not just the 13 syllables. A directed triangle in span $b$ has physical diameter at most $2b$: for ordered positions $i<j<k$, either the edge $k\to i$ is backward (and $k-i\le b$), or both $k\to j$ and $j\to i$ are backward. For a nonempty fixed period $W$, inspect the bi-infinite repetition locally. If it contains a triangle, translated copies by multiples of $|W|$ yield $\Omega(k)$ vertex-disjoint triangles inside $PW^kS$, for example use spacing $|W|(\lfloor 2b/|W|\rfloor+1)$. The shuffle injection gives at least $3^{\Omega(k)}$ paths. If the bulk has no triangle, every triangle in the finite word is within one of the fixed end windows (the corresponding end plus $2b$ vertices), because it must touch $P$ or $S$. A vertex of a nonsingleton strong tournament is on a triangle: otherwise every in-neighbor beats every out-neighbor, a strong-connectivity contradiction. Thus all nontrivial strong components lie in those end windows, and none can use both windows once they are far apart and their arcs point strictly forward. The graphs within each end window are independent of $k$ for large $k$ because the end phases are fixed. Their component counts multiply, with all middle vertices contributing one. Hence the path count is eventually **constant** in the triangle-free case. This finite local dichotomy excludes every nonconstant polynomial count from repeating a literal bounded-span motif (although it still does not exclude a normalized polynomial channel). The syllable cone above is stronger in its numeric growth constant, not necessary for the qualitative obstruction.

Nothing here proves a density bound on how many syllables arbitrary strong prime words admit. That is a separate covering or cone-extension problem. In particular the bound for arbitrary words is conditional on their certified runs, not the full span-three (1.7699) theorem.

## Exact obstruction to a single-syllable-period radix-1000 ladder

The number pattern $1001,1000001,1000000001$ equals $1000^k+1$ at $k=1,2,3$. But there do **not** exist a fixed nonempty motif $W$ made of the 13 units and fixed literal span-three ends $P,S$ for which $H(PW^kS)=1000^k+1$ at every sufficiently large $k$. This does not rule out a different host at each height, a non-syllable motif, or a larger width.

First, a rank-one power of $U_W$ would make the output eventually geometric, ruling out this sum of two exponentials. Thus every power has full rank at least two. The ambient cone lift forces $\rho^{|W|}\le1000$, so $|W|\le12$ because $17699^{13}>1000\cdot10000^{13}$. Second, the filtration $0\subset K\subset E\subset\mathbb Q^{62}$, where $K=E\cap\ker Q$, shows that the full characteristic polynomial of $U_W$ is $t^{54}$ times the characteristic polynomial of $R_W$. The ambient quotient by $E$ is killed by every unit, the restriction to $K$ is nilpotent (even two units kill it), and the induced map on $E/K$ is exactly $R_W$ by the projection. The scalar sequence $1000^k+1$ requires both $1000$ and $1$ as roots of this characteristic polynomial: by Cayley–Hamilton the scalar sequence obeys that recurrence, and the two independent terms force both roots, even from an eventual tail only.

There are exactly 99,420 nonempty weighted 13-unit words of cost at most 12 (each two-letter unit costs two, each three-letter unit three). The replay checks this entire bounded family, keeping physical order in the products. Over $\mathbb F_{101}$, of these 945 pass first at eigenvalue $1000\bmod101$; 21 also pass at $1$. Those 21 exact rational matrices have determinant pairs $\bigl(\det(1000I-R_W),\det(I-R_W)\bigr)$ in the two-element set
$$
\{(-4324134008000000000000000,-6464),\;(-3231798000000000000000000,-3030)\}.
$$
Neither determinant vanishes. Modular screening is sound since every entry denominator in the fixed rational generators is prime to $101$, which the script checks. This is a finite, exact *exclusion of a grammar*, not a heuristic extrapolation from the three given values.

## Replay

Run
```
python verify_syllable62.py path/to/span3_intrep.txt
```
with SymPy. The repository file is `04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt`; its SHA-256 is pinned in the script. The finite checks use exact rational linear algebra. The all-length claims follow from the identities and the independent cone proof, not a word enumeration.

For the radix replay additionally run:

```
g++ -std=c++17 -O3 radix_sieve.cpp -o /tmp/radix_sieve
python verify_radix.py path/to/span3_intrep.txt /tmp/radix_sieve
```
