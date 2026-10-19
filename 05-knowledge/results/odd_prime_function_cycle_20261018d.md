# Odd-order arithmetic, signed fibres, and odd-cover descent

Status: **PROVED** elementary statements with **FINITE-EXACT** independent controls. This is a source-bound interface for fixed ternary labels and binary shades, not a new bound on the chromatic number of the plane. No priority claim is made for classical graph-cover, signed-incidence, or quadratic-character facts.

The useful conclusion is stronger than an analogy: a loopless matching thickening with an odd number of vertices in each source fibre has **exactly the same fixed-label, uniformly pinned six-colour support** as its source graph. Majority gives a constructive descent. Oddness does not preserve unrestricted colourability, extension counts, or a metric realization.

## 1. Inheritance and the live coordinates

The closest proved mechanisms are:

* [LEM-004, tournaments are odd functions](../../01-canon/theorems/LEM-004-tournaments-are-odd-functions.md): a circulant tournament on an odd cyclic group uses an antisymmetric sign on nonzero **differences**; multiplicative Legendre signs are a distinct construction.
* [THM-4644, Fermat primes, rigid circulants, and Wieferich orbit explosion](../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md): prime-power residue actions must retain unit versus nonunit strata. This note does not reuse a prime-field calculation as a prime-power calculation.
* [Plane six reframes, §2–3](plane_six_reframes_parity_20261011b.md): after fixing a three-way colour label, the remaining two-way choice is a signed graph problem, with actual cycle and pin-path obstructions. Its classical balance attribution remains there.
* [Reanchoring kernel, §7](reanchoring_kernel_20261018c.md): in the current CRT chart, a proposed ternary label must agree with every prescribed colour modulo 3 **before** testing binary parity. The isolated pinned vertex is the minimal hostile. A fixed cycle basis also misses label-dependent monochromatic cycles.

The concept board is: translation orbits; involution pairs; valuation layers; signed edge equations; fixed-label covers; retained boundary pins. The canonical losses are respectively an odd orbit, a fixed point, a nonunit residue, a modular Gram matrix, a nonconstant label fibre, and a pin imposed on only one lift.

Here the six colours are explicitly

\[
c_v=3b_v+4\ell_v\pmod6,\qquad b_v\in\{0,1\},\quad\ell_v\in\{0,1,2\}.
\tag{1}
\]

Thus $b_v=c_v\bmod2$ and $\ell_v=c_v\bmod3$. On edges with different labels, the colours are automatically different. On the **actual** subgraph $G_\ell$ of equal-label edges, one needs $b_u\ne b_v$. This chart is fixed throughout; it must not be silently identified with the different chart $c=2\ell+b$.

## 2. Translation oddness is not involution oddness

Let values lie in $\mathbb Z/q\mathbb Z$, $q\ge2$. Translation by $a$ on $\mathbb Z/m\mathbb Z$ has $g=\gcd(m,a)$ orbits of common length $L=m/g$. Along an orbit the equation

\[
f(x+a)=-f(x)
\tag{2}
\]

determines all values from one value $t$. Closing the orbit imposes no restriction if $L$ is even, and imposes $2t=0$ if $L$ is odd. Consequently the exact number of solutions is

\[
\begin{cases}q^g&L\text{ even},\\ \gcd(q,2)^g&L\text{ odd}.\end{cases}
\tag{3}
\]

In particular, odd $m$ and odd $q$ force the zero function. The same argument over any field of characteristic different from 2 gives zero on every odd translation orbit. An odd-order additive group has no nontrivial homomorphism to $\{\pm1\}$, since the image of any element has order dividing both an odd integer and 2.

In contrast, the involution $x\mapsto-x$ has $f=\gcd(m,2)$ fixed points and $(m-f)/2$ two-element orbits. Its anti-invariant functions number

\[
\gcd(q,2)^f q^{(m-f)/2}.
\tag{4}
\]

For odd $m$ and characteristic different from 2 there is one forced zero at 0 and $(m-1)/2$ free pairs. Restricting nonzero residues to signs gives the $2^{(m-1)/2}$ intrinsic circulant orientations from LEM-004. These statements involve different permutations; the adjective “odd” alone supplies no map between them. Characteristic 2 erases the distinction between a sign and its negative.

## 3. Prime powers: the missing valuation coordinate

Let $p$ be an odd prime, $r\ge1$, and $U_r=(\mathbb Z/p^r\mathbb Z)^\times$. Every sign character $U_r\to\{\pm1\}$ factors through reduction modulo $p$: its restriction to the kernel, of odd order $p^{r-1}$, is trivial. There are exactly two such characters, the trivial one and the quadratic character $\chi_p$ pulled back from $\mathbb F_p^\times$.

An elementary square-lift proof also retains the branch. Given $x^2\equiv u\pmod{p^j}$ with $p\nmid x$, the unique next digit satisfies

\[
x'=x+p^j t,\qquad
t\equiv\frac{u-x^2}{p^j}(2x)^{-1}\pmod p.
\tag{5}
\]

Every square unit therefore has exactly two square roots at each level, and every principal unit is a square. This assertion excludes nonunits: for example, 9 modulo 27 has six square roots. The square subgroup of the unit group has index 2. The production `square_lift` retains the selected residue root instead of silently changing square-root sheets.

Now require a sign on **all nonzero residues**, not just units, to be quadratic-unit-covariant:

\[
\sigma(ud)=\chi_p(u)\sigma(d)\quad(u\in U_r, d\ne0).
\tag{6}
\]

The complete classification is

\[
\boxed{\sigma(p^v a)=\varepsilon_v\chi_p(a\bmod p),\quad
0\le v<r,\quad p\nmid a,\quad\varepsilon_v\in\{\pm1\}.}
\tag{7}
\]

Indeed units act transitively on each valuation layer. The stabilizer of $p^v$ is contained in $1+p\mathbb Z$, so the prescribed quadratic sign is well-defined. The $r$ choices $\varepsilon_v=\sigma(p^v)$ are independent and determine everything. Since $\chi_p(-1)=(-1)^{(p-1)/2}$, these signs are anti-invariant under $d\mapsto-d$ exactly when $p\equiv3\pmod4$. There are then precisely $2^r$ such symbols; for $p\equiv1\pmod4$ there are none.

For example, $(p,r)=(3,2)$ and both layer signs positive give

\[
(\sigma(0),\ldots,\sigma(8))=(0,1,-1,1,1,-1,-1,1,-1).
\]

This is a complete tournament orientation by differences on $\mathbb Z/9\mathbb Z$. The naive $\chi_3(d\bmod3)$ instead gives zero at 3 and 6 and cannot orient those pairs. This is not a claim about the additive group or quadratic character of $\mathbb F_9$: the ring and field are different marked sources.

The source/target map is now explicit: a valuation-layer frame $(\varepsilon_v)$ produces an intrinsic binary relation on residue differences. It preserves covariance and antisymmetry; forgetting the frame loses orientation on nonunit differences. It does not give a plane unit-distance embedding or a colouring certificate.

## 4. Signed incidence is the exact odd-characteristic filter

For a finite simple graph $G$, let $B$ have one row $e_u+e_v$ per edge. The equations $Bf=0$ say $f_u=-f_v$. A spanning tree determines a connected component from one value $t$. On a bipartite component there is no additional condition; on a nonbipartite component an odd cycle imposes $2t=0$. Thus, if $b$ and $o$ count bipartite and nonbipartite components, including isolated vertices among the former,

\[
\#\ker(B\bmod q)=q^b\gcd(q,2)^o.
\tag{8}
\]

Over every odd-characteristic field, $\dim\ker B=b$ and $\operatorname{rank}B=|V|-b$. A component has an everywhere nonzero normalized alternating sign precisely when it is bipartite. Merely exhibiting the zero vector proves nothing. Over $\mathbb Z/6\mathbb Z$, even the nonzero vector $(3,3,3)$ solves the triangle equations; those values are not units or the signs $\pm1$.

Prescribed values add exactly the tree-path constraints. For (1), first check pinned ternary labels, then bipartiteness, then the relative binary parity of all pins in each component. If consistent, there are $2^f$ lifts, where $f$ is the number of components with no binary pin. `fixed_fibre` computes this exact count and an actual colour word. It never interprets a failed branch as a refutation of all ternary labellings.

There is a sharp spectral boundary. Over the reals,

\[
f^{\mathsf T}(D+A)f=\|Bf\|^2=\sum_{uv\in E}(f_u+f_v)^2,
\tag{9}
\]

so $\ker(D+A)=\ker B$. On a $d$-regular graph this is the usual adjacency eigenvalue $-d$ criterion. **The Gram replacement is invalid over finite fields.** For $K_4$, $B$ has rank 4 modulo 3, while

\[
D+A=2I+J\pmod3
\]

has rank 3 and kills the all-ones vector. But $B\mathbf1$ is the nonzero all-twos vector. Positivity in (9), not just formal matrix multiplication, was the missing hypothesis. Retaining the edge equations gives a reliable odd-prime filter; a modular eigenvalue alone does not.

## 5. Odd matching thickenings preserve the fixed-label support

Let $G$ be a finite simple source graph and partition the vertices of a host $H$ into fibres $F_v$, each of the same **odd** size $k$. Retain a perfect matching between $F_u$ and $F_v$ for every source edge $uv$. Assume every host edge projects to an actual edge of $G$; in particular, no host edge is collapsed within a fibre. Extra edges over old source pairs are allowed.

Fix labels $\ell_v\in\mathbb Z/3\mathbb Z$ and prescribed source colours. Pull the labels back constantly on each fibre and replicate every pin on **all** of its fibre. Then

\[
\boxed{\text{The prescribed fixed-label source fibre is nonempty}
\iff\text{the corresponding host fibre is nonempty}.}
\tag{10}
\]

**Proof.** A source colouring pulls back because every host edge projects to an old edge. Conversely let $s_w=(-1)^{b_w}$ be the signs of a proper host lift and put $S_v=\sum_{w\in F_v}s_w$. For an equal-label base edge, each retained matching edge has opposite signs. Summing over the matching gives $S_u=-S_v$. Each $S_v$ is an odd integer, hence nonzero. Its sign determines a base binary shade satisfying every equal-label edge. A uniformly replicated pin has $S_v=\pm k$ with the required sign. Different-label base edges are already safe. Equation (1) now reconstructs the required source colouring. □

The API `Cover` stores the source, $k$, and each matching as an actual permutation. `audit_thickening` permits additional host edges only after checking their source images and retaining every matching edge. `thickening_descend` checks the entire actual colour word before majority; `thickening_pullback` checks every host edge after pulling back. Neither routine assumes that arbitrary similarly sized layers form a cover.

A genuine graph cover is the special case with exactly the matching edges. For an odd simple cycle in the base, its monodromy is a permutation of an odd number of sheets and therefore has an odd-length orbit $d$. Traversing the base cycle $d$ times gives an actual simple odd cycle of length $dL$ upstairs. `odd_cycle_lift` returns and authenticates that vertex word. This is also a concrete negative certificate for a monochromatic-label branch. No infinite covering argument is used.

### Sharp losses

* **Even sheet count:** the canonical double cover of a triangle is a six-cycle. Alternating shades exist upstairs; the fibre sums can vanish. An odd prime power works because it is odd, not because primality itself supplies a new mechanism.
* **Counts:** a single edge has two unpinned binary colourings; its trivial three-cover has eight. Its complete bipartite matching thickening has two. Equation (10) preserves support, not multiplicity.
* **Pins:** in the trivial three-cover of an edge, the word $(0,3,0;3,0,3)$ pins one lift of each endpoint to colour 0 and remains proper. Pinning both base endpoints to 0 is impossible. One selected pin per fibre cannot replace uniform pin replication.
* **Variable labels:** take base $K_7$, with vertices $0,\ldots,6$, and three sheets. For each $i<j$, match $(i,t)$ to $(j,t+1)$ modulo 3. This is a genuine three-cover, with proper three-colouring $c(i,t)=t$. Its base is not six-colourable. Thus unrestricted six-colourability does not descend through odd covers. The label fibre is not constant in this example.
* **Projected contacts:** a collapsed edge or an edge over a new source pair invalidates the pullback half of (10). A failed chosen projection is not itself an equality gadget or a six-colouring refutation.

For the current HN search there is an even simpler one-way consumer: any graph retraction $H\to G$ fixing the old vertices pulls back **every** old boundary colouring, regardless of parity or matchings. The odd matching theorem supplies a two-way fixed-label statement, but requires stronger data. Increasing the number of vertices or contacts over old source edges cannot by itself defeat a certified retraction. A useful next geometry must obstruct the actual retained repair maps, and must still authenticate a separating-support loss or a genuine equality relation. No metric realization follows from the cover construction alone.

## 6. A second actual odd-degree transfer: square classes and unit circles

Let $L/K$ be a finite field extension of odd degree $d$, with characteristic different from 2. If $a\in K^\times$ becomes a square $a=b^2$ in $L$, then

\[
N_{L/K}(b)^2=N_{L/K}(a)=a^d,
\qquad
\left(\frac{N_{L/K}(b)}{a^{(d-1)/2}}\right)^2=a.
\tag{11}
\]

Thus $K^\times/(K^\times)^2\to L^\times/(L^\times)^2$ is injective. This is a genuine descent map with a retained field norm; an odd prime-power extension is one special case. Even degree fails already for $2\in\mathbb Q$ becoming a square in $\mathbb Q(\sqrt2)$.

In the marked basis $1,\theta,\theta^2$ with $\theta^3=2$, the checker independently compares the multiplication determinant with

\[
N(x+y\theta+z\theta^2)=x^3+2y^3+4z^3-6xyz.
\]

The square has coefficients $(x^2+4yz,\ 2xy+2z^2,\ y^2+2xz)$. Finite coefficient controls retain the entire vector; comparing the scalar coefficient alone would not authenticate the premise $b^2\in K$.

There is a related, distinct geometric map. Suppose $K\subseteq L\subseteq\mathbb C$ are conjugation-stable, $[L:K]$ is finite and odd, and two distinct centres $a,b\in K$ satisfy

\[
(z-a)(\bar z-\bar a)=(z-b)(\bar z-\bar b)=1
\]

for $z\in L$. Subtraction expresses $\bar z=\lambda z+\mu$ with $\lambda\ne0$ and $\lambda,\mu\in K$. Substitution gives a nonzero quadratic equation for $z$ over $K$. Since $[K(z):K]$ divides an odd integer, it must equal 1. Hence an exterior point $z\in L\setminus K$ can have at most one old $K$-centred unit neighbour. This is the two-circle mechanism used by the concurrent [surviving-repair geometry packet](surviving_repair_geometry_20261018d.md).

Equations (10) and (11) share an odd-degree obstruction to losing binary information, but their maps are different: graph matching sums versus field norms. The two-circle statement requires actual coordinates, conjugation, and a finite extension. None of these follow merely from a graph cover or a residue-ring sign symbol. Sparse attachment alone also does not prove a colouring extension when exterior vertices have additional mutual constraints.

## 7. Reproduction and bounded scope

Owned [checker](../../04-computation/experiments/odd_prime_function_cycle_20261018d.py) and [JSON receipt](odd_prime_function_cycle_20261018d.json). Reproduce with:

```text
python -B 04-computation/experiments/odd_prime_function_cycle_20261018d.py
python -B -O 04-computation/experiments/odd_prime_function_cycle_20261018d.py
```

The declared controls are all 1,099 simple graphs through five vertices for incidence rank over 3, 5 and 7; all graphs through four vertices for literal modular kernels and every ternary labelling with three pin contexts; all 216 three-sheet triangle covers, with every ternary label and three pin contexts; all 1,296 three-sheet square covers; direct finite permutation-kernel enumeration; quadratic-unit-covariant symbols for $(p,r)=(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),(7,2)$, with exhaustive odd-symbol comparison through modulus 9; 343 marked cubic-field coefficient vectors; and exact positive, loss, and malformed-source controls.

The triangle cover cycle decoder selects 144 length-3 and 72 length-9 witnesses. The all-parameter proofs are separate from these finite receipts. Public validators check exact integer types, source edges, permutation bijections, label charts, and pins; booleans and float aliases are rejected before arithmetic comparisons.

Normal and `-B -O` runs agree exactly: **58,682 checks**. The LF-normalized JSON SHA256 is `c4560eb125d30736932138e0c7ce752f72c0a350df59c4ec82b0c2fbe09e7f67`. Independent full proof/API review and optimized-versus-saved replay by `tiling_fold` passed. Its sole prose clarification was applied before freeze: the two-square-root statement concerns units, not nonunit squares. There is no claim of an all-label geometric search, a new high-chromatic unit graph, or a universal arithmetic-to-colouring correspondence.
