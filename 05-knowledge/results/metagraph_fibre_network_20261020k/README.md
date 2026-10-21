# Marked deletion fibres, quadratic orders, and their electrical quotient

2026-10-20. Research packet; no theorem ID reserved.

**PROVED, using the cited classical class-group facts:** $3m^2$ is idoneal
exactly for $m=1,2,4$. This is a classical restricted classification, with a
short proof below; it does not use completeness of the general 65-number list.
The three central-deletion fibre types admit explicit **marked** regular-action
identifications with the corresponding class groups. Their choice-state graph
has an exact weighted electrical quotient and an exact Markov quotient with
holding probabilities. Neither is an identification of the tournament metagraph.

**FINITE-EXACT:** literal path/deletion, quadratic-form, electrical and typed
receipt controls in [check.py](check.py), reproduced by [receipt.json](receipt.json).
**OPEN:** an electrical classification or a connectivity proof for the full
black tournament metagraph. The missing port information is specified in §6.

## 1. Inheritance, source types, and correction boundaries

The worktree was read at `5df3c24e1f25`. The closest proved mechanism is the
pointed Boolean product in
[zero-heavy microcosm, §§2–5](../zero_heavy_ports_20261020g/microcosm.md), extended
to arbitrary interior thinning in
[balanced deletion, §§1–3](../balanced_fibre_endpoint_extensions_20261020/BALANCED_DELETION_EXTENSIONS.md).
The source is a **literal Hamiltonian path in a labelled two-return tournament**;
the map physically deletes the specified vertices. The $4,2,1$ fibres and the
collapse of some black paths onto blue paths are inherited, not new counts.

The canonical hostile is an ordinary deletion that leaves a non-path orphan;
the balanced decoder avoids this by authenticating each joining arc. A second
near miss is replacing the full endpoint matrix by its rank-one bulk: the
inherited full matrix is unimodular, and its exceptional ports matter. The
least-used sidecar here is the deleted finite-bit address, including whether
either coordinate is the distinguished point at infinity.

[The trijection graph audit](../trijection_reset_20261020/graph_audit.md) retains
ordered vertices and literal ports. The current
[session17 synthesis, §3](../metagraph_trijections_reset_session17_20261020.md)
warns that its aggregate cut-class Cayley graph is not an equitable quotient of
the actual metagraph: successive aggregate edges can require different
representatives. Its aggregate edges lie in commuting squares, while an actual
transitive tournament has a degree-one bridge. We preserve this distinction.

The existing [idoneal information packet](../idoneal_information_channels_20261021/IDONEAL_INFORMATION_CHANNELS.md)
keeps the order conductor and the square subgroup of the class group. Its
discriminant $-252$ example has a cyclic group of order four, so class-number
four alone does not mean two independent bits. The proof below checks exponent,
not just cardinality. Incoming commits `32579ff9c2` and `4efed09d67` reserve
THM-4739/4740 and HYP-9351 for Kirchhoff/Hodge and consecutive-port work; those
reservations are not proved dependencies of this packet.

The concept board is: physical path words; deleted-bit torsors; quadratic-order
classes; within-fibre edges; energy versus random-walk clocks; port-compatible
path lifting. Every map below specifies which of these data it retains.

## 2. Why $3m^2$ gives exactly $3,12,48$

We use two classical inputs. Grube's criterion identifies idoneality of $d$
with the proper class group of discriminant $-4d$ having exponent dividing two;
see Kani, **Theorem 6 and Remark 7(b), printed p.202**. Kani's Theorem 12 also
already implies the restricted classification below. Thus no priority or
general-list completeness claim is made. [Kani (2011)](https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf)
The ring-class-number formula is stated with its unit index in **Theorem 9** of
[Chao Li's complex multiplication notes](https://www.math.columbia.edu/~chaoli/docs/MinorThesis3.html).
We also use classical reduction: each positive primitive proper form class has
one reduced representative with $b\ge0$ when $|b|=a$ or $a=c$; its inverse is
represented by $(a,-b,c)$. An ambiguous class has $b=0$, $|b|=a$, or $a=c$.

### 2.1 Odd prime factors force nonbinary class information

The order with discriminant $-12m^2$ lies in $\mathbb Q(\sqrt{-3})$, has
conductor $2m$, maximal class number one, and unit index three. Consequently

$$h(-12m^2)=m\prod_{\substack{p\mid m\\p\text{ odd}}}
 \left(1-\frac{\chi_{-3}(p)}p\right). \tag{1}$$

Write $m=2^a3^b\prod_{p\ge5}p^{e_p}$. Then

$$h(-12m^2)=2^a3^b\prod_{p\ge5}p^{e_p-1}(p-\chi_{-3}(p)). \tag{2}$$

For $p\ne3$ odd, $\chi_{-3}(p)=1$ for $p\equiv1\pmod3$ and $-1$ for
$p\equiv2\pmod3$. Thus $6\mid p-\chi_{-3}(p)$. If $m$ has any odd prime
factor, $3\mid h$. A finite group of exponent dividing two has power-of-two
order, so idoneality forces $m$ to be a power of two.

### 2.2 An explicit inverse pair excludes every larger power of two

Let $m=2^a$, $a\ge3$. Choose the unique
$b\in\{-6,-4,-2,2,4,6\}$ with $b\equiv3m\pmod7$. The form

$$f_m=\left(7,b,\frac{b^2+12m^2}{28}\right) \tag{3}$$

is integral: the numerator is divisible by four because $b$ is even, and by
seven because $b^2+12m^2\equiv21m^2$. It is primitive since $7\nmid b$.
For $m\ge8$ it satisfies $0<|b|<7<c$. Therefore it and its inverse are
distinct strictly reduced forms. The class has order greater than two, which
excludes idoneality. The first witness is $(7,-4,28)$ at $m=8$.

The remaining complete lists are

| $m$ | $d=3m^2$ | Discriminant | Reduced forms | Class group |
|---:|---:|---:|---|---|
| 1 | 3 | $-12$ | $(1,0,3)$ | trivial |
| 2 | 12 | $-48$ | $(1,0,12),(3,0,4)$ | $C_2$ |
| 4 | 48 | $-192$ | $(1,0,48),(3,0,16),(4,4,13),(7,2,7)$ | $C_2^2$ |

All listed forms are ambiguous. The reduction bound
$a\le\sqrt{|D|/3}$ proves the lists complete by a finite check. This proves
the iff for every positive integer $m$, not just the tested range.

## 3. The actual $4,2,1$ fibres and a marked arithmetic bridge

For $k\ge2$, let $T_k$ have vertices $0,\ldots,2k$, all arcs increasing
except $k\to0$ and $2k\to k$. Put

$$X_k=\mathcal P(\{1,\ldots,k-1\})\sqcup\{\infty\},
\qquad q=|X_k|=2^{k-1}+1.$$

Represent $\infty$ by the full set $A=\{0,\ldots,k-1\}$. With
$B=\{k+1,\ldots,2k\}$ and $\sigma(v)=2k-v$, the inherited decoder is

$$P(I,K)=I^\uparrow,\ [\sigma(A\setminus K)]^\uparrow,\ k,
 (A\setminus I)^\uparrow,\ [B\setminus\sigma(A\setminus K)]^\uparrow. \tag{4}$$

It is a bijection $X_k^2\to\mathrm{HP}(T_k)$, with blue predicate $I=K$.
Physical central deletion from $T_{k+1}$ removes vertices $k,k+2$, relabels
the survivors increasingly, and on coordinates forgets the last finite bit
while fixing infinity. The one-coordinate weights are

$$w_x=2\ (x\ne\infty),\qquad w_\infty=1,$$

so the product fibre over $(x,y)$ has $w_xw_y$ points.

A fibre with $r$ infinity coordinates has a regular action of
$\mathbb F_2^{2-r}$, by flipping the deleted finite bits. Choose its zero-bit
lift and order its active coordinates. Set $m=2^{2-r}$. The table in §2 then
gives a marked equivariant bijection from the fibre to
$\mathrm{Cl}(-12m^2)$. Explicit generators are $(3,0,4)$ at $m=2$, and
$(3,0,16),(4,4,13)$ at $m=4$; their product is the remaining nonidentity class
because the latter group is $C_2^2$.

**Preserved:** bit-flip composition and the regular action. **Lost unless
retained separately:** actual path, ports, blue status, quadratic-form
representation values, and any metagraph adjacency. There is no canonical
identity in an unmarked fibre. In particular, equal sizes alone would not
distinguish a $C_4$-torsor from a $C_2^2$-torsor.

For a literal hostile, the path

$$P(\varnothing,\{2\})=(5,6,3,0,1,2,4)\quad\text{in }T_3$$

is black, while deleting vertices $2,4$ yields the blue path
$(3,4,2,0,1)$ in $T_2$. A finite diagonal parent has four lifts, two blue and
two black. A black parent has only black lifts; the converse fails.

## 4. Exact electrical and stochastic quotients of the choice graph

Modulo coordinate swap, black paths are the edges of the complete graph on
$X_k$. These graph vertices are **choice states**, not tournament vertices or
isomorphism classes of tournaments. No intrinsic orientation is imposed.

For any finite positive-conductance graph with Laplacian $L$ and a surjection
$\pi$, let $P_{v,x}=1_{\pi(v)=x}$. Wiring each fibre gives

$$L_{\rm wired}=P^{\!T}LP,\qquad
 (Pu)^{T}L(Pu)=u^TL_{\rm wired}u. \tag{5}$$

Each coarse conductance is the sum of the actual fine conductances between
the two fibres. Edges inside a fibre vanish from (5), so their labels must be
kept if their count or role matters. This construction composes exactly under
successive surjections. It restricts potentials to be fibre-constant; it is
not a Schur complement that minimizes over internal potentials.

Apply this to the unit-conductance choice graph $K_{X_{k+1}}$. Put
$W=2q-1=\sum_xw_x$. The between-fibre conductances and Laplacian are

$$c_{xy}=w_xw_y,\qquad L=W\operatorname{diag}(w)-ww^T. \tag{6}$$

Exactly $q-1$ unoriented black edges are internal to fibres and collapse.
They account for the inherited black correction. For a potential $u$,

$$\mathcal E(u)=W\sum_x w_xu_x^2-\left(\sum_xw_xu_x\right)^2. \tag{7}$$

Inject one unit of current at $x$ and remove one at $y$. Equation (6), or a
weighted-mean-zero solution, gives

$$R_{xy}=\frac1W\left(\frac1{w_x}+\frac1{w_y}\right). \tag{8}$$

Thus finite/finite and finite/infinity resistances are $1/W$ and $3/(2W)$.
These are responses between **wired fibres**. They are not the fine graph's
resistance $2/W$ between two specified individual vertices.

There is an exact stochastic statement as well. The ordinary discrete simple
walk on $K_W$ strongly lumps to

$$Q_{xy}=\frac{w_y-1_{x=y}}{W-1}. \tag{9}$$

Every fine state in fibre $x$ has precisely $w_y-1_{x=y}$ neighbours in fibre
$y$, which proves the assertion and retains the original step clock. The
stationary distribution is $w_x/W$. At a finite parent there is holding
probability $1/(W-1)$; at infinity there is none. The loopless conductance
walk of (6) instead jumps from $x$ to $y$ with probability
$w_y/(W-w_x)$. It has suppressed these holds and is a different clock.

Finally, a uniform fine **ordered path state** pushes forward to parent mass
$w_xw_y/W^2$, not $1/q^2$. Its conditional lift entropy is

$$H(\text{fine path}\mid\text{parent path})
 =\frac{4(q-1)}{2q-1}\ \text{bits}. \tag{10}$$

This uses only $\log_2 2=1$: each finite coordinate has two equiprobable
lifts. It is an information loss of the marked deletion channel, not a
thermodynamic interpretation of tournament intransitivity.

## 5. Why a weighted quotient does not settle metagraph connectivity

Let the fine graph have edges $a-b$ and $c-d$, all with conductance one, and
identify $b,c$ into one coarse vertex. The wired graph is a two-edge path,
although the fine graph is disconnected. The wired resistance from $a$ to
$d$ is two while the actual resistance is infinite. Add the fine edge $b-c$
with conductance $g>0$: the wired graph is unchanged, but the actual resistance
is $2+1/g$. Even a connected fibre does not make wiring preserve resistance.

This is the first failed implication: an aggregate edge gives an existential
representative, while two consecutive edges may require incompatible
representatives in the shared fibre. Neither fibre cardinality nor the sum of
conductances determines that compatibility. General graph aggregation is not
strong lumpability either: fine vertices in one fibre may have different
transition rates to another fibre. Equation (9) holds because complete-graph
neighbour counts are constant, not merely because a surjection exists.

The strongest elementary survivor is useful: if every fibre induces a
connected graph, every coarse path lifts to a fine walk. Select the actual
fine edge realizing each coarse step and join its endpoints within the
intermediate fibre. Thus connectivity is equivalent under this additional
hypothesis. With quantitative internal connector costs one gets path-length
bounds, but not equality with the zero-cost wiring energy.

## 6. A concrete remaining metagraph interface

For a proposed aggregation of an actual black metagraph, store:

1. every representative tournament/path state and its literal black-move edge;
2. for each coarse incident edge, the fine endpoint inside the fibre;
3. the internal connectivity relation on these fine ports, or explicit
   connector paths; and
4. if electrical responses are desired, the fibre's boundary response matrix,
   obtained by eliminating only its genuine internal vertices.

A coarse walk is realizable exactly when representative edges can be chosen
whose adjacent fine ports are connected by these retained within-fibre paths.
This is a finite compatibility test, not an appeal to the coarse degree. The
$a,b,c,d$ hostile proves the port relation is necessary. Its conductance-$g$
variant proves that connectivity alone is insufficient for exact resistance.

The next positive obligation is therefore to prove a source-bound internal
connector theorem, or compute these boundary relations on the relevant black
fibres. No such connector theorem for the full metagraph is assumed here.
The inherited choice graph happens to have complete fibres and hence meets the
connectivity condition; that fact does not transfer to arbitrary metagraph
representatives by matching their fibre sizes.

## 7. Reproduction and finite scope

Run from the repository root:

```text
python -B -X utf8 05-knowledge/results/metagraph_fibre_network_20261020k/check.py --verify
python -B -O -X utf8 05-knowledge/results/metagraph_fibre_network_20261020k/check.py --verify
```

The receipt contains **6,770 always-active checks**, LF SHA256
`4179549c14b73beec5672fc7511a7daa23143a6c3ad4664add1450bb61ab85bd`.
The finite universe comprises reduced forms for $1\le m\le64$, odd-factor
divisibility for $m\le2000$, the explicit power-of-two hostile for exponents
$3\le a\le100$, literal balanced source/target deletion for $2\le k\le5$,
and exact choice networks at $q=3,5,9$. Independent DFS verifies the complete
small path decoder; rational nodal elimination verifies (8) independently of
its displayed formula. All failures use exact integers or fractions; public
inputs reject bool/float aliases and wrong source words.

The proof, rather than these finite ranges, supplies the all-$m$ arithmetic
classification and the arbitrary finite-surjection energy theorem. No
Hamiltonian-path search is used to assert a theorem about all tournament
metagraphs.

Independent peer audit: **PASS** (family algebra), full note/core/API and a
separate optimized replay equal to the saved receipt. The audit checked the
marked torsor, all-$m$ arithmetic argument, held Markov quotient, and
connected-fibre walk-lift boundary; it did not re-prove the cited classical
class-group inputs. The root also independently reproduced the 6,770-check
receipt and checked the arithmetic proof. The three owned files are frozen.
