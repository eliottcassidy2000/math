# Matchings determine slot demand; triangle wraps test retention

**Proved.** In a one-port tournament substitution, the maximum number of nonempty module pieces permitted by the surroundings in an ordered $r$-cover is exactly $r+\nu$, where $\nu$ is a bipartite matching size. For a specific module $A$ the exact maximum is $\min(|A|,r+\nu)$. Separately, reading the two-slot count through a bounded number of triangle wraps detects every unequal same-order slot polynomial. An equal-response pair with Hamiltonian count gap 6 cannot preserve its response through even one common triangle wrap.

## Random branch and the transfer that survived

Use local tree 1b91742ce55a202b73c80d56561ba56d00b503a6. Sort the Markdown files under 05-knowledge/results, exclude paths containing chanpak, chan-pak, mistake, audit, index or readme, and sample 12 with Python Random(740689):

    session25_spectrum_wiring_20261021/session25-odd-girth-review/REVIEW.md
    collatz_guard_cover_defects_20261007e.md
    lrc14_literature_kps.md
    hyp9310_is_klmsy_corollary_reflection_proof_20261018.md
    overnight2_20260906_no3line_third.md
    session26_random_bridges_20261022/session26-wildcard-review/REVIEW.md
    planar_jc48_sep06_dg_surface.md
    collatz_mod6_20260922_lead_reflection3.md
    long_frontier_sep06_dimension.md
    three_problems_marked_defects_20261012a.md
    modular_return_debt_20261004.md
    crossroads10_20260926_graph.md

The [Walecki odd-girth audit](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session25_spectrum_wiring_20261021/session25-odd-girth-review/REVIEW.md) uses disjoint repair pairs to isolate compatibility costs. This suggested using a matching rather than counting all the triangles at a port. The [Collatz guard ledger](../collatz_guard_cover_defects_20261007e.md) distinguishes an interchangeable guard from interchangeable future obligations. Its precise analogue here is that equal $\Pi_2$ does not imply equal continued response. The [233/223 edge-flip note](../crossroads10_20260926_graph.md) makes the same point about losing mixed interactions. The planar-Jacobian boundary example and the reflection/CM note were exploratory lenses, not lemmas used below.

## 1. Exact demand of a single port

Let $K$ be a tournament with marked vertex $x$, $O=N^+_K(x)$, $I=N^-_K(x)$. Let $\nu$ be the maximum matching size among arcs $o\to i$ with $o\in O,i\in I$. These are the triangles $x\to o\to i\to x$ with the common vertex removed. Substitute a nonempty module $A$ for $x$. Let $\Pi_r$ count ordered $r$-slot path covers allowing empty slots, and $P_j$ count ordered covers by $j$ nonempty paths. Then

$$\Pi_r(K[A])=\sum_{j=1}^{r+\nu} c_{r,j}(K,x)P_j(A),\qquad c_{r,j}>0\quad(1\le j\le r+\nu).$$

Compress each maximal visit to $A$ to $x$. Outside vertices appear once; no slot has adjacent copies of $x$. A template with $j$ appearances of $x$ has exactly $P_j(A)$ lifts by assigning the paths in order. This is the inherited visit-matrix argument; the exact demand and its sharpness are new here. Every $P_j(A)$ with $1\le j\le |A|$ is positive by cutting a Hamiltonian path; thus restricting the template bound by $|A|$ is sufficient.

* **Upper bound.** Each gap between consecutive appearances of $x$ in one slot is an outside path from $O$ to $I$. It contains an arc $O\to I$. Gaps consume disjoint outside vertices, so choosing one crossing per gap gives a matching. At most $r$ appearances are not paired with a previous appearance in their slot.
* **Lower bound, $j\ge r$.** Choose $j-r$ matched edges, concatenate the excursions $x,o,i,x$ in the first slot, and put one more copy of $x$ in each of the other slots. Treat copies of $x$ as distinct vertices with identical arcs to outside vertices and arbitrary arcs among themselves. Insert each remaining outside vertex into the first path by the elementary tournament insertion lemma. Insertion cannot make two copies adjacent. Compress the copies back to $x$.
* **Lower bound, $j<r$.** Put a Hamiltonian path of $K$ in one slot, put $j-1$ extra singleton copies in other slots, and leave the rest empty.

When $|K|-1=2\nu$, the top coefficient is the number of perfect matchings times $\nu!\binom{\nu+r-1}{r-1}$: each gap has length two, every slot contains $x$, and ordering a perfect matching then making a weak division into $r$ strands specifies the cover. With extra outside vertices the same factorial/binomial factor is a lower bound by extending one fixed matching; erasing the unused vertices recovers its division and order.

**Consequence for repairs.** For $A=TT_s$, $P_j(A)$ counts surjections from $s$ vertices to $j$ labels, so the nonzero top term gives exponential base $r+\nu$. If $\nu>0$, no universal linear expression in the first $r$ coordinates $\Pi_1,\ldots,\Pi_r$ can update $\Pi_r$, since those coordinates are $1^s,\ldots,r^s$ on the test modules. This need not rule out a specially chosen family or a nonlinear update. If $\nu=0$ then every vertex of $I$ beats every vertex of $O$. The substitution is the ordinal join $K[I]\Rightarrow A\Rightarrow K[O]$, so $H$ and $\Pi_r$ simply multiply by their side factors. Under nested one-port substitutions the maximal demand is $r$ plus the sum of the matching sizes, because every top coefficient is positive.

## 2. What a single triangle costs an equal-response pair

Use the inherited triangle reader $C_3[A,1,1]$. Its operator on $f(r)=\Pi_r(A)$, with $f(0)=0$, is

$$(Lf)(r)=r[f(r+1)+r f(r)-f(r-1)].$$

For a direct count, retain a variable $z$ marking visits to the port. The two singletons occupy different slots in $r(r-1)(1+z)^r$ templates; their same-slot polynomial is $r(1+3z+z^2)(1+z)^{r-1}$. Summing gives $r[(1+z)^{r+1}+r(1+z)^r-(1+z)^{r-1}]$, and $\Pi_j=\sum_k\binom jkP_k$ supplies the operator.

For two modules of the same order $n$, $D(r)=\Pi_r(A)-\Pi_r(B)$ is an even-coefficient integer polynomial of degree at most $n-2$, with exponents of parity $n$. By the ordinary odd-cycle identity on each of the $r$ color classes,

$$\Pi_r(T)=\sum_{\text{vertex-disjoint odd cycle collections }C}2^{|C|}r^{n-|V(C)|+|C|}.$$

The empty collection contributes $r^n$ and all others contribute even coefficients in degrees of the stipulated parity. This is [THM-4756](../../../01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md), step 2.

Set $\epsilon=1$ at odd order and $\epsilon=2$ at even order. If $\Pi_2$ agrees then $D(r)=2r^\epsilon(r^2-4)G(r^2)$ for integral $G$. If it still agrees after one triangle, $(LD)(2)=0$ implies $D(3)=D(1)$. But $D(1)=-6G(1)$ and $D(3)=10\cdot 3^\epsilon G(9)$. Thus the initial $H$-gap must be divisible by **30**. At odd order it must be divisible by **120**: $G(1)+5G(9)=0$ and $G(9)\equiv G(1)\pmod8$, so $4\mid G(1)$. In particular an initial gap 6 preserves that ordinary gap but forces a next response difference of $-12$ modulo 60 (odd order) or modulo 180 (even order).

More generally, matching $\Pi_1,\ldots,\Pi_R$ at equal order forces $D$ divisible by $2r^\epsilon\prod_{i=1}^R(r^2-i^2)$. Therefore the next coordinate difference is divisible by $2(R+1)^{\epsilon-1}(2R+1)!$. For $R=2$ at even order this is sharp: order-eight codes 87558995 and 249732159 have $(H,\Pi_2)=(345,5640)$ but $\Pi_3=45225,44505$. In this order $D(r)=2r^2(r^2-1)(r^2-4)$. A triangle wrap preserves both counts at 5985 but separates their response by 1440. Codes use increasing lexicographic pairs, first pair least significant, bit 1 meaning the smaller vertex wins.

## 3. Bounded-depth test for the entire same-order profile

Suppose $(L^iD)(2)=0$ for $i=0,\ldots,k$. These equations determine $D(2),\ldots,D(k+2)$ from $D(1)$ and $D(0)=0$: the coefficient at $D(i+2)$ in equation $i$ is $2\cdot3\cdots(i+1)\ne0$. Let $g$ solve the formal eigen-sequence $Lg=g$ with $g(0)=0,g(1)=1$. The equation at 1 gives $g(2)=0$ and the later ones give

$$g(r+1)=g(r-1)-(r-1/r)g(r),\qquad g(3)=1,\qquad g(4)=-8/3.$$

The data must agree with $D(1)g$ through index $k+2$ by triangular uniqueness. Starting at $3$, signs of $g$ alternate and magnitudes increase, since each next magnitude is at least $(r-1-1/r)|g(r)|>|g(r)|$. If $D(1)\ne0$, the values at $3,\ldots,k+2$ alternate and force $k-1$ positive roots there, in addition to the root at 2, for $k\ge1$. If $D(1)=0$ then all the integers through $k+2$ are roots instead. Parity reflects the positive roots and contributes the zero at the origin with multiplicity at least $\epsilon$. A nonzero $D$ therefore needs degree at least $2k+\epsilon$ (also true at $k=0$). Thus $2k+\epsilon>n-2$ forces full profile equality. A fixed physical test family distinguishes all different same-order profiles in bounded depth.

Run python verify_matching_retention.py. The verifier independently enumerates all compressed-cover supports on 1,099 labeled pointed small tournaments, checks the saturated leading coefficients and literal triangle wraps, then recounts the colliding order-eight profiles and their separation. The open step is a repair with the specific extra data its next reader needs, not merely one with the same current $\Pi_2$.
