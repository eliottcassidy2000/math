# The next Paley run: a wheel in A10 and a six-vertex terminal in R10

2026-10-14. **PROVED** guarded reduction identities and terminal classification;
**FINITE-EXACT** symbolic coefficient extraction and independent signed-walk
controls; **CITED** the inherited prime wheel evaluation. No uniform asymptotic
remainder or new modular interpretation of the six-vertex terminal is claimed.

## 1. Inheritance and the new statements

The closest mechanism is the exact degree-one/two injective reducer in
[the first wheel/run interface](wheel_run_interface_20261014h.md), which proves
$A_8=2q(q-1)(q-3)(q-7)(7q-55)$. The current correction lineage includes
MISTAKE-061: a cycle-only classification misses even-series/theta patterns.
The next hostile is more specific: even a correct classification of **raw
path quotients** does not classify the graphs produced by subsequent
degree-two elimination. The sidecars are the original graph, vertex count,
orientation signs, injectivity, and the finite-field congruence $q\equiv3\pmod4$.

The concept board is: marked runs, occupied-value subtractions, signed graph
symmetry, polynomial factors, arithmetic terminals, and the independent
subset-DP/three-matrix evaluators. The anchor is $A_{10}$; the broader terminal
compiler is the niche; the factor interpretation is the wildcard.

Let $q\equiv3\pmod4$ be an odd prime power, and extend its quadratic character
by $\chi(0)=0$. For an oriented graph $G$ with labelled vertex set $[n]$, put

$$I_G(q)=\sum_{x:[n]\hookrightarrow\mathbf F_q}
                  \prod_{(a,b)\in E(G)}\chi(x_b-x_a). \tag{1}$$

An isolated vertex remains part of the injection domain. A repeated directed
edge is a repeated factor, not a distinct assignment constraint. Let $A_k$
be (1) for the directed path of $k$ edges. Let $U_q$ be the **unanchored,
injective** four-wheel in the preceding packet: hub 0, spokes $0\to i$,
and rim $1\to2\to3\to4\to1$.

The first new single-run formula is

$$\boxed{A_{10}(q)=
 2q(q-1)(q-3)(21q^3-576q^2+5300q-16198)-34U_q.} \tag{2}$$

This holds for every odd prime-power field $q\equiv3\pmod4$, including fields
too small to admit the injection. For an odd prime, the independently proved
[wheel-to-Apéry interface](wheel_apery_extension_20261014h.md) gives
$U_q=-q(q-1)(a_q+4)$, where
$\sum a_nz^n=z\prod_{m\ge1}(1-z^{2m})^4(1-z^{4m})^4$. Thus

$$A_{10}(q)=q(q-1)
 \left[2(q-3)(21q^3-576q^2+5300q-16198)+34(a_q+4)\right]. \tag{3}$$

The coefficient of this arithmetic term is nonzero. We do **not** infer from
that fact alone that no polynomial in $q$ could agree on some chosen domain;
such a non-polynomiality assertion would need an additional theorem about
the modular coefficients. Equation (2) is the exact, unconditional transfer.

There is also a new terminal in the complete degree-ten run coefficient.
Define $V_q=I_{K_6-P_6}(q)$ using vertices $0,\ldots,5$, deleting the five
edges $01,12,23,34,45$, and orienting every retained edge from its lower
label to its higher label. If

$$R_q(t)=\frac1{q!}\sum_{\sigma\in S(\mathbf F_q)}
  \prod_{i=1}^{q-1}(1+t\chi(\sigma_{i+1}-\sigma_i)),$$

then for $q\ge11$,

$$\boxed{[t^{10}]R_q(t)=
 \frac{P_{10}(q)+(7q^2-173q+1166)U_q+64V_q}{(q)_{10}},} \tag{4}$$

where $(q)_{10}=q(q-1)\cdots(q-9)$ and

$$\begin{aligned}
P_{10}(q)=\frac{q(q-1)(q-3)}{120}\bigl(&q^7-71q^6+2308q^5-44330q^4\\
 &+539149q^3-4096079q^2+17639622q-32287320\bigr).
\end{aligned}$$

$V_q$ is a necessary nonzero **formal terminal for the present reduction
rules and this coefficient**. We do not claim it is algebraically independent
of the wheel, or that no stronger character identity could reduce it further.
That is the next arithmetic question, with its exact graph now retained.

## 2. Guarded reduction and complete terminal classification

The earlier identities extend to coefficients that are polynomials multiplying
each terminal, rather than just constant wheel coefficients.

* A loop gives zero. Under injectivity, two copies of the same directed
  factor cancel, while reversing one copy costs minus one.
* An odd number of surviving edges gives zero by the simultaneous bijection
  $x_v\mapsto-x_v$.
* Removing an isolated vertex from an $n$-vertex source multiplies **every**
  remaining coefficient by $q-n+1$.
* A degree-one vertex has full-field character sum zero. Subtract all values
  occupied by non-neighbours; neighbour values already give zero.
* At degree two, orient both factors outwards, retaining the sign. Distinct
  neighbours give full-field sum $-1$; subtract all other occupied values as
  in the degree-one case.

Each branch decreases the number of vertices; parity reduction never increases
the edge count. Consequently every finite source graph with at most ten factors
reduces to terminal graphs of minimum degree at least three and at most ten
edges, together with polynomials from eliminated vertices. A terminal has at
most six vertices. The following classification is complete.

| Vertices, edges | Possible terminal | Reason for zero or survival |
|---|---|---|
| 4,6 | $K_4$ | Transposing two vertices reverses the Vandermonde sign. |
| 5,8 | $W_4$ | Survives; its complement is two disjoint edges and an isolated vertex. |
| 5,9 | $K_5-e$ | Odd edge count. |
| 5,10 | $K_5$ | Sign-reversing transposition. |
| 6,9 | $K_{3,3}$ or triangular prism | Odd edge count. The complements are $C_3\sqcup C_3$ or $C_6$. |
| 6,10 | complement of $C_5\sqcup K_1$ | A reflection reverses the five-cycle rim; the spokes are merely permuted. |
| 6,10 | complement of $C_4\sqcup K_2$ | Interchanging the ends of one retained matching edge reverses that factor and permutes all others. |
| 6,10 | complement of $C_3\sqcup P_3$ | Interchanging the ends of the sole edge in the three-vertex part reverses its factor. |
| 6,10 | complement of $P_6$ | Survives the sign test; the path reversal reverses ten retained edges, an even number. |

For the last four cases, the complement has five edges and maximum degree two.
Its degree sequence is either $(2,2,2,2,2,0)$ or $(2,2,2,2,1,1)$, giving exactly
the displayed cycle/path possibilities. On five vertices the same degree
argument gives the stated three cases. This also explains why there are
exactly ten isomorphism types in the census. The executable independently
enumerates all **634 labelled** terminal graphs and recovers those ten types.

Therefore every source in this edge range has an exact expression

$$I_G(q)=P_G(q)+B_G(q)U_q+C_G(q)V_q,\qquad P_G,B_G,C_G\in\mathbf Z[q]. \tag{5}$$

The source graph can have isolated vertices; those force polynomial, rather
than scalar, terminal coefficients. For example, an isolated vertex added to
the wheel gives $(q-5)U_q$. The public compiler is explicitly capped at twenty
source vertices and ten factors, covering every union of positive runs with
total length at most ten. The proof of the reduction identities itself does
not impose the twenty-vertex cap.

The code's canonical six-vertex representative has value **$-V_q$** in the
natural orientation used in (4). This sign is stored and checked, not inferred
from unsigned graph isomorphism. The canonical wheel has value $U_q$.

The same proof gives a finite hierarchy at every fixed edge budget $E$.
After repeated-edge cancellation and degree-zero/one/two elimination, each
nonempty terminal is simple, has minimum degree at least three, and has at
most $E$ edges. Therefore it has at most $\lfloor2E/3\rfloor$ vertices.
There are finitely many such signed graph types; odd-edge and sign-reversing
automorphisms remove further types. Every finite source with at most $E$
factors is consequently a $\mathbf Z[q]$-linear combination of this finite
terminal list and the empty graph. The coefficient polynomials retain the
source's eliminated vertices. This is recursive compression into a finite
core hierarchy, not a claim that one polynomial law evaluates every level.
Only the $E\le10$ level is implemented and classified in this packet.

## 3. Exact coefficient extraction

All odd-length marked runs vanish by reversing that run's ordered vertex
assignment. The seven multisets of even run lengths totaling ten are listed
below. The program retains the full polynomial part of each joint sum in
the [JSON](paley_next_cores_20261014i.json); this table displays the arithmetic
coefficients in the natural $U,V$ orientations.

| Runs | Weight in $S_{10}$ | Coefficient of $U$ | Coefficient of $V$ |
|---|---:|---:|---:|
| $(10)$ | 1 | $-34$ | 0 |
| $(8,2)$ | 1 | $1744-112q$ | 0 |
| $(6,4)$ | 1 | $1396-92q$ | 0 |
| $(6,2,2)$ | $1/2$ | $-24720+2084q-28q^2$ | $-384$ |
| $(4,4,2)$ | $1/2$ | $-20248+1418q-6q^2$ | 0 |
| $(4,2,2,2)$ | $1/6$ | $261504-23040q+384q^2$ | $3072$ |
| $(2,2,2,2,2)$ | $1/120$ | $-2764800+254400q-4800q^2$ | $-30720$ |

The weights are $1/\prod_\ell m_\ell!$, where $m_\ell$ is the multiplicity
of a run length. To justify them, a prescribed ordered list of $r$ runs with
$E$ edges has $\binom{q-E}{r}$ placements, and its normalized injection sum
has denominator $(q)_{E+r}$. Their ratio is $1/(r!(q)_E)$. Sum over distinct
orders of the same multiset. If $E+r>q$, the injection sum is zero, so these
terms contribute zero to the simplified formula as well. For $E=10$ this
proves the weighting in (4). The six-vertex coefficient is
$-384/2+3072/6-30720/120=64$, so it does not cancel.

The coefficient tables come from a terminating recursion with exact integer
polynomial arithmetic, not interpolation in $q$. Its 660,689 reduced states
for the complete ten-edge run calculation compress to the three source-aware
channels in (5). The graph/type classification and the induction on the number
of vertices prove the universal validity of each reduction; the saved exact
calculation supplies the displayed finite coefficient data.

### Why raw Euler parity was insufficient

A quotient of one oriented trail has at most two odd-degree vertices; the
wheel and $K_6-P_6$ each have four. Thus neither is a direct ten-edge trail
quotient in that sense. But a degree-two full-sum term deletes two incident
edges and changes the parity of both neighbours. It need not remain a trail
quotient. The nonzero coefficient $-34U_q$ in (2) is the decisive hostile to
using raw Euler parity as an invariant of the entire reduction process.

No such surprise is hidden by the compiler: every occupied-value subtraction
and orientation sign remains in its exact source receipt.

## 4. What the factors mean

The exact formulas for $A_2,A_4,A_6,A_8,A_{10}$ have leading coefficients
$1,2,5,14,42$ in degrees $2,3,4,5,6$ respectively (for the last assertion,
along odd primes, use the inherited bound $a_q=O(q^{3/2})$). These are the
first five Catalan numbers. The familiar rooted binary split
$C_{k+1}=\sum_{i=0}^k C_iC_{k-i}$ is a precise addition/multiplication
recursion: choose a split, then independently choose its two subtrees.
This identifies the observed leading coefficients; it does not identify
the entire signed graph sum with a tree count, or prove an all-k leading
law from five cases. The retained wheel is an explicit lower-order
coordinate that the leading coefficient does not remember.

For any source with at least two vertices and an even number of factors,
translation and nonzero scaling act freely on its injective assignments and
preserve the weight. This explains the integer divisibility by $q(q-1)$.
The symbols $q=0,1$ are formal roots of displayed polynomials, not actual
finite fields in this domain.

In the $A_8$ formula, $q=3$ and $q=7$ are real finite-field boundary cases:
nine distinct path vertices cannot fit. Once polynomiality is established
for that run, these zeros give factors $(q-3)(q-7)$. The remaining factor
$7q-55$ is an exact result of the elimination; no field or graph threshold
at $55/7$ is asserted.

At degree ten, there is no valid demand that the **polynomial part alone**
have every small-field zero. At $q=7$ it is $-39984$, while $U_7=-1176$;
the full sum is $-39984-34(-1176)=0$. Forgetting the terminal would violate
injectivity even at this smallest useful test.

Likewise, zero signed sum is not absence of assignments. The directed
triangle sum vanishes over these fields, but it has $q(q-1)(q-2)$ injective
assignments. Replacing a squared character by one in an **all-map** sum is
also invalid: at $q=7$, the two-vertex doubled-edge all-map sum is 42, not
49. Finally, the $q\equiv3\pmod4$ sign sidecar is essential: $A_2(5)=-20$,
whereas the present polynomial $q(q-1)$ would give $+20$ outside its domain.

## 5. Independent numerical evaluators and reproduction

The script independently enumerates anchored injections for every terminal
representative at $q=3,7,11$, and every simple source on at most four vertices
at $q=3,7$. It compares the older eight-edge formulas without changing that
module. An independent signed self-avoiding-walk subset DP gives

| $q$ | $A_8$ | $A_{10}$ |
|---:|---:|---:|
| 11 | 154880 | 478720 |
| 19 | 10243584 | 226059264 |

This DP fixes its first vertex at zero and retains the used subset and last
vertex. Every transition multiplies by the actual quadratic character. Its
NumPy/Numba `int64` scope is explicitly $q\le19$, where all absolute partial
and total sums are bounded by $19!<2^{63}$. These are exact integer controls,
not floating-point sampling. The complete $q=11$ Hamiltonian product DP also
gives $[t^{10}]R_{11}=34/2835$, agreeing with (4).

The parent [three-matrix trace compiler](pairwise_core_trace_20261014i.md)
and its [script](../../04-computation/experiments/pairwise_core_trace_20261014i.py)
independently retain every cross-block inequality and character factor.
Its separate trace and anchored-injection computations agree on
$V_7=336$, $V_{11}=-12320$, $V_{19}=-76608$, and $V_{31}=52080$.
This is a precise computational compression of a signed graph sum; it is not
a replacement of the six-vertex source by an unmarked scalar.

```text
python -B -X utf8 04-computation/experiments/paley_next_cores_20261014i.py
python -B -O -X utf8 04-computation/experiments/paley_next_cores_20261014i.py
```

`compile_receipt`/`audit_receipt` retain the vertex count and every oriented
factor, reject Boolean/float aliases and altered/source-replaced packets,
and recompute the exact reduction. Loops, doubled/opposite edges, isolates,
odd-edge vanishing, small-field cancellations and the residue-class boundary
are explicit controls. The remaining arithmetic target is a source-preserving
evaluation of $V_q$ beyond the finite-field summation, followed separately by
the still-open uniform control of the full Hamiltonian run expansion.

Normal and optimized executions agree on all **278 always-active checks** and
the complete saved JSON. Its LF-normalized SHA256 is
`db3a490b59be10b127f0c7b0920079879c09512a39e922d453350935621e5670`.
