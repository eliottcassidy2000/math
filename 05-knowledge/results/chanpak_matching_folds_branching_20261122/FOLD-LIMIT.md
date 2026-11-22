# Typical signed folds: exact second moment and a uniform CLT

**Status: PROVED.** This is about reversal-fixed paths, not the ordinary
Hamiltonian-path count. It complements the repository's extremal work
(Session29 `session29-halfpath-asymptotics/REPORT.md`).

Fix an anti-involution on 2m or 2m+1 vertices and write F for the number
of directed Hamiltonian paths fixed by applying the anti-involution and
reversing path order. On the m antipodal pairs, fix an owner tournament D:
for each owner edge i→j, the sign `x_ij` at i beats both signs of j,
and its opposite loses to both. The terminal sign `b_i` of pair i is
the winner of the internal arc at even order, or the sign beating the
fixed center at odd order. The sign variables x and b are independent
fair signs; optional odd-order internal arcs do not affect F.

We use the [signed-fold Walsh coefficient rule](SIGNED-FOLD.md#2-exact-walsh-coefficients-alternating-linear-forests). For completeness, an order of the pairs contributes the product
over runs of owner-backward steps: a run ending at interior j,i,k yields
`1−x_ij x_ik`, and one ending at terminal j,i yields `1−x_ij b_i`.
Each backward run frees one sign and imposes a simultaneous pair of
constraints at its end. On expansion, a selected term prescribes a
centered turn with its two neighbors (or one neighbor and the terminal).
The required adjacency graph must be a linear forest; contract its
components into blocks, order the blocks, and allow two orientations
per nontrivial block unless it is pinned to the last position by a
terminal. With r selected turns, c ordinary blocks (i.e. counting
singletons but not the special terminal block), and p ordinary
nontrivial blocks, the exact coefficient is `(−1)^r 2^p c!`. The
empty term gives `E F=m!`. A graph without a terminal has only paths
with an even number of edges. The terminal component has an odd
number of edges and its endpoint indicates its single terminal variable.

## 1. A uniform conditional CLT

Let the owner and terminal signs be independent fair signs. The owner
tournament D may be arbitrary and fixed. As m tends to infinity, put
`d_i = outdeg_D(i)` and `V_D = sum_i (2d_i²−d_i)`. Then
$$
   \frac{F-m!}{(m-2)!\sqrt{V_D}} \quad\Longrightarrow\quad N(0,1).
$$
This holds along **any** sequence of owner tournaments. In particular,
the fixed count is typically asymptotic to `m!`, despite the nontrivial
extremal behavior studied in Session29 (halfpath-asymptotics).

**Proof.** By the exact degree-two coefficients in the signed-fold formula, the degree-two
Walsh part divided by `−(m−2)!` is the independent sum over owners
$$
 \sum_i (S_i^2-d_i+b_iS_i), \qquad S_i=\sum_{i\to j} x_{ij}.
$$
The summands have zero mean and variances `2d_i²−d_i`. The fourth
moment of each is `O(d_i⁴)` when `d_i≥1` (expand the eighth moment of
a Rademacher sum: only tuples with all multiplicities even contribute).
For all D, the degree sum `m(m−1)/2` and Cauchy–Schwarz give
`V_D=Ω(m³)`; the sum of the fourth moments is `O(m⁵)`. Thus the
Lyapunov ratio is `O(1/m)`, proving a CLT for this part.

For the remainder use Parseval and count the squared coefficients.
Each nonterminal nontrivial path component has `2l` edges for some
`l≥1`, with its owners at alternating interior positions; a terminal
component has `2l−1` edges and one distinguished (owning) terminal.
For a given D, upper-bound the number of usable components by allowing
every required orientation. A nonterminal component on k≥3 vertices
has k!/2 undirected path orderings and squared orientation multiplier
4; together with singletons its labelled generating series is
`G(z)=z(1+z²)/(1−z²)`. A terminal component has k! ordered paths
toward the distinguished terminal, with no orientation multiplier;
its series is `J(z)=z²/(1−z²)`. A forest with c ordinary components
has squared coefficient factorial `(c!)²`; unordered components give
the denominator c!. Thus the squared remainder scale,
divided by `(m!)²`, is bounded above by the terms with r≥2 in
$$
 \sum_{r=0}^{\lfloor m/2\rfloor}
 \frac{(m-2r)!}{m!}[t^r] \left(\frac{1+t}{1-t}\right)^{m-2r+1}.
$$
In detail, multiply `sum_c c!G(z)^c` by `1+J(z)`, replace
`1+J(z)` coefficientwise by `(1+z²)/(1−z²)`, and divide the
coefficient of `z^m` by `m!`. The index r is the number of marked
owners: ordinary size has deficit `2r` relative to counting each
ordinary component as one vertex, and a terminal component is ignored
in that component count. Degree two corresponds exactly to r=1.
Coefficientwise, `(1+t)/(1−t)≤1/(1−2t)`. The r-th bound is therefore
at most `2^r/(r! (m)_r)≤(4/m)^r/r!`. Summing r≥2 gives `O(m^−2)`.
The remainder is consequently negligible in probability at the
standard deviation scale of the leading part, which is of order
`m! / sqrt(m)`. ∎

## 2. Exact second moment in the fully uniform model

For a uniformly random labelled tournament D as well, an exact formula
is available:
$$
 \frac{\mathbb E F^2}{(m!)^2}
 =\sum_{r=0}^{\lfloor m/2\rfloor}
   \frac{(m-2r)!}{m!}[t^r]
   \left(\frac{1+t/4}{1-t/4}\right)^{m-2r+1}.
$$
Indeed, the required owner edges are independently present in D with
probability 1/2. Weight each nonterminal component on `2l+1` vertices
by `2^(-2l)`, and each terminal component on `2l` by `2^(1−2l)`.
The corresponding series in the preceding forest argument become
`G(z)=z(1+z²/4)/(1−z²/4)` and `1+J(z)=G(z)/z`. The r=1 term is
`1/(2m)` exactly. Since `(1+t/4)/(1−t/4)≤1/(1−t/2)`
coefficientwise, the remaining terms are bounded by `1/(m^r r!)`.
Thus `Var(F)/(m!)²=1/(2m)+O(m^−2)` in the uniform anti-invariant
labelled ensemble at either parity. The unused mate-arc choices at
odd order do not affect F. For a fixed D, the leading conditional
variance is exactly the degree-two variance and depends on D only
through its triangle count:
`V_D = binom(m,2)+4(binom(m,3)−c_3(D))`, since each transitive triple
has exactly one vertex which dominates its other two vertices.

For uniform D, put `mu=(m−1)/2`. The deterministic sum of the
`d_i−mu` is zero, and
`E sum_i(d_i−mu)²=m(m−1)/4`, because each degree is binomial. Markov's
inequality gives `sum_i(d_i−mu)²=o_P(m³)`; hence `V_D/m³→1/2` in
probability. The uniform Lyapunov and remainder bounds from §1 give
the unconditional normalization
$$
   \sqrt{2m}\,(F/m!-1)\;\Longrightarrow\; N(0,1).
$$

Here “uniform” means uniform labelled orientations invariant under a
**specified** anti-involution, not uniform self-converse isomorphism
classes. The normal law provides no conclusion about complete
integer-valued support or the total-path conjecture.
