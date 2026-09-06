# Freezing every row label doubles the diagonal-defect exponential rate

Status: **PROVED ANALYTICALLY + INDEPENDENTLY AUDITED**; the specified exact
controls are **FINITE-EXACT**. This is a separate extension of the frozen,
audited tenth theorem. No theorem identifier or priority claim is requested.

The [independent referee](overnight11_20260906_no3line_rowfreeze_audit.md)
accepts the complete analytic proof with no correction requested. Its separate
edge-subset matching and row-mask zeta engine passes 134,492 exact gates in
normal and optimized modes. All four larger conditional mean fractions also
agree exactly with the producer's component-polynomial engine.

## 1. Inheritance and the corrected conditioning step

The closest proved mechanism is
[the tenth excess-defect theorem](overnight10_20260906_no3line_defect.md),
with its [independent audit](overnight10_20260906_no3line_exponential_audit.md).
It proves a uniform exponential probability bound by exposing two independent
label permutations and uses the exact matching polynomial of the full cycle
skeleton. Its row/column transposition range four is sharp at n=5.

The least-used sidecar is the retained conditional-label viewpoint of
[the palette theorem](overnight6_20260906_no3line_palettes.md).
Here its concrete replacement for the full matching polynomial is the family
of matching polynomials of the **row-induced prefix and suffix subgraphs**.
The concept board is: fixed row labels; induced row subgraphs; their shared
columns; bounded occupancy excess; one remaining random permutation;
the full-success versus selected-direction zero event.

The cheap hostile is that conditioning does change the exact mean. For the
ten-cycle with row i adjacent to columns i and i+1 modulo 5, row labels
rho=(0,1,2,3,4) give conditional mean 11/10, while rho=(0,2,3,1,4) gives 7/6.
The earlier unconditional mean is 17/15. Thus the old exact mean cannot be
substituted unchanged after freezing rows. This is a new conditioning
near miss, not a retraction of the tenth theorem.

The repaired claim is stronger in the useful direction: the conditional mean
has the same linear term **uniformly for every row ordering**. This permits
concentration over only the column permutation, halving the denominator.
No independence between diagonals or geometric triples is asserted.

## 2. Conditional theorem and its probability consequence

Let G be any simple bipartite 2-regular graph on n+n vertices. Fix any
deterministic row labelling rho, and choose only the column labelling sigma
uniformly. Keep the tenth notation

```text
D_d={(r,c):c-r=d},  L_d=n-|d|,  Y_d=|B intersect D_d|,
F=sum_d (Y_d-2)_+,  alpha=1-5/e^2,
mu(rho)=E_sigma F.
```

Let X count every nonaxis integer-collinear triple. As before, F=0 is
equivalent to no triple on the selected slope-one diagonals, and X=0 implies
F=0. It is not equivalent to full no-three-in-line success.

**Theorem.** For every such G and every fixed rho, if n>=4 then

```text
alpha*n-21 <= mu(rho) <= alpha*n+20,                     (1)
mu(rho) >= 2(n-9)/15 + 2/[n(n-1)].                       (2)
```

For n>=2, concentration conditional on rho gives

```text
P_sigma(X=0 | rho) <= P_sigma(F=0 | rho)
 <= exp(-mu(rho)^2/[8(n-1)]).                           (3)
```

In particular define, for n>=4,

```text
b_n^*=max(0, alpha*n-21, 2(n-9)/15+2/[n(n-1)]).
```

Then (3) is at most exp(-(b_n^*)^2/[8(n-1)]), uniformly in both G and rho.
Consequently

```text
limsup_(n->infinity) (1/n) log sup_(G,rho)
                    P_sigma(X=0 | rho)
 <= -(1-5/e^2)^2/8
 = -0.013067267481528453... .                            (4)
```

This strengthens the tenth limiting coefficient by a factor of two. It holds
even if the row ordering is chosen adversarially in advance; the remaining
column permutation must still be uniform conditional on that choice. It does
not apply to choosing rows after inspecting the random columns in order to
find a successful board.

Averaging over an independent uniform row permutation proves the same improved
rate for the original two-permutation model. For that original model one may
combine it with the tenth finite bound:

```text
P_G(X=0) <= min(exp(-n^2/[900(n-1)]),
                exp(-(b_n^*)^2/[8(n-1)])).              (5)
```

The first bound in (5) is not newly claimed for each fixed row ordering.
Uniform consequences pass to arbitrary mixtures over G and fixed row orders
when the columns remain conditionally uniform. As in the tenth theorem, a
bound around an arbitrary mixture's own mean does not follow from (3).
No extremal nonexistence or optimal exponential rate is asserted.

## 3. The exact conditional matching sidecar

For a diagonal let S be its set of row vertices, now fixed by rho, and let
L=|S|. Write G[S] for the graph retaining all edges incident to rows in S,
together with their column endpoints, and m_k(S) for its unordered k-matchings.
For 0<=k<=n,

```text
E_sigma(Y_d)_k = k! m_k(S)/(n)_k.                        (6)
```

Indeed a matching chooses k distinct rows and k distinct source columns.
There is exactly one prescribed target column label for each of these rows
on the diagonal. A uniform column permutation meets these k distinct label
requirements with probability 1/(n)_k. Ordering the k chosen cells gives k!.
For k>L the matching number is zero; for k>n the factorial moment is defined
directly to be zero without evaluating the displayed quotient.

By the tenth positive-part expansion, the exact conditional mean is

```text
mu(rho)=sum_d sum_(k=3)^(L_d)
            (-1)^(k+1)(k-2) m_k(S_d)/(n)_k.             (7)
```

When d>=0, S_d is the row prefix of length n-d; when d<0 it is the suffix
of length n+d. The full-length graph occurs once. Formula (7) explains why
the complete skeleton cycle profile alone no longer determines this mean.

The additional data are inexpensive: every nontrivial component of G[S] is
a path or an even cycle. A path of q edges has matching coefficients
binom(q+1-k,k); a cycle of q edges has coefficients
binom(q-k,k)+binom(q-k-1,k-1). Taking the product over components computes
all m_k(S). Isolated columns are irrelevant. This is an exact algorithmic
sidecar, not a hypothesis about the distribution of row orderings.

## 4. A uniform conditional factorial bound

Set theta=L/n, lambda=2theta. For every S, n>=2 and k>=2,

```text
0 <= lambda^k-E_sigma(Y_d)_k
 <= 2^k k(k-1) theta^(k-1)/n.                            (8)
```

To prove it, first take 2<=k<=L. Choose k ordered distinct rows in S,
then independently choose one of the two incident edges at each row. If
q_S(k) is the probability that their columns are distinct, then

```text
k!m_k(S)=(L)_k 2^k q_S(k).                              (9)
```

Let c(S) count columns whose two neighbor rows both lie in S. There are 2L
row-edge incidences, so c(S)<=L. A selected pair of distinct rows has column
collision probability exactly

```text
2c(S)/[4L(L-1)] = c(S)/[2L(L-1)] <= 1/[2(L-1)].         (10)
```

This formula allows two columns with the same neighbor-row pair, as happens
in C4; it counts column identities, not just distinct row pairs. A union bound
over the k positions gives

```text
1-q_S(k) <= k(k-1)/[4(L-1)].                            (11)
```

The same row-sampling collision bound as in the tenth theorem is

```text
0<=theta^k-(L)_k/(n)_k <= binom(k,2)theta^(k-1)/n.
```

Combining with (6), (9), and (11), then using
theta^k/(L-1)<=2theta^(k-1)/n for L>=2, proves (8).
For k>L the moment is zero; the row-sampling collision term already bounds
lambda^k because k(k-1)/2>=L. This also handles k>n without a zero denominator.
For L=0 the statement is direct, and for L=1 every k>=2 is in the previous
case. The first factorial moment is exactly lambda and is treated separately.

Thus the Poisson leading term is uniform over **all** row subsets, not just
over a randomly selected subset or over a bounded bank of cycle types.

## 5. Summing the conditional error

Let g(lambda)=lambda-2+(lambda+2)e^(-lambda), the Poisson expectation of
(Y-2)_+. The positive-part expansion is uniformly absolutely dominated by
sum_(k>=3)(k-2)2^k/k!, so its full expectation can be subtracted termwise.
Separate its odd and even terms and use (8). Exactly as in the tenth proof,

```text
-(2/n)lambda^2 cosh(lambda)
 <= E_sigma(Y_d-2)_+-g(lambda)
 <= (2/n)lambda^2 sinh(lambda).                         (12)
```

The trapezoid identity for the diagonal lengths is unchanged by conditioning:

```text
alpha*n <= sum_d g(2L_d/n) <= alpha*n+4e^(-2).
```

For increasing h>=0 with h(0)=0, the same sum is at most
n*integral_0^2 h+h(2). Write S=sinh(2), C=cosh(2). Since

```text
integral_0^2 lambda^2 cosh(lambda) d lambda =6S-4C,
integral_0^2 lambda^2 sinh(lambda) d lambda =6C-4S-2,
```

the lower error in (12), summed over diagonals, is at most

```text
2(6S-4C)+8C/n <=12S-6C =3e^2-9e^(-2)<21,              (13)
```

for n>=4. Including the trapezoid excess, the upper error is at most

```text
2(6C-4S-2)+8S/n+4e^(-2)
 <=12C-6S-4+4e^(-2)=3e^2+13e^(-2)-4<20.               (14)
```

The strict constants follow from the rational Taylor enclosure 7<e^2<37/5.
This proves (1). The actual expressions in (13)-(14), with their n dependence,
can be retained for sharper finite estimates.

For completeness there is also a direct fourth-order finite lower bound.
At k=3 two different column-collision events cannot occur simultaneously:
they would force all three chosen edges onto a column of degree three.
Thus the union bound is exact and, for L>=3,

```text
m_3(S)=(4/3)(L)_3-2c(S)(L-2)
       >=(2/3)L(L-2)(2L-5).                            (15)
```

For k=4, merely ignoring column collisions gives
`m_4(S)<=(2/3)(L)_4`. Divide these by (n)_3 and (n)_4, respectively,
and use f(y)>=binom(y,3)-2binom(y,4). Summing over the diagonal lengths,
with L=1,2 contributing zero to the third moment, gives

```text
sum_d E binom(Y_d,3) >=2(n-3)/3+2/[n(n-1)],
sum_d E binom(Y_d,4) <=2(2n-3)/15.
```

Their difference with coefficient two proves (2). This bound becomes positive
at n=9; for smaller sizes the trivial lower bound zero is retained. No claim
that (15) is attained simultaneously at every diagonal is needed.

## 6. Only one permutation remains to concentrate

Conditional on rho, a column transposition still moves at most four points,
and the net defect change is at most four. Expose the first n-1 values of
sigma. The completion-transposition bijection from the tenth proof gives
conditional interval length at most four at every step. Its self-contained
tilted-variance exponential lemma therefore yields

```text
E_sigma exp[-s(F-mu(rho))] <= exp[2(n-1)s^2].
```

Optimizing at s=mu(rho)/[4(n-1)] proves (3). Uniformity of (1) before any
row averaging proves (4). The improvement removes n-1 random reveals from
the concentration proof; it does not replace the sharp four-point bound or
assume that the conditional mean is constant.

## 7. Exact controls and the first conditioning hostile

The standalone
[source](../../04-computation/overnight11_20260906_no3line_rowfreeze.py) and
[output](overnight11_20260906_no3line_rowfreeze.out)
import no producer or repository code, use no floating-point gates, and keep
all checks active under optimized Python.

Two different rook engines are compared: a column-mask dynamic program
chooses at most one incident edge per selected row; a separate graph-component
traversal multiplies path and cycle matching polynomials. Their results agree
on every one of the **2,668 row subsets** of all **21 cycle profiles** with
shore size 2..8. Those same subsets check the exact double-column numerator
and every factorial deficit through order n+3, including zero/full subsets
and orders exceeding their row capacity.

The literal conditional universe has **312 row orderings** and **41,512 column
labellings**. Every row ordering is checked at n<=5 for every cycle type; at
n=6 each type uses four explicitly indexed row orderings and all 720 column
permutations. Each literal conditional mean equals the induced-rook formula.

All row orderings have the same mean for every type at n<=4. At n=5, C10
has conditional means 11/10,17/15,7/6, each on 40 of the 120 row orderings.
Thus n=5 is the first possible size at which exact row-independence fails
in this model, with explicit witnesses in Section 1. This finite minimality
claim is about the conditional mean, not about a failure of the uniform
linear asymptotic term.

Four larger exact targets avoid permutation enumeration: n=96 with 48 C4
components, and n=97 with one C194, each at the identity row order and at
rho(i)=5i modulo n. Induced-component rook polynomials compute their exact
conditional means. They lie within the proved envelope, and both the
fourth-order bound and alpha*n-21 are strictly positive there. Rational
Taylor arithmetic separately verifies the rounded constants 21 and 20.

```text
python 04-computation/overnight11_20260906_no3line_rowfreeze.py
python -O 04-computation/overnight11_20260906_no3line_rowfreeze.py

PASS 54906 always-active gates
normal and optimized output: byte-identical, LF

source SHA256:
12e962f404d849ecba8b2acb1db0dd252d276bb7df2d7dec9c5f1d70dfa360f8
output SHA256:
c48ea4c4a0755e08b69eee212390fc7ca75fb6742c60c90bc120ca60a92ca94d
```

The full proof is analytic; the finite bank checks normalization, conditioning,
and the positive target rather than extrapolating a decay law. The next open
question is whether a more precise conditional variance or exposure scheme
can improve this rate further. It must retain the prefix/suffix row data and
the coupling source; exact row-independence is already refuted.

**Filing:** root integrated these audited artifacts in the eleventh checkpoint;
reproduction commands above are relative to the repository root. Outside-worktree
locations preserve author provenance, not the present reproduction location.
