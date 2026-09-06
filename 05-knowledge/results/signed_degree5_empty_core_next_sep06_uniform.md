# Sharp second-coefficient duplication margins in every degree

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
Separate structural addendum to the frozen, independently audited
[degree-five theorem](signed_degree5_empty_core_next_sep06.md).
No theorem ID or external-priority claim.

## 1. Exact finite-degree and uniform statements

For every integer `n>=4`, define

```text
t_n=[2+(n-4)sqrt((n-1)/(2(n-2)))]/n,
c_n=(5-3t_n)/(3(1+t_n)).                                  (1)
```

Let `r_1,...,r_n` be real, with `e_2(r)=0` and
`E=e_2(r_1^2,...,r_n^2)>0`. Zero entries are allowed. Put

```text
G(s)=prod_i(1+r_i s),        D=-[s^4]G(s)^2.
```

Then

```text
D >= c_n E.                                               (2)
```

This constant is sharp for each `n`. Equality holds exactly at permutations
and common nonzero real scalings of

```text
(1,1,-q_n,...,-q_n),
q_n=1/[n-2+sqrt((n-2)(n-1)/2)],                            (3)
```

with `n-2` negative entries. In particular equality has exact degree `n`.
The recovered constants are `c_4=7/9`,
`c_5=(81-8sqrt(6))/87`, and `c_10=61/99`.

For all finite real lists satisfying these cancellation and positive-energy
hypotheses, including `n=3`, there is the stronger uniform strict inequality

```text
D/E > c_*=(13-8sqrt(2))/3 >1/3.                           (4)
```

The constant `c_*` is a sharp infimum over degrees, approached by (3), and
is never attained by a finite positive-energy list. At `n=3`, every such
list has `D/E=1`; with at most two entries there is no positive-energy list
satisfying `e_2=0`.

For a real-rooted ordinary core with a nonzero constant coefficient, multiply
`D,E` by the square of that constant. The same statements then concern its
actual ordinary-square coefficient. The monomial-shift and complex-gauge
normalization from the degree-five note remains valid, with effective index
two; no real sign is assigned before choosing a real gauge.

## 2. Inheritance and a checked primary-source antecedent

The source mechanism is the full signed-subset SOS of
[THM-4440 — signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md).
Its terminal level gives `1/3` at index two. The
[quartic stability note](signed_duplication_stability_empty_core_sep06.md)
and degree-five note recover stronger constants by retaining root geometry.
The present addendum evaluates the finite multiplicity reduction in arbitrary
degree; it does not reuse an arbitrary-row or channel-polynomial positivity
claim.

The low-value reduction has a directly applicable literature antecedent.
[Riener, *On the degree and half degree principle for symmetric polynomials*,
arXiv:1001.4464v2, Lemma 4.2](https://arxiv.org/pdf/1001.4464), building on
Timofte's principles, minimizes a linear coefficient functional over
real-rooted monic polynomials with the first `s` coefficients fixed at a
polynomial with at most `s` distinct roots. Here use
`Q(z)=prod(z-r_i)`, fix `e_1=1,e_2=0`, take `s=2`, and observe that our
objective is `12[(c-1)e_3-(c+1)e_4]`. Thus its hypotheses match. The
degree-five note supplies a self-contained tangent proof and the stronger
all-local-minima statement used for equality. This reduction is credited
as recovered mathematics, not claimed as a new general principle.

Bounded searches for the exact uniform constant and the normalized
third/fourth-power-sum inequality below did not recover a checked matching
formula. This is a retrieval report, not a claim of absence or priority.

The retained board is: actual square coefficient; complete root-product
energy; compact moment sphere; integer multiplicity; zero-energy boundary;
finite-root-plus-exponential limit. The map to two values preserves the
minimum and retains their integer multiplicities. Discarding those integers
would admit a different optimization problem, including an invalid approach
to the zero-energy one-root orbit.

## 3. Eliminate the two root values before optimizing multiplicity

As in the degree-five proof, normalize by the nonzero root sum. The list
then satisfies `p_1=p_2=1`, and

```text
E=(1-p_4)/2,
D=(5-8p_3+3p_4)/6,
F_c=6(D-cE)=3(1+c)p_4-8p_3+5-3c.                         (5)
```

For every `c>-1`, the proved minimum lemma says that every local minimum
of `F_c` on this smooth compact sphere has exactly two distinct values.
Zero entries remain in this optimization. The zero-energy orbit consists
of permutations of `(1,0,...,0)` and has `F_c=0` for every `c`.

Suppose the two values are `x,y`, with multiplicities `m,n-m`. Set `t=x+y`.
Since each root satisfies `r_i^2=t r_i-xy`, summing gives

```text
xy=(t-1)/n,
p_3=((n-1)t+1)/n,
p_4=((n-1)t^2+1)/n.                                     (6)
```

At positive energy `|t|<1`; substituting (6) into (5) gives the especially
simple exact ratio

```text
D/E=R(t)=(5-3t)/(3(1+t)),       R'(t)=-8/[3(1+t)^2]<0.     (7)
```

Thus only the largest eligible two-value sum matters. Interchange the two
values to assume `1<=m<=floor(n/2)`. The two constraint solutions give

```text
t_(m,+/-)=[2 +/- (n-2m)sqrt((n-1)/(m(n-m)))]/n.            (8)
```

The `m=1,+` choice has values `1,0`, hence zero energy, and is excluded
from the ratio but retained in `F_c`. The `m=1,-` choice has
`t=(4-n)/n<=0`. All other minus choices have `t<=2/n`.
For plus choices with `m>=2`, the square

```text
(n-2m)^2/[m(n-m)] = n^2/[m(n-m)]-4                        (9)
```

decreases strictly as `m` grows through `2,...,floor(n/2)`, because
`m(n-m)` increases. Therefore the largest eligible value is exactly
`t_(2,+)=t_n`. At `n=4` the two signs give the same unordered pair; this
introduces no extra equality shape.

Since `R` decreases, every positive-energy two-value orbit satisfies
`D/E>=c_n`. The zero-energy orbit has `F_(c_n)=0`. The compact minimum
lemma therefore proves `F_(c_n)>=0` on the entire constraint sphere,
establishing (2) without minimizing a singular quotient.

If equality holds with positive energy, it is a global minimum of
`F_(c_n)`. The all-local-minima statement forces two values, and the strict
multiplicity comparison forces `m=2,+`. Undoing the normalization gives
(3): solving `e_2(1,1,-q,...,-q)=0` chooses the smaller positive root

```text
1-2(n-2)q+binom(n-2,2)q^2=0,
q=1/[n-2+sqrt((n-2)(n-1)/2)].                             (10)
```

The larger quadratic solution gives the other branch and is not an extra
equality case, except for the already identified unordered `n=4` symmetry.

## 4. Strict uniform bound and its equivalent power-sum inequality

For `n>=4`,

```text
sqrt(2)n t_n=2sqrt(2)+(n-4)sqrt(1+1/(n-2)).
```

The elementary estimate `sqrt(1+u)<=1+u/2` for `u>=0` gives

```text
sqrt(2)n t_n
 <=2sqrt(2)+n-4+(n-4)/(2(n-2))
 <n+2sqrt(2)-7/2<n.                                     (11)
```

The signs are explicit: `n-4>=0` and `7/2>2sqrt(2)`. Hence
`t_n<1/sqrt(2)`, so (7) gives `c_n>R(1/sqrt(2))=c_*`.
Also `t_n->1/sqrt(2)`, proving `c_n->c_*` and sharpness of the infimum.
For `n=3`, `e_4=0` in the exact identities
`D=-2e_1e_3-2e_4`, `E=-2e_1e_3+2e_4`, so `D=E`. The smaller cases
have zero energy under cancellation. This completes all finite degrees.

Equivalently, for every finite real list with `p_1=p_2=1`,

```text
p_3 <=(2-sqrt(2))p_4+sqrt(2)-1.                           (12)
```

Equality in (12) holds exactly at a permutation of `(1,0,...,0)`.
Indeed its difference is `F_(c_*)/8`; (4) makes it strictly positive at
positive energy, while the zero-energy orbit gives equality. The constraint
`p_1=p_2=1` is part of this statement; no unspecified probability-measure
or weighted-multiplicity generalization is used.

## 5. The sharp limit keeps a linear exponential tail

Put `L=n-2` and use the sharp carriers

```text
G_n(s)=(1+s)^2(1-q_n s)^L.
```

Here `L q_n->2-sqrt(2)` and `L q_n^2->0`. Consequently, coefficient by
coefficient,

```text
G_n(s)->(1+s)^2 exp[-(2-sqrt(2))s].                        (13)
```

For example this follows by expanding the logarithm of the last factor
through any fixed order: only its linear coefficient survives. Its second
coefficient is zero, because (10) is exact for every `n`. The root-product
energy satisfies

```text
E_n=1+2L q_n^2+binom(L,2)q_n^4 ->1.                       (14)
```

The negative fourth coefficient of the square of (13) is exactly `c_*`,
which is also independently checked by finite formal multiplication in
`Q(sqrt(2))`.

This is a different limit from the inherited Hermite approximation. Two
reciprocal roots stay macroscopic; the remaining roots go to zero, retaining
a nonzero sum but vanishing square sum. Simply deleting the small roots
would give `(1+s)^2`, whose second coefficient is one and does not cancel.
The exponential factor records precisely the first-moment information that
such a deletion loses. In the ordinary-root picture the small reciprocal
roots escape to infinity. No entire-function or arbitrary Laurent claim is
needed to prove (2); (13) only explains the sharp uniform limit.

## 6. Exact controls and audit manifest

[Source](../../04-computation/signed_degree5_empty_core_next_sep06_uniform.py)
and [output](signed_degree5_empty_core_next_sep06_uniform.out).

```bash
python3 -B 04-computation/signed_degree5_empty_core_next_sep06_uniform.py
python3 -B -O 04-computation/signed_degree5_empty_core_next_sep06_uniform.py
```

The exact universe is all 68 oriented two-value configurations for
`n=4,...,12`, `1<=m<=floor(n/2)`, both signs, plus the nine independently
multiplied sharp carriers (3). The checker uses rational quadratic fields,
retains square-radicand and zero-energy cases, verifies every power-sum
elimination, checks the rational squared multiplicity comparison and the
strict uniform root-sum bound, recovers `n=4,5`, and verifies the formal
linear-drift limit and the separate `n=3` case. All **519 explicit gates**
pass. The proof is the compact minimum and exact multiplicity argument;
the small finite controls are not an all-degree census.

Frozen normal/optimized equality and raw hashes are recorded below.

```text
source SHA256 733478d5ec2f3289affeb894448fa9b98131c7a0b29bff979065bae9ec1df867
output SHA256 9385ca3e4455a9887495ec7f44e2b4b78f92471431c4802ca9ab2e58d5444f14
trace  SHA256 954ef412205a73dcb4efb38af18a3906764df33e66ce8788ea1bbbab23c70f62
```

Independent final referee `certificate_audit`: **PASS**. The referee
re-derived the two-value power sums, fractional-linear ratio, integer
multiplicity maximizer, zero-energy and degree-four symmetry cases,
equality shape, strict uniform bound, coefficientwise exponential limit,
and equivalent power-sum inequality. They also checked the cited primary
Riener lemma and its exact fixed-coefficient map, replayed both this source
and the degree-five source, and verified all four raw hashes. No mathematical
correction was required; source/output remain frozen.
