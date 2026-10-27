# Parity halves the cyclic-context observability budget

## Result

For a nonempty tournament $T$ of **known order** $n$, put
\[
\Psi_T(q)=\sum_{d=1}^{n}\operatorname{pc}_T(d)(q)_d,
\qquad(q)_d=q(q-1)\cdots(q-d+1).
\]
The path-cover convention is unordered spanning directed paths, with
singletons allowed. Let $\Phi(T)=C_3[T,K_1,K_1]$, oriented
$T\to a\to b\to T$, and $H_k=H(\Phi^k(T))$.

**Theorem.** With $r=\lfloor(n-1)/2\rfloor$, the $r$ observations
\[
H_0,H_1,\ldots,H_{r-1}
\]
recover the complete path-cover profile. For $n=1,2$, no observations
are needed: the known order determines the polynomial and profile.
For $n\ge3$, the largest probe has $n+2r-2\le2n-3$ vertices.
This is exact reconstruction from counting observations, not an efficient
algorithm for obtaining those observations.

The polynomial framework has classical context: Grinberg and Stanley's
[Rédei–Berge symmetric function paper](https://arxiv.org/abs/2307.05569)
studies weakly increasing colorings with strict increases on designated adjacent
arcs. Their convention marks forward arcs; the formula here marks backward
arcs and therefore specializes their function for the opposite tournament.
Reversing each vertex word reconciles the principal specialization. No
claim of novelty for polynomial parity or of global priority is made.

## 1. Counting polynomial, monicity, and parity

For a permutation $w=(w_1,\ldots,w_n)$, let $b_T(w)$ be the number of
backward adjacent arcs: $w_{i+1}\to w_i$. Assign colors in
$\{1,\ldots,q\}$ weakly increasingly along $w$, requiring strict
increases at every backward arc. Each equal-color run is then a directed
path. Conversely, an unordered path cover and an injection of its $d$
components into the $q$ colors give exactly one such colored word, by
ordering the components by color. Thus, for positive integral $q$,
\[
\Psi_T(q)=\sum_{w\in S_n}\binom{q+n-1-b_T(w)}{n}.              \tag{1}
\]
This is also a polynomial identity, taking binomial coefficients in their
polynomial sense. Since the all-singleton cover is unique, $\Psi_T$
is monic of degree $n$; it has zero constant term.

Reversing a word complements its backward-arc count:
$b_T(w^{\mathrm{rev}})=n-1-b_T(w)$. The identity
\[
\binom{-q+n-1-b}{n}=(-1)^n\binom{q+b}{n}
\]
in (1) therefore gives
\[
\boxed{\Psi_T(-q)=(-1)^n\Psi_T(q).}                            \tag{2}
\]
Consequently the possible polynomials lie in the affine space
\[
q^n+\operatorname{span}\{q^{\epsilon},q^{\epsilon+2},\ldots,q^{n-2}\},
\quad \epsilon=1\text{ for odd }n,\quad\epsilon=2\text{ for even }n,
                                                                  \tag{3}
\]
of dimension $r$.

## 2. Cyclic attachment becomes a finite-difference operator

The exact substitution kernel gives, with zero padding,
\[
p'_d=p_{d-2}+(2d-1)p_{d-1}+d(d+2)p_d+d(d+1)p_{d+1}.           \tag{4}
\]
For completeness, (4) follows by taking $[x^ayz]$ in
\[
W^d,\qquad W=\frac{x+y+z+xy+yz+zx+3xyz}{1-xyz},
\]
and multiplying by $a!/d!$. The needed coefficient identity is
\[
[yz]W^d=d x^{d-1}(1+3x+x^2)+d(d-1)x^{d-2}(1+x)^2;
\]
for $d=1$ the second term is absent. Actual distinct path components
decorate quotient occurrences before dividing by $d!$; bare words do
not justify that division. See the inherited
[THM-3121](https://github.com/eliottcassidy2000/math/blob/02d1218a9a48986564693cae1ae65a050fb9f0e6/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).

Substituting (4) into the falling-factorial expansion gives
\[
\boxed{\Psi_{\Phi(T)}(q)=\mathcal L\Psi_T(q),\qquad
\mathcal Lf(q)=q^2f(q)+q\bigl(f(q+1)-f(q-1)\bigr).}           \tag{5}
\]
One can check this coefficientwise: the contribution of $(q)_a$ to
$\mathcal Lf$ is
\[
(q)_{a+2}+(2a+1)(q)_{a+1}+a(a+2)(q)_a+a(a-1)(q)_{a-1}.
\]
The last term is zero for $a=1$.

Since $\Psi_T(1)=p_1=H(T)$, we have
$H_k=(\mathcal L^k\Psi_T)(1)$. Starting at evaluation point 1, the
operator accesses only evaluations $0,1,\ldots,k+1$, and the coefficient
of $\Psi_T(k+1)$ is
\[
1\cdot2\cdots k=k!.
\]
Indeed only consecutive upward shifts reach $k+1$ in $k$ steps.
The value at zero is known to vanish, and the zero row cannot propagate
into negative arguments. Triangular exact division therefore recovers
$\Psi_T(1),\ldots,\Psi_T(r)$ from the $r$ observations.

Subtract $q^n$, divide by $q^{\epsilon}$, and interpolate at
$q^2=1^2,2^2,\ldots,r^2$. These are distinct nonzero points, and the
unknown polynomial in $q^2$ has degree at most $r-1$. Vandermonde
invertibility proves unique recovery. Finally the falling factorials are
a monic triangular basis, recovering all $p_d$.

## 3. The affine dimension is attained by tournaments

Let $B_s$, $s\ge3$, be the increasing transitive tournament on
$0,\ldots,s-1$, except that its long edge is reversed to $s-1\to0$.
Then
\[
\boxed{\Psi_{B_s}(q)=q^s+q\bigl((q+1)^{s-2}-(q-1)^{s-2}\bigr).} \tag{6}
\]
To prove (6), interpret $\Psi$ as the sum over all vertex colorings
of the product of Hamiltonian counts within their color classes. This is
exactly the colored-cover interpretation above. Normally every class is
transitive and contributes 1. Only a class containing both endpoints and
$t\ge1$ intermediate vertices changes the count. Its Hamiltonian paths
are the increasing path, plus $2^t$ paths using the reversed edge: each
intermediate vertex is assigned before $s-1$ or after 0, in increasing
order within each segment. With $t=0$ there is just one path.
Choose the endpoints' common color, the $t$ intermediates in that color,
and the other colors. The excess is
\[
q\sum_{t=1}^{s-2}\binom{s-2}{t}2^t(q-1)^{s-2-t},
\]
which is (6).

Joining $B_s$ transitively to $n-s$ singletons multiplies its
polynomial by $q^{n-s}$: for every color the factor-paths concatenate
uniquely across a transitive join. For $s=2j+1$, the perturbation from
$q^n$ has lowest term $2q^{n-2j}$. Taking $j=1,\ldots,r$ supplies
triangularly independent perturbations, plus the transitive tournament's
polynomial $q^n$. Thus the affine span of actual tournament polynomials
has dimension **exactly $r$**.

Accordingly $r$ real-valued linear measurements are necessary for injective
linear recovery on this full real affine span, and the probes above attain
that budget. This is not a lower bound for arbitrary nonlinear encoding
of the finite set of tournament profiles. Fewer measurements might separate
some finite families; such a claim requires a different argument.

## 4. Small case and executable checks

Orders 3 and 4 have $r=1$, explaining why $H$ alone determines the
profile there. At order 5, $r=2$: the earlier codes 69 and 76 have
polynomials respectively
\[
q^5+8q^3+6q,\qquad q^5+10q^3+4q.
\]

Both give $H_0=15$, but their $H_1$ values are 123 and 135.
For the general order-five polynomial $q^5+a q^3+bq$,
\[
H_0=1+a+b,\qquad H_1=33+9a+3b,
\]
so $a=(H_1-3H_0-30)/6$ and $b=H_0-1-a$. This recovers
$(a,b)=(8,6)$ and $(10,4)$ exactly.

Run `python check_parity.py` and `python -O check_parity.py`.
The standalone standard-library checker verifies 77 actual tournaments,
the descent formula at four positive integers for each, parity, the
operator identity, and reconstruction from the reduced number of probes.
It checks 15 reversed-long-edge or joined-basis cases through order 8,
including nonsingularity of the actual tournament affine bases. Ten abstract
falling-factorial basis tests verify the operator identity beyond realizable
profiles, and the converse pair below is checked explicitly. It rejects
insufficient observations, empty order, and fractional interpolation.
All checks use explicit exceptions and remain active under optimization;
normal and optimized output agree byte-for-byte.

## Boundary: a count carrier, not a tournament classifier

Together with THM-3121, equal profiles are indistinguishable by raw-H
observations in nested **disjoint uniform module-substitution contexts**:
there is one distinguished tournament slot, the other factors are fixed,
and interblock arcs are uniform. The cyclic probes prove the converse at
known order. This does not cover overlapping blocks, shared ports, or
nonuniform attachments.

Reversing every arc reverses every path and preserves the profile.
For instance a source joined to a directed triangle, and its converse
(a triangle joined to a sink), are nonisomorphic tournaments with the same
profile $(3,9,6,1)$ and polynomial $q^4+2q^2$. The first has a source
and no sink; the second a sink and no source. No uniform-substitution raw-H
context distinguishes them. Thus this carrier does not classify tournament
isomorphism or enumerate A000568; it is not the full carrier required for
more general guarded or shared-port attachments. Counts here are raw
Hamiltonian counts, not counts divided by automorphism orders.


## Exact integers matter: an all-depth modulus-12 obstruction

For the two five-vertex witnesses, the polynomial difference is
\[
\Psi_B(q)-\Psi_A(q)=2(q^3-q)=2(q-1)q(q+1).
\]
Three consecutive integers have product divisible by 6, so this difference
is divisible by 12 at every integer argument, positive, zero, or negative.
The operator $\mathcal L$ takes integer linear combinations of function
values at $q-1,q,q+1$ and therefore preserves the set of integer-valued
functions whose every value is divisible by 12. Induction proves
\[
H_k(B)\equiv H_k(A)\pmod {12}\qquad\text{for every }k\ge0.
\]
Nevertheless the profiles differ and $H_1(B)-H_1(A)=12$. Thus even infinitely
many of these probes modulo 12 cannot replace exact-integer observations.
This is an all-depth proof, not extrapolation from a finite test. The checker
verifies the symbolic initial difference and illustrates twelve iterations;
its interpolation and division guards deliberately retain exact arithmetic.
