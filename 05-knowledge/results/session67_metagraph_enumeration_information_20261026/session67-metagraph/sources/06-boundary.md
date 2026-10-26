# Fixed-point-free cycle boundaries and live divisor transfer

Date: 2026-10-20. Research result, not a new canon ID.

**PROVED below:** exact elimination of expired divisibility coordinates; a
signed ordered-list transform for pair-free fixed-permutation counts; and a
compatible family of differential equations determining every finite cycle
alphabet boundary. **FINITE-EXACT:** the accompanying implementation computes
A000568 through 100 without constructing tournament isomorphism classes,
and checks the operator recurrence against independent small core data.
These formulas retain the Davis/Burnside fixed-permutation input. They are
not a Burnside-free solution, a scalar recurrence for A000568, or a Collatz
certificate.

## Inheritance and scope

The closest proved mechanisms are:

- THM-514, [unit-tail core kernel and prime-tail ladder](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md):
  separate fixed points from nontrivial odd automorphism cycles.
- THM-516, [Burnside core kernel phi reframe](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md):
  replace pairwise gcd interactions by a quadratic divisibility profile.
- [Weighted core count packet](../weighted_core_transport_20261020/core_count.md):
  the canonical homogeneous-chain quotient and
  \(T(p)=C(p_i/(1-p_i))\), \(C(p)=T(p_i/(1+p_i))\).
- The less-used sidecar THM-3079,
  [Laguerre row transform](../../../01-canon/theorems/THM-3079-laguerre-pf-row-transform-and-strict-integer-mesh-terminal-minus-one.md),
  supplies the neighboring factorial coefficient array in section 5.

The corrected near miss is to call the *full* divisibility profile a
compression: it is invertible, so no cycle-multiplicity information was
lost. The hostile control below shows why discarding a still-live prime
coordinate loses exact continuation weights. The positive replacement
for that failed shortcut is to discard a coordinate precisely when no
future operation can observe it.

Live concept board: automorphism cycles; canonical chains; future
continuation observables; gcd/divisor transforms; sparse backward arcs.
The first four give the exact transfer here. The fifth gives a finite
boundary truncation when a feedback budget is imposed.

## 1. Exact future-observable transfer

Let \(\lambda\) be an odd partition, interpreted as the cycle lengths of one
permutation. Write \(m_d\) for the multiplicity of \(d\),
\(M_r=\sum_{r\mid d}m_d\), and \(t=M_1\). Define

\[
 z_\lambda=\prod_d d^{m_d}m_d!,\qquad
 e(\lambda)=\sum_{\ell\in\lambda}\frac{\ell-1}{2}
              +\sum_{i<j}\gcd(\ell_i,\ell_j).
\]

There are \(2^{e(\lambda)}\) tournaments fixed by this permutation.
THM-516 gives

\[
 e(\lambda)=\frac12\left(\sum_{r\ {\rm odd}}\varphi(r)M_r^2-M_1\right).
 \tag{1}
\]

The complete profile contains the original multiplicities:
\(m_d=\sum_{j\ge1}\mu(j)M_{dj}\), with finite support understood.
Equation (1) alone therefore changes coordinates without reducing the state.

Instead process the allowed odd lengths in decreasing order. Just before
length \(d\), the processed lengths exceed \(d\). Retain their total mass
\(m\), and only profile entries \(M_r\) with odd \(r\le d\). The exact
interaction of one new \(d\)-cycle with the past is

\[
 J_d=\sum_{r\mid d}\varphi(r)M_r
     =\sum_{\ell\in\text{past}}\gcd(d,\ell).
\]

Adding \(k\) cycles of length \(d\) multiplies the stored coefficient by

\[
 \frac{2^{\,k(d-1)/2+d\binom{k}{2}+kJ_d}}{d^k k!},
 \tag{2}
\]

adds \(dk\) to the mass, and adds \(k\) to every \(M_r\) with \(r\mid d\).
After this step delete \(M_d\). No remaining smaller length is divisible
by \(d\), so no future interaction can use it. Sum coefficients of states
that now have the same mass and remaining profile.

This is an exact quotient of partial computations. More generally, for a
fixed future alphabet \(F\), equality of mass and all \(J_f\), \(f\in F\),
is sufficient for identical normalized continuation weights. A one-copy
continuation distinguishes unequal \(J_f\). The current absolute weight
can be aggregated because later updates are linear in that weight.
Under a finite mass bound, some continuations may be unavailable, so this
criterion is sufficient rather than claimed globally minimal.

After processing 3, only \(M_1=t\) remains. The result is exactly the
THM-514 kernel

\[
 B_{m,t}=\sum_{\substack{\lambda:\ |\lambda|=m,\ \ell(\lambda)=t\\
                         \text{all parts}\ge3}}
                   \frac{2^{e(\lambda)}}{z_\lambda}.
\]

Fixed points can then be restored by

\[
 A_n=\sum_{\substack{m\le n\\t\ge0}}B_{m,t}
       \frac{2^{\binom{n-m}{2}+(n-m)t}}{(n-m)!}.
 \tag{3}
\]

The new reduction is the *intermediate* aggregation while building this
kernel, rather than enumerating every partition and aggregating afterward.
It is variable elimination over continuation observables.

### Hostile control: an unexpired coordinate cannot simply be omitted

Take past partitions \((7,15)\) and \((9,13)\). Both have

\[
 (\text{mass},t,M_3,e)=(22,2,1,11).
\]

Adding a 5-cycle increases \(e\) by 8 and 4 respectively:
their gcd interactions are \(1+5=6\) and \(1+1=2\), and the new cycle's
internal contribution is 2. The missing coordinate is \(M_5\).
Thus mass, cycle count, one small-prime coordinate, and even the current
edge exponent do not determine the next weight.

The failure persists for any permanently fixed finite divisor bank. Let
\(L\) be its odd least common multiple, choose an odd prime \(p\nmid L\),
and choose sufficiently large odd \(a,b\equiv1\pmod L\), with
\(a\equiv0\pmod p\) and \(b\not\equiv0,2L\pmod p\).
Then \((a,b)\) and \((a+2L,b-2L)\) have the same mass and all tracked
divisibility counts, but their interactions with a new \(p\)-cycle differ.
This general argument does not additionally claim equality of their
current exponents; the explicit 5-cycle hostile does have that stronger
property.

## 2. Finite-exact state reduction through order 100

The universe is all odd cycle types of mass at most 100. Empty type is
included; parts 1 are handled only by (3). Exact rational coefficients
are used throughout.

| Just processed length | Retained states | Unmerged multiplicity states |
|---:|---:|---:|
| 99 | 2 | 2 |
| 51 | 26 | 26 |
| 31 | 143 | 153 |
| 21 | 547 | 608 |
| 15 | 1,568 | 2,104 |
| 11 | 2,568 | 6,580 |
| 9 | 3,663 | 13,568 |
| 7 | 3,486 | 32,648 |
| 5 | 2,049 | 98,937 |
| 3 | 834 | 444,793 |

The computation peaks at 3,663 retained states and uses 66,404 weighted
transitions. The final 834-state kernel recovers all \(A_0,\ldots,A_{100}\).
The 2,049-state row at length 5 recovers the earlier 3-tail residual state
count without first generating all 98,937 partial partitions.

An independent explicit-partition calculation using pairwise gcds checks
every count through 40 and the corresponding kernel entries. Recorded
A000568 values through 20 and the independently recorded \(A_{100}\) in
`04-computation/a000568_core_tail_ladder_codex.py` also agree.
The latter has 1,333 decimal digits; the SHA-256 of its decimal expansion is
`1d290d31b74fb08a9e73a6de3e842fc4a1d5bca8c1657107e7f14c0a0471e9c4`.

This establishes exact state and transition savings for this bound, not
an asymptotic complexity theorem or a runtime comparison: rational
arithmetic transitions cost more than counting states alone suggests.

## 3. Clone removal as a signed ordered-list transform

Fix distinct odd cycle lengths \(D=(d_1,\ldots,d_s)\), allowing 1 if desired.
For a multiplicity vector \(m\), let \(a_m\) count tournaments fixed by a
permutation of that type, and \(c_m\) count the fixed tournaments having
no homogeneous pair. Thus \(a_m=2^{e(m)}\). Empty type has count 1.

Define

\[
 L(m,k)=\frac{m!}{k!}\binom{m-1}{k-1}\quad(1\le k\le m),
 \qquad L(0,0)=1,
\]

with all other cases zero. The cycle-index clone-removal substitution gives

\[
 c_m=\sum_{k\le m}a_k
       \prod_i(-d_i)^{m_i-k_i}L(m_i,k_i).
 \tag{4}
\]

Its inverse is

\[
 a_m=\sum_{k\le m}c_k
       \prod_i d_i^{m_i-k_i}L(m_i,k_i).
 \tag{5}
\]

Proof: the coefficient of \(p^m\) in
\((p/(1+p))^k\) is
\((-1)^{m-k}\binom{m-1}{k-1}\). Multiply the cycle-index coefficient
by \(d^m m!\); its source normalization was \(d^k k!\). Tensor over cycle
lengths. The inverse substitution uses \(1-p\) and positive coefficients.
The factors \(L(m,k)\) count partitions into unordered ordered lists;
the additional \(d^{m-k}\) records the phase factor in this normalization.

In particular, a cycle type having no repeated length has \(c_m=a_m\):
there is no nonidentity contraction in (4). A single odd \(d\)-cycle has
\(c_1=2^{(d-1)/2}\). Repeated cycle lengths are the entire source of clone
corrections. For example, two 3-cycles have \(a_2=32\), \(c_2=20\),
giving cycle-index coefficients \(16/9\) and \(10/9\).

## 4. A closed operator system for the whole finite-alphabet boundary

Use exponential generating functions in cycle multiplicities,
with the cycle-length powers absorbed into the variables:

\[
 F_D(t)=\sum_m a_m\frac{t^m}{m!},\qquad
 G_D(t)=\sum_m c_m\frac{t^m}{m!}.
\]

Equation (4) is
\(G_D(t)=F_D(t_i/(1+d_i t_i))\).
Put

\[
 A_i=2^{(d_i-1)/2},\quad B_{ij}=2^{\gcd(d_i,d_j)},\quad
 H_i(t)_j=\frac{B_{ij}t_j}{1-d_j(B_{ij}-1)t_j}.
\]

Adding one \(d_i\)-cycle yields
\(a_{m+e_i}=A_i\prod_jB_{ij}^{m_j}a_m\). Consequently

\[
 \partial_iF_D(t)=A_iF_D(B_{ij}t_j),\qquad
 (1+d_i t_i)^2\partial_iG_D(t)=A_iG_D(H_i(t)).
 \tag{6}
\]

These identities hold as formal power series; no analytic convergence is
needed. Together with \(G_D(0)=1\), the equations for every \(i\) determine
all coefficients. For \(D\subseteq\{3,5,\ldots\}\), this computes the
fixed-point-free core boundary using cycle multiplicities only.

For an explicit integer recurrence, the coefficient with target \(m+e_i\)
is

\[
\begin{aligned}
 c_{m+e_i}={}&A_i\sum_{k\le m}c_k
   \prod_j\left[
      L(m_j,k_j)B_{ij}^{k_j}
      \bigl(d_j(B_{ij}-1)\bigr)^{m_j-k_j}\right]\\
 &-2d_i m_i c_m-d_i^2m_i(m_i-1)c_{m-e_i}.
\end{aligned}
\tag{7}
\]

The last term is omitted if \(m_i<2\). Every right-hand coefficient has
strictly smaller total multiplicity than the target. The zero conventions
for \(L\) handle absent cycle lengths.

The maps in (6) are elementary conjugates of scalar multiplication.
For \(h_{d,b}(t)=bt/[1-d(b-1)t]\),
\(h_{d,b}\circ h_{d,c}=h_{d,bc}\).
Thus all \(H_i\) commute. If
\(D_i=(1+d_i t_i)^2\partial_i\), then

\[
 D_iH_j^*=B_{ji}H_j^*D_i.
\]

Compatibility of mixed derivatives follows from
\(B_{ij}=B_{ji}\). This is a precise algebraic connection between the gcd
matrix and commuting fractional-linear changes of variables.

The previous equation involving only the fixed-point coordinate leaves
\(C(0,p_2,p_3,\ldots)\) unspecified. Equations (6) in the nontrivial
cycle coordinates supply that missing boundary. They do not eliminate
the need to enlarge \(D\) as the desired tournament order increases,
nor do they make the automorphism arithmetic disappear.

## 5. The older factorial-transform sidecar

THM-3079 section 2 uses
\(q(m,k)=m!\binom mk/k!\) and the row polynomial
\(Q_m(z)=m!L_m(-z)\), where \(L_m\) here is the Laguerre polynomial.
For \(m>0\),

\[
 L(m,k)=\frac{k}{m}q(m,k),\qquad
 \sum_kL(m,k)z^k=\frac zm Q'_m(z).
\]

Thus clone removal uses the derivative-weighted neighboring row of that
older factorial transform, followed by the sign and phase factor
\((-d)^{m-k}\). This connection recovers an existing computational object
instead of introducing an unrelated coefficient table. Positivity or
PF conclusions for the unsigned row are not automatically conclusions
about this alternating transform.

## 6. Sparse defects make the nontrivial boundary finite

Suppose a tournament has at most \(k\) backward arcs in some total order.
For any automorphism, every nontrivial cycle length is odd, say
\(\ell=2s+1\). Its induced tournament is regular of outdegree \(s\):
the cycle acts transitively on those vertices. In the restricted total
order, the last \(s\) vertices require at least \(s,s-1,\ldots,1\)
backward outgoing arcs in total. Summing over disjoint automorphism cycles,

\[
 k\ \ge\ \sum_{\ell\ {\rm nontrivial}}\frac{\ell^2-1}{8}
       \ \ge\ \frac{\#\{\text{moved vertices}\}}3.
 \tag{8}
\]

Consequently feedback budget \(k\) permits nontrivial automorphism support
of size at most \(3k\), and a length-\(d\) multiplicity at most
\(\lfloor8k/(d^2-1)\rfloor\). These are necessary support filters.
They do not say that every surviving cycle type is realized at that budget.

This couples the boundary calculation to a finite sparse-defect grammar:
for fixed \(k\), only finitely many nontrivial cycle coordinates can matter,
even as the total number of fixed vertices grows. The unrestricted
coefficients in (4)--(7) count all fixed tournaments; they must not be
substituted directly for counts satisfying a feedback budget.
Likewise a backward-arc enumerator tied to one chosen vertex ordering is
not a class invariant and cannot be used as a weighted Burnside statistic
without an additional invariant construction.

## Reproduction and remaining question

Run:

```text
python 05-knowledge/results/defect_boundary_20261020/boundary_check.py
python -O 05-knowledge/results/defect_boundary_20261020/boundary_check.py
```

The script uses explicit exception checks, not Python assertions.
Controls include 8,697 explicit odd partitions through mass 40; 364
inverse-transform entries; 322 operator states with 604 independent
differentiation-direction checks; rational fractional-linear composition;
and 18 actual small-core cycle types through order 7. The last control
imports the frozen previous core census solely for independent verification;
the count through 100 constructs no individual tournament cores.
Normal and optimized runs produce identical output.

The remaining compression question is concrete: can the live vector of
continuation interactions itself be replaced by a smaller exact statistic
for a constrained future alphabet? A fixed finite bank cannot handle
unrestricted new primes, as section 1 shows. A finite feedback budget does
bound the alphabet by (8), and therefore supplies the extra hypothesis
under which a finite boundary operator is mathematically justified.
