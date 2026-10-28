# Exact future responses recover incidence states

**Status: PROVED (elementary, all finite alphabets); FINITE-EXACT audit below.**
Derived in session 88 from the June Burnside/divisor bridge, with the July
conditional-expectation and continuation-state mechanisms. No priority claim.

## 1. A general incidence-algebra theorem

Let `P` be a finite poset, let `w_d>0` for every `d in P`, and let
`m=(m_a)_(a in P)` be a vector of nonnegative integer atom multiplicities.
Define

```
Z_(d,a) = 1_(d<=a),
M=Zm,
C=Z^T diag(w) Z,
Q(m)=1/2 m^T C m + ell^T m.
```
Here `ell` is any fixed real vector. The operations are adding one atom:
`m -> m+e_a`. The corrected response is

```
R_a(m)=Q(m+e_a)-Q(m)-Q(e_a)=(Cm)_a.
```

**Theorem.** The full one-step response vector `R(m)` determines `m` exactly.
Consequently, if a quotient state determines the present value `Q(m)` and the
value of `Q` after every allowed addition (and therefore after every future
addition word), that quotient must distinguish every different multiplicity
vector. There is no nontrivial exact behavioural quotient for this particular
observable and this full operation alphabet.

**Proof.** In a linear extension of the poset, `Z` is triangular with diagonal
ones, hence invertible. For `v!=0`,

```
v^T C v = sum_d w_d (Zv)_d^2 >0.
```

Thus `C` is invertible and `m=C^(-1)R(m)`. More explicitly, `R=Z^T diag(w)M`;
apply poset Möbius inversion to the upper incidence transform, divide by the
positive weights, and apply inverse `Z` to obtain `m`. If two quotient-equal
states have the same current and all one-step values, their responses agree,
so their multiplicities agree. Conversely, the multiplicities update by
addition and retain all values. QED.

**Equivariant symmetry refinement.** The theorem compares fixed, individually
labelled probes. It does not prohibit genuine relabelling symmetries. If a
finite group `G` acts by poset automorphisms preserving the weights (and, to
preserve `Q`, preserving `ell`), then `C` commutes with the induced permutation
action and

```
R(gm)=gR(m),
Stab_G(R(m))=Stab_G(m).
```

The first identity follows by equivariance of the incidence relation and
weights; the second uses injectivity in both directions. Thus the response
map is a faithful equivariant embedding, and `G`-orbits of multiplicity states
correspond exactly to `G`-orbits of response vectors. Quotienting *after*
retaining the coupled probe labels is lawful; forgetting their coupling
before comparing states need not be. This is a positive symmetry statement,
not merely a no-compression statement.

**Boundary.** A restricted operation menu need not separate states. If only
responses in a subset `A` are demanded, the linear invisible directions are
exactly `ker(C_(A,*))`. A fixed yes/no target may require vastly less than all
exact future values. The theorem is not an information-theoretic lower bound
for every enumeration algorithm, a prohibition on aggregate recurrences, or
a claim that a convenient state needs a long written certificate.

## 2. Burnside odd-cycle cores: a faithful arithmetic future state

The inherited source is
`01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md` (June 15).
A core `mu` is a multiset of odd integers at least three: **a permutation
cycle type**, not a tournament isomorphism class. Write `m_p` for the number
of parts `p`, `t=sum_p m_p`, and

```
M_d=sum_(p:d|p) m_p.
e(mu)=sum_i (p_i-1)/2 + sum_(i<j) gcd(p_i,p_j).
```

This is the number of unordered-edge orbits of a permutation with cycle
type `mu`. Every such odd-order permutation fixes exactly `2^e(mu)` labelled
tournaments. Its contribution to the *unlabelled* Burnside sum additionally
has denominator `z_mu=product_p p^m_p m_p!`; that denominator is not being
silently discarded here.

Let the allowed part alphabet be the odd numbers in `[3,N]`. No total-mass
bound is imposed. Set `C_(p,q)=gcd(p,q)`. The identity

```
gcd(p,q)=sum_(d|p,d|q) phi(d)
```

writes `C` as the Gram matrix of divisor indicators with positive Euler-phi
weights. The columns indexed by `d>=3` form a triangular incidence matrix
with diagonal one; the `d=1` column adds a positive rank-one term. Thus `C`
is positive definite even though the part `1` is excluded.

Direct expansion gives

```
e(mu)=1/2 m^T C m - t/2,
Delta_p e = (p-1)/2 + sum_q gcd(p,q)m_q.
```

**Corollary.** The one-step responses
`Delta_p e-(p-1)/2`, for every odd `3<=p<=N`, determine the entire core.
Two distinct cores cannot remain indistinguishable for all allowed future
add-part operations while preserving exact `e`. The current `e`, mass and
number of parts do not change that conclusion.

There is also a transparent divisor-by-divisor decoder when `t` is retained.
Set `A_1=0` and, for odd `p>=3`,

```
A_p=Delta_p e-(p-1)/2-t
   =sum_(d|p,d>=3) phi(d)M_d.
M_d = 1/phi(d) sum_(k|d) mobius(d/k) A_k       (d>=3),
m_j = sum_(k odd,jk<=N) mobius(k) M_(jk)     (j>=3).
```

The first inversion recovers divisor counts from the *lower-divisor*
response sum; the second inverts the *upper-multiple* counts. Confusing these
two directions is an easy wrong formula. This is the precise bridge from
Burnside cycle-index symmetry to incidence algebra, not a numerical analogy.

## 3. Smallest hostile: same current exponent, different future response

The two cores

```
mu=(3,11),   nu=(5,9)
```

have identical `(mass,t,e)=(14,2,7)`, hence each corresponding permutation
fixes `128` tournaments. Their `z_mu` values are `33` and `45`, so their
weighted Burnside summands are **not** equal.

Their first available probe, adding a 3-cycle, gives the same increment `5`.
The next probe, adding a 5-cycle, gives

```
Delta_5 e(mu)=2+gcd(5,3)+gcd(5,11)=4,
Delta_5 e(nu)=2+gcd(5,5)+gcd(5,9)=8.
```

The resulting exponents are `11` and `15`, and the fixed-tournament counts
are `2^11=2048` and `2^15=32768`. The missed coordinate is visible directly:
`M_5(mu)=0`, `M_5(nu)=1`. Agreement on the 3-cycle probe did not authorize a
5-cycle update.

**Minimality is elementary.** At length one, mass determines the core. At
length two below mass 14, the only masses with two possibilities are 10 and
12: `(3,7)` versus `(5,5)` have gcds `1,5`; `(3,9)` versus `(5,7)` have gcds
`3,1`. Their exponents differ. At length three below 14, only mass 13 has two
possibilities: `(3,3,7)` and `(3,5,5)`, with exponents `10,12`. At length four
below 14, only `(3,3,3,3)` occurs; larger lengths cannot occur. Thus mass 14
is the first coarse collision, not merely the first in a tested prefix.
The numerical equality of this mass with another problem's parameter 14
has no claimed mathematical significance.

## 3A. Restoring the unit part changes the boundary

The minimality claim above is **only for unit-free cores**. If ordinary odd
cycle types include the part `1`, the first collision is already

```
(1,7), (3,5): (mass,t,e)=(8,2,4).
```

Adding `1` gives increment `2` in both; adding `3` gives increments `3,5`.
The separate exhaustive unit-containing audit through mass 20 confirms this
minimal boundary. It must not be mixed with the unit-free mass-14 claim.

The response-reconstruction scopes are:

1. **No units, support bounded by N:** queries at every odd `3<=p<=N`
   alone determine all multiplicities by the positive-definite restricted
   GCD matrix. The length `t` need not be an extra input.
2. **Units allowed, support bounded by N, t known:** the same queries and
   the displayed Möbius decoder recover all `m_j` for `j>=3`; then
   `m_1=t-sum_(j>=3)m_j`.
3. **Units allowed and query 1 permitted:** `Delta_1 e=t`, so case 2 applies.
4. **Units allowed, no t supplied, all odd queries >=3 permitted:** choose
   an odd prime larger than every part; its corrected response is `t`, so
   case 2 again applies. For comparing two finite-support cores, choose a
   prime beyond the union of their supports.
5. **Units allowed, no t supplied, only queries 3..N:** the *response vector
   alone* may not suffice. Already with `N=3`, `(1,1,1)` and `(3)` have
   the same mass 3 and corrected add-3 response 3, but lengths `3,1`
   and current exponents `3,1`. With `N=5`, `(1,1,1,1,1,1,1)` and `(3,3,5)` both
   have corrected responses `(7,7)` at queries `(3,5)`. Their lengths are
   `7,3` and their present exponents are `21,9`. This is not a counterexample
   to reconstruction from responses **plus** present exponent; no such
   stronger failure is asserted.

## 4. A rank-one conditional-expectation obstruction

For a finite fibre `B` carrying positive probabilities `pi_x`, the averaging
operator is the rank-one projection

```
(E_B f)(x)=sum_(y in B) pi_y f(y).
```

For any operation-response observable `r`, its loss on this quotient is

```
||r-E_B r||_pi^2
 = 1/2 sum_(x,y in B) pi_x pi_y (r(x)-r(y))^2.
```

It vanishes exactly when that response descends to one value on the fibre.
This is the finite form of the conditional-variance/complete-graph Dirichlet
identity in `THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary.md`.
With equal weights on the hostile pair and `r=Delta_5 e=(4,8)`, the mean is
`6` and the loss is exactly `4`. The 3-cycle probe has loss zero. Thus a
quotient can be perfectly faithful for one operation and fail for the next;
the operation alphabet is part of the symmetry claim.

For a map `f_a:X->X'`, coarse partitions on `X,X'`, conditional expectations
`E,E'`, and pullback `U_a h=h o f_a`, exact descent of the operation is

```
(I-E) U_a E'=0.
```

Indeed, `E'` ranges over all target-block-constant functions. Its pullbacks
are source-block constant iff every pair in one source block reaches one
target block. Indicator functions of target blocks supply the converse.
This is a general criterion; it does not assert that a particular root-switch
quotient satisfies it. A short certificate for one target is weaker than
commutation with every allowed future operation.

## 5. Exact audit and reproducibility

```
python3 divisor_future_state_audit.py
python3 -O divisor_future_state_audit.py
```

Universe: all odd-part-at-least-three cores of total mass at most 40, empty
core included. Results:

- 1,113 cores; direct permutation/edge-orbit counts independently match `e`.
- 21,147 add-part checks for all odd `3<=p<=39`.
- 1,113 divisor and multiplicity reconstructions, with exact integers.
- 225 `(mass,t,e)` collision cells; all 1,419 collision pairs separated by an
  allowed one-step probe.
- The minimal hostile and its exact conditional-variance loss are checked.
- A separate unit-containing universe through mass 20 checks the mass-8
  boundary and the finite-query/no-length response-only hostile.
- Normal and optimized executions produce byte-identical JSON.

The proof is all-alphabet/all-mass; these computations are finite checks of
its specialization, not a replacement for the proof.
