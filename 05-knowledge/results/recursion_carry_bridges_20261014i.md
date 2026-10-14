# A triangular quotient of Mahler's ceiling map, with its lost sheet restored

2026-10-14. **PROVED:** the quotient, exact precision law, carry derivative,
marked restoration, and period-doubling statements below. **FINITE-EXACT:**
the stated rational, residue and word controls. Earlier divisor and clock
results are **INHERITED**. Mahler's Z-number problem and Collatz remain **OPEN**.

## 1. Recover the objects before joining their numbers

The closest mechanism is the source-retaining carry coordinate of
[THM-2228, Mahler carries and integral stabilization](../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md).
The hostile is a formal safe rational source which is not an ordinary positive
integer. The corrected near miss is a common finite clock read as a common
evaluator. The least-used sidecar here is the **reflection sheet**. Our board
is divisor counts / full recurrence state / quadratic observation / carry
orientation / precision / strict real safety.

The older references have several distinct meanings:

| Object | Exact meaning and retained boundary |
|---|---|
| `F=U+S` | For an integer N>=2, F counts proper nontrivial divisors, U counts prime divisors strictly below N, and S counts squarefree proper nontrivial divisors. U's underlying set lies inside S's. This is a count equality, not a disjoint partition or a recurrence. |
| `10=7+3` | At N=p^2qr, the counts are F=10,S=7,U=3. The full balanced shapes are p,p^3,p^2qr, giving F=0,2,10. The 1 at p^3 is U=S=1. Thus a list `{10,2,1}` mixes component roles unless another definition is specified. |
| Squaring a marked prime factor | p^2qr to p^4qr changes the counts to18,7,3, a defect of8. The valid18=15+3 repair changes S to fourth-power-free divisors. |
| Fibonacci and triangular clocks | Full Fibonacci states modulo10 and triangular states modulo30 are pointed60-cycles. Their scalar pair has only48 values; adding index parity restores60. Triangular residues modulo8 have period16, not8. |
| Triangular Fibonacci values | The old bounded calculation and the present independent check through Fibonacci index40 give positive common values1,3,21,55. The family(2n-1)(4n-1)=T_(4n-2) is triangular for every integer n, but its next value105 is not Fibonacci. No all-index classification is proved here. |

The divisor definitions and complete classification are retained in
[predecessor prime shapes, §1](predecessor_prime_shapes_20261014b.md), tracing
back to [divisor balance DB1--DB3](arithmetic_braids_20260917_divisors.md).
The clock theorem is
[THM-4035, sixty-clock separation](../../01-canon/theorems/THM-4035-sixty-clock-separation-and-finite-kakeya-spine.md).
The earlier finite triangular computation is
[the saved Heron/boat calculation](triangular_fibonacci_heron_boat_kps.out).

In particular,55=F_10=T_10 does not identify the evaluator of the Paley factor
7q-55. The current [Paley run packet](wheel_run_interface_20261014h.md) obtains
A_8=2q(q-1)(q-3)(q-7)(7q-55) by injective character-sum reduction. No map from
triangular Fibonacci values to those signed graph sums is supplied here.

## 2. An actual quotient, rather than a shared value

Write

\[
M(n)=\frac{3n+e}{2},\qquad e=n\bmod2,\qquad
T(n)=\frac{n(n+1)}2,\qquad R(n)=-1-n.                 \tag{1}
\]

For ordinary integers, M is ceil(3n/2). These formulas also define continuous
maps on the2-adic integers; for a rational with odd denominator, its parity
means its residue modulo2. We have

\[
M R=R M,\qquad T R=T.                               \tag{2}
\]

Indeed parity reverses under R. Moreover T(n)=T(m) in the2-adic field iff
(n-m)(n+m+1)=0, hence iff m=n or m=R(n). Therefore there is a well-defined
quotient map

\[
\Phi(T(n))=T(M(n)).                                 \tag{3}
\]

To show that every2-adic target has exactly this two-point fibre, use the
even representative n=2z and the chart

\[
Q(z)=T(2z)=z(2z+1),\qquad
Q(z)-Q(w)=(z-w)(2z+2w+1).                           \tag{4}
\]

Thus Q is an isometry. On every finite ring modulo2^k it is injective,
hence bijective; compatible inverses prove that Q is a bijection of Z_2.
There is exactly one even representative and one odd reflected representative
above every T-coordinate. In this chart the quotient is

\[
Q^{-1}\Phi Q(z)=C(z)=
\begin{cases}3z/2,&z\text{ even},\\-(3z+1)/2,&z\text{ odd}.\end{cases}    \tag{5}
\]

This is a specified signed map; it is not the ordinary positive Collatz map.
The two branches each map their parity coset bijectively onto Z_2, with
inverses2w/3 and-(2w+1)/3 respectively. Consequently every finite binary
branch word occurs in the formal2-adic dynamics. Additional safety tests
can still exclude words; formal realizability alone supplies no real-safety
or ordinary positive-source certificate.

For example, the M-cycle(-2/5,-3/5) becomes the single fixed point-3/25
under T. The complete fixed-point set of Phi in Z_2 is{0,-3/25}: solve
M(n)=n or M(n)=-n-1, obtaining n=0,-1,-2/5,-3/5. This is a change of
period under a quotient, not the discovery of a new positive-integer cycle.

## 3. The missing sheet is exactly one bit

Retain a bit sigma with

\[
n=\begin{cases}2z,&\sigma=0,\\-2z-1,&\sigma=1.\end{cases}
\quad
(z,\sigma)\longmapsto(C(z),\sigma\mathbin{\mathrm{xor}}(z\bmod2)).       \tag{6}
\]

Equation(6) is conjugate to the full marked M dynamics. If c_t=n_t mod2
is its original carry word and a_t=T(n_t), then

\[
d_t:=a_t\bmod2=c_t\mathbin{\mathrm{xor}}c_{t+1}.       \tag{7}
\]

Thus the triangular parity word is the binary derivative of the Mahler carry
word. It identifies exactly a word and its complement. Given d and c_0,
successive XOR recovers every c_t. The restoration bit is necessary even
if the entire future triangular orbit is retained.

If d has least period m, summing(7) over one period gives
c_(t+m)=c_t xor(sum d_j mod2). Therefore the least period of c is exactly
m when the period has even weight, and2m when it has odd weight. A smaller
period of c would induce a smaller period of d. This is an elementary
two-sheet doubling law with a specified cocycle. It is not Bott periodicity.
The [marked Bott packet](collatz_bott_marked_periodicity_20261007d.md) retains
the distinct Clifford/stable-homotopy categories and their own hypotheses.

The sheet is indispensable to Mahler safety. For the inherited real-tail
functional

\[
Y_t(c)=\sum_{j\ge0}c_{t+j}(2/3)^{j+1},\qquad
Y_t(1-c)=2-Y_t(c),                                  \tag{8}
\]

strict safety Y_t<1 does not descend through T. The exact rational hostile is

| Source | Carries | Three tails | Triangular orbit |
|---|---|---|---|
| -9/19 | (100) repeated |18/19,8/19,12/19| -45/361,-30/361,-39/361 |
| -10/19 | (011) repeated |20/19,30/19,26/19| the same orbit |

One is strictly safe at every suffix and the other fails at every suffix.
Neither is an ordinary positive integer. By THM-2228, the separate ordinary
positive stabilization premise is still required for a Z-number. A future
triangular-orbit statistic cannot repair the lost sheet, because these two
entire orbits agree.

## 4. Exact precision and the directed graph interface

For distinct a,b in Z_2 with a=b mod2, (4)--(5) give

\[
v_2(\Phi(a)-\Phi(b))=v_2(a-b)-1.                    \tag{9}
\]

The correct finite interface is therefore

\[
\Phi_k:\mathbb Z/2^{k+1}\longrightarrow\mathbb Z/2^k.  \tag{10}
\]

Each output has exactly two inputs, one in each source parity class. This
also proves that a length-m triangular parity word determines exactly one
source residue modulo2^m. The script computes Q's inverse by lifting one bit
at a time, and exports `phi_residue` and `parity_cylinder` with canonical
residue and precision checks.

The loss of one bit is sharp at every k>=1. Compare n=0 with n=2^(k+1):

\[
T(n)=2^k(2^{k+1}+1)=0\pmod{2^k},\qquad
T(M(n))=3\,2^{k-1}(3\,2^k+1)=2^{k-1}\pmod{2^k}.     \tag{11}
\]

So Phi is not a well-defined fixed-level map on residues modulo2^k.
Wrapping a truncation into a finite functional graph invents transitions.
The legitimate graph connects adjacent precision levels; reversing arrows
gives its two inverse branches. This is distinct from reversing an injective
ordinary integer function graph, and from taking a transpose after forgetting
edge multiplicities.

## 5. Three linear coordinates and the unbounded cost of a tiny carry

For the authenticated branch e=n mod2, direct expansion gives

\[
T(M(n))=\frac{18T(n)+(6e-3)n+3e}{8}.                 \tag{12}
\]

Thus the column(1,n,T(n)) has the exact branch matrix

\[
\begin{pmatrix}1&0&0\\e/2&3/2&0\\3e/8&(6e-3)/8&9/4\end{pmatrix}.       \tag{13}
\]

It is the degree-two affine lift, with the divisibility/parity guard retained.
The span{1,n,T} is the smallest rational vector space of polynomial
observables containing1,T and closed under both branch pullbacks: the e=0
formula already forces n. This is a statement about linear polynomial-state
representations, not a lower bound on arbitrary encodings.

There is no rational polynomial P with P(T(n))=T(M(n)) for all nonnegative
integers n. Growth forces P to have degree one and leading coefficient9/4;
on even n the residual in(12) is-3n/8, not a constant. The exact scalar
quotient(3) therefore exists without a polynomial scalar recursion.

For doubling D(n)=2n, the operation order has the exact discrepancy

\[
M D(n)=3n,\qquad D M(n)=3n+e,
\]
\[
T(DM(n))-T(MD(n))=e(3n+1).                          \tag{14}
\]

One forgotten parity bit causes an unbounded error in the quadratic
observable. Conversely, retaining the branch and the three coordinates
permits exact concatenation. This is the useful compression contract;
neither a shared numeral3 nor an octal identity supplies it.

## 6. Reproduction and retained scope

The script checks signed rational numerators-128..128 over denominators
1,3,5,7,19; all dyadic residue graphs through output precision11; ordinary
sources0..256 through12 triangular parity steps; all primitive binary words
of lengths1..8; and the complete small arithmetic controls stated above.
All arithmetic is exact. Boolean/float aliases, even denominators, invalid
bits and noncanonical source residues are rejected.

```text
python -B 04-computation/experiments/recursion_carry_bridges_20261014i.py
python -B -O 04-computation/experiments/recursion_carry_bridges_20261014i.py
```

No computation is an all-index triangular-Fibonacci classification, an
all-depth safety proof, a change to the Paley A8 character-sum theorem, or
a Collatz stopping certificate. The next genuine obligation is to control
the **marked** orbit's real tails together with ordinary positive source
stabilization. Deleting the sheet already makes the strict-tail component
undetermined by the quotient alone, even with its full infinite future
supplied. The displayed rational hostile does not distinguish two actual
Z-number candidates: both fail the separate positive-integer condition.

Normal and optimized runs reproduce **17,180 checks**, with identical saved
JSON. LF-normalized SHA256:
`0757f91168ab6c0bb320f7c1a06af0761eede9541fe180317ad342a727c9b76e`.
Independent proof/API and optimized-versus-saved audit: **PASS**.
