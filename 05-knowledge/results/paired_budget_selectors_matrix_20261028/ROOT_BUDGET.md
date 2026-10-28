# The complete one-step paired ROOT budget

**Status: PROVED, by an explicit finite induction certificate.**
The inherited half-budget question is resolved, and the same mechanism gives

\[
\boxed{D_{\rm pair}\le\min\{S/2+3,\ S/3+4\}.}
\]

Both intercepts are sharp at target `u=1`, cutoff `N=3`, where `D=S=6`.
The coefficient `1/3` is optimal among bounds `D <= c S + O(1)` for these
complete one-step fibres. This does **not** prove a uniform multilevel ROOT
bound or universal Collatz entry. Nor does it permit arbitrary replacement
of the original source weights by new weights.

The same grammar also proves the sharp joint amplitude bound
`Q(sum_(P(s)=u) (-1)^t A_N(s)) <= 13`. Triple extensions leave the total
unchanged; the 702 base states have maximum 13. The inherited target-seven
family attains it. Thus the older merged-charge hostile of norm 13 is now
matched by an all-cutoff upper bound.

## Inheritance and exact scope

Use the primitive positive odd sources `s not congruent to 5 mod 8`, native
map `P`, integer clock `t`, and Eisenstein charge `A_N(s)` of
`../root_ports_closed_carriers_20261028/ROOT_RAYS.md`. In particular, with
`B(s)=4s+1`, `omega^2+omega+1=0`, and `ell=omega-omega^2`,

\[
A_N(s)=\sum_{B^j(s)<2N}\omega^{B^j(s)}\ell+
\sum_{2N\le B^j(s)<4N}\omega^{B^j(s)+1},\qquad
z_s=(-1)^{t_s}A_N(s),\quad Q(a+b\omega)=a^2-ab+b^2.
\]

Each source has `Q(z_s) in {0,1,3}`. Keep its actual address, multiplicity,
terminal phase and clock. Source `1` keeps its separate ROOT-stop tag and
clock zero; its even sign agrees with the native self-edge's sign.

For each complete inverse fibre `P^{-1}(u)` and positive integer cutoff `N`,
let `S=sum Q(z_s)`. For a ray word `(z_k)`, put

\[
T_k=\sum_{j\ge k}z_j,\qquad
\delta_k=Q(T_k)-Q(z_k)-Q(T_{k+1}).
\]

The initial odd/even split costs the positive part of
`Q(T_odd+T_even)-Q(T_odd)-Q(T_even)`. Pair the adjacent odd/even internal
splits at the same inverse-ray index and charge `max(0,delta_odd+delta_even)`.
When `u=2 mod 3`, first remove the unmatched odd head and charge its positive
increment separately. Their sum is `D_pair`. This is the pairing and the
actual dyadic comb of the inherited theorem, not a new order selected to
reduce the debt.

The closest proved mechanism was the single-ray `D <= 2S/3` bound. Its
sharp families rule out an isolated-ray half bound. The target-seven hostile
family also rules out a constant paired debt. The new ingredient is the
**joint order of the two rays' cutoff events**; independent ray estimates
lose this coordinate.

## Normalize the two rays without moving their clocks

Write the even-clock ray as

\[
e_k=(C64^k-7)/9,\quad C=2^{\tau_e}(3u+1),\quad
C\equiv16\pmod{72},\ C\ge16.
\]

Its matching odd-clock sources are exactly

\[
o_k=(2C64^k-5)/9=2e_k+1.
\]

For `u=0,1 mod 3`, both start at `k=0`. For `u=2 mod 3`, the odd ray also
has `k=-1`, its clock-one head; the paired part begins at odd clock seven
and even clock six. In that case `C>=2176`, and the head is integral. All
subsequent clocks increase by six. Set

\[
\alpha=e_0\pmod3,\quad\beta=2\alpha+1\pmod3.
\]

Then `e_k=alpha+k mod 3` and `o_k=beta-k mod 3`. Do not average over
`alpha`: the proof tests all three phases separately.

For shell `i=3k+j`, the actual `B^j` addresses are

\[
E_{i,k}=\frac{C4^i-4\,4^{i-3k}-3}{9},\qquad
O_{i,k}=\frac{2C4^i-2\,4^{i-3k}-3}{9}.
\]

## The two event clusters

As `N` grows, all charge changes occur at the following two alternating
clusters for `i=0,1,...`. A threshold `x/c` is crossed when `x<cN`, so
integer-rounding can merge events but cannot reverse their order.

* **A:** `E_(i,k)/2`, then `O_(i,k)/4`, alternating in increasing `k`.
* **B:** `E_(i+1,k)/4`, then `O_(i,k)/2`, alternating in increasing `k`.

In the head case its odd event comes first in each cluster. In B there is
one extra final even event when `i=2 mod 3`.

These order assertions follow directly by subtraction. For example, writing
`j=i-3k`, the within-pair gaps in A and B are respectively

\[
(3\,4^j+3/2)/18>0,\qquad (3\,4^j-3/4)/9>0.
\]

The gap from the odd event to the next even event in A is
`((15/16)4^j-3/2)/18>0` whenever that next index exists (`j>=3`);
in B it is `((15/16)4^j+3/4)/9>0`. The head-to-index-zero gaps are positive
by the same formulas with `k=-1`. Distinct A/B clusters are separated:
without a head this follows from `C>=16>8`; with a head it follows from
`C>=2176>128`. The largest event in A is below the first in B, and the
largest in B is below the first in the next A.

Between clusters nothing changes. After B's extra even event the word is
exactly the initial word of the next A cluster. Before the first full paired
shell there are at most two nonzero sources (the even source zero and the
optional odd head). Their individual energies are at most three. Since their
only possible positive debt is a cross term bounded by their total energy,
`D<=S<=6` gives both stated bounds in this initial regime.

## A finite grammar for every cutoff

Let `r=i mod 3`, `U_r=1+omega+...+omega^(r-1)`, so
`U_0=U_3=0`, `U_1=1`, `U_2=1+omega`.
For each common source index, even and odd coefficients before/after its
cluster event are `b_e,f_e,b_o,f_o` below:

| Mode | `f_e` | `b_e` | `f_o` | `b_o` |
|---|---|---|---|---|
| A | `ell U_r + omega^(r+1) - omega^(r+2)` | `ell U_r + omega^(r+1)` | `ell U_r + omega^(r+1)` | `ell U_r` |
| B | `ell U_(r+1) + omega^(r+2)` | `ell U_(r+1)` | `ell U_r + omega^(r+1) - omega^(r+2)` | `ell U_r + omega^(r+1)` |

This is obtained by summing the full earlier `B` shells. Their length differs
by `3k`, which contributes zero to the geometric phase sum. In A an even
event changes an upper contribution to a lower one, while an odd event adds
an upper contribution. In B these roles reverse.

Every event cut consists of `m` completely passed pairs, at most one partly
passed pair (`epsilon=0,1`, even passed first), and `n` unpassed pairs. Thus

\[
z^e_k=\omega^{\alpha+k}
\begin{cases}f_e&k<m+\epsilon\\b_e&k\ge m+\epsilon,\end{cases}
\quad
z^o_k=-\omega^{\beta-k}
\begin{cases}f_o&k<m\\b_o&k\ge m.\end{cases}
\]

Here `0<=k<m+n+epsilon`. If present, the head is
`-omega^(beta+1) f_o`; before its event it is instead
`-omega^(beta+1) b_o`, in which case `m=epsilon=0`. Zero trailing sources
may be retained: they change neither energy nor debt. The last extra B event
was already assigned to the next A word. This proves coverage of the grammar
for every target and every integer cutoff, not merely the sampled ones.

## Finite induction certificate

Adding three leading fully passed pairs leaves each ray's total charge
unchanged because `1+omega+omega^2=0`. It also leaves the old paired
increments, initial cross term and unmatched-head term unchanged. Only the
three new paired increments contribute. Their cost depends only on
`r,alpha,epsilon,m mod 3,n mod 3`. Adding three trailing unpassed pairs has
the same zero-sum property and costs no positive paired debt.

The exact complete local table is:

| Mode, r | Energy added by three front pairs | Possible added front debt | Energy added by three back pairs | Added back debt |
|---|---:|---|---:|---:|
| A,0 | 12 | 0 | 3 | 0 |
| A,1 | 12 | 0,1 | 12 | 0 |
| A,2 | 3 | 0,1 | 12 | 0 |
| B,0 | 12 | 0 | 12 | 0 |
| B,1 | 12 | 0,1,2 | 12 | 0 |
| B,2 | 3 | 0 | 3 | 0 |

Consequently both `2D-S` and `3D-S` are nonincreasing under either triple
extension. It remains to check `m,n in {0,1,2}`. There are **702** base
states, including all three phases, both modes, the partly passed pair,
absent head, passed head, and the allowed unpassed-head states. Direct
integer arithmetic in `Z[omega]` gives

\[
\max(2D-S)=6,\qquad\max(3D-S)=12.
\]

The executable certificate is `certificate()` in `check_paired_budget.py`.
It uses no floating arithmetic or disabled Python assertions. The zero-sum
argument proves the recurrence at arbitrary lengths; the script checks its
finite state table and repeated extensions. This closes both inequalities.

## Sharpness and why an unweighted estimate still fails

For target seven take

\[
N_h=\frac{44\,4^{9h+2}-29}{9},\qquad h\ge0.
\]

These integer cutoffs lie after the last even A event and before the last
odd event with `C=88`, `i=9h+2`. They have grammar state
`(A, r=2, alpha=0, m=3h, n=0, epsilon=1, no head)`.
Every even charge is zero. The odd word is

\[
(1+\omega,-\omega,-1)^h,\ 1+2\omega.
\]

Direct suffix summation gives `D_pair=h` and `S=3h+3`. Hence no bound with
coefficient less than `1/3` and a fixed additive constant is possible.
The intercept four is separately forced by `(u,N)=(1,3)`.

This family also explains why pairing is not magically contractive: a ray
can carry zero charge at the chosen cutoff while its partner retains repeated
positive refinement cost. The joint event grammar bounds that cost by the
actual source energy; it does not erase it.

## Remaining depth obligation

At the next ROOT depth the charge of a child is a transported sum from its
own fibre, rather than the original `A_N` leaf charge. Applying this theorem
again requires proving that those **weighted, phase-retaining child words**
belong to a closed grammar with a summable intercept budget. Neither fact
follows from one-step coverage. The next useful target is a source-preserving
grammar for two successive native fibres, with explicit shared cutoff and
integer-clock data. A global averaged norm or an independent rescaling of
the two rays is not a substitute.

That boundary is now concrete. `check_depth_two.py` enumerates the original
sources and keeps both clocks. At cutoff `N=19`, target one has outer words
`(-omega)` and `(-omega,-omega)`: outer paired debt is six but the child
energy is three. The one-step bound cannot simply be applied to these new
words. At `N=57`, target one has outer debt 27, child energy 12 and original
source energy 16; its actual endpoint norm is 39. Even replacing the child
energy with the original source energy fails to repair the naive iteration.
The receipt retains every original address and clock for both witnesses.
This refutes an extension, not the proved one-step theorem.
