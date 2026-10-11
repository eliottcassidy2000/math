# Atomic source weights pay unbounded ranks; credit becomes an integer rank

2026-10-11. **PROVED:** the countable common-rank implication, its sharp
atomic integrability threshold, and an explicit integer rank for every
legal H/G/run credit episode. **FINITE-EXACT:** the declared controls.
**OPEN:** a ROOT-valid rank on all positive sources, a native bridge out
of the unresolved components, and universal Collatz.

## Incoming work and the operation being reused

The checkpoint `cd090a8e54` was fetched and its
[integration](collatz_shapes_and_completion_20261011.md),
[run receipts](collatz_projective_receipts_20261009.md),
[credit interface](collatz_run_credit_interface_20261011.md),
[policy rank](collatz_receipt_flow_geometry_20261009.md),
[marked geometry](simplex_sphere_cube_geometry_20261009.md), and
[arithmetic boundary](collatz_boundary_fibres_20261011.md) were read before
this continuation. Their actual guards, common-future clocks, and terminal
limitations are inherited, not independently rediscovered here.

The nearest proved mechanism is the finite common-rank refinement theorem.
The hostile is the fixed boundary point (e_x,-1,-1), which rules out a
uniformly decreasing continuous forward rank on that compact closure.
The corrected near miss is reading a terminating paid episode as a ROOT
proof. The least-used sidecar is a source measure that can integrate an
unbounded rank. Board: **source atom / common rank / native guard / credit /
boundary growth / terminal**.

The map takes the incoming exact potential into an integer rank; it
preserves legal action progress while discarding physical step cost. The
actual receipt and clock remain sidecars. Its cheap hostile is the
neutral C/G loop, whose rank cannot decrease around the whole cycle.
The second map weights source-indexed ranks by the explicit rational
prior from [R6](collatz_renewal_source_20261011.md). It preserves every
source atom but supplies neither an available action nor its ROOT endpoint.

## A1. A common rank need not be bounded across all sources

Let a countable authenticated proof-state system have genuine terminal
states, and let h be a nonnegative, finite-valued function, zero at those
terminals. Every nonterminal state must have an available authenticated
probability row of total mass one; the process chooses such a row after
each history until it terminates. Suppose every permitted conditional row at every
nonterminal state z satisfies

\[
                \mathbb E[h(Z_{t+1})\mid\mathcal F_t]
                    \le h(z)-1.                         \tag{1}
\]

Transitions may depend on the entire past. The inequality includes finite
conditional expectation; it is not merely a statement about support.
Then for each initial state z,

\[
 \mathbb E_z\tau\le h(z),\qquad
 \Pr_z(\tau>t)\le \min(1,h(z)/(t+1)).                   \tag{2}
\]

**Proof.** Kill h at the terminal. Iterating conditional expectations in
(1) gives E h(Z_T)+sum_{t<T}Pr(tau>t)<=h(z). Nonnegativity and monotone
convergence prove the expectation bound; Markov's inequality applied to
tau>=t+1 proves the second claim. No global upper bound on h is used.
Moreover a supported successor has rank at most h(z)-1. Selecting such
successors yields a finite deterministic terminal path; with rational
ranks and an effective enumeration of positive-probability authenticated
successors this is an effective extraction. Its physical Collatz word still needs the receipt
exporter and its clock bounds. QED.

The same proof works under finite or countable hidden-state refinement
when each projected conditional row retains (1) and the terminal fibre is
genuine. This is the elementary unbounded extension of the incoming
finite-rank proof. It loses its **uniform geometric** estimate: countdown
states j->j-1 with h(j)=j have survival1 at every fixed horizon if the
start is chosen large enough. The zero-terminal condition cannot be
assigned to an unresolved label.

This evades the *form* of the compact continuous-rank obstruction by
allowing growth at the boundary. It does not construct the missing
Collatz rank or relax native transition authentication.

## A2. Which source weights can pay that growth?

Write n=2m+1 and use

\[
 \nu_r(m)=(m+1)^{-r}-(m+2)^{-r},\qquad r=1,2,\ldots.
\]

The elementary comparison

\[
 \frac{r}{(m+2)^{r+1}}\le\nu_r(m)\le
 \frac{r}{(m+1)^{r+1}}
\]

shows the sharp threshold

\[
       \sum_{m\ge0}\nu_r(m)(m+1)^\beta<\infty
                         \quad\Longleftrightarrow\quad\beta<r.       \tag{3}
\]

Thus if an independently proved common rank has initial bound
h(z_{2m+1})<=C(m+1)^beta with beta<r, its average terminal time under
nu_r is finite by(2). This is a **conditional payment theorem**, not a
proof that such a ROOT rank exists. The bound concerns the initial rank;
it does not assume an unproved bound on the entire physical orbit height.

There is a second exact zeta connection. For N=m+1, the law has
Pr(N>=k)=k^{-r}. Hence for integer1<=s<r,

\[
 \mathbb E_{\nu_r}N^s
    =\sum_{k\ge1}\frac{k^s-(k-1)^s}{k^r}.                \tag{4}
\]

In particular E_{nu_r}N=zeta(r) for r>1 and
E_{nu_r}N^2=2zeta(r-1)-zeta(r) for r>2. The finite identity includes the
boundary term-N^s/(N+1)^r; omitting that term at a finite cutoff is wrong.
Under nu_2 the linear initial-rank budget is zeta(2), while the refinement
alias constant is2zeta(3). Under nu_3 the linear budget is zeta(3).
All individual weights remain rational. These are series identities, not
an irrationality argument or a theorem about every odd zeta value.

## A3. An actual all-height integer rank for the H/G/run library

Inherit the authenticated actions and guards from the credit packet:

\[
 \mathcal E(x,k)=(x+5)(9/8)^k,\quad
 \kappa=2349/2560,\quad\kappa_R=57/64<\kappa.
\]

At completed nonempty legal action prefixes x>=111,k>=0, so E>=116.
H multiplies E by at most kappa and earns two credits; R multiplies it
by at most kappa_R and earns one; G preserves E and spends one.
The episode starts from immutable source N with zero credit; guards and
the original-source payment inequality are unchanged.

For rational E>=116 define the computable nonnegative integer

\[
 b(E)=\max\{j\ge0:116\le E\kappa^j\},\qquad
                  \boxed{\mathcal R(x,k)=3b(\mathcal E(x,k))+k.}      \tag{5}
\]

If E'<=kappa E and E'>=116, then b(E')<=b(E)-1. The proof is immediate:
if b(E')=j then116<=E'kappa^j<=E kappa^{j+1}. Therefore

| Action | Credit change | Guaranteed rank decrease |
|---|---:|---:|
| H | +2 | at least1 |
| R | +1 | at least2 |
| G | -1 | exactly1 |
| Incoming bundled GR action, zero net credit | 0 | at least3 |

The last action is the explicitly authenticated ray
x=251+288t -> m=223+256t, t>=0, from the incoming integration. Its
energy ratio is at most57/64. It is a declared atomic action, not a
retroactive permission to spend an unavailable credit at its intermediate
state. It also adds no new source-guard coverage, as the incoming note
explains.

Thus every legal episode has at most R at entry many completed actions.
This repackages the existing finite-episode bound as one common integer
rank on an unbounded state space. Any refinement that only narrows native
guards inherits it. A random or history-dependent selector inherits the
same rank because every allowed edge decreases it. Internal states of a
receipt are not being counted as action endpoints.

The actual counter-control is C(n)=(8n-5)/9, whose proposed extra credit
keeps E constant. On the authenticated loop283->251->283, C increases R
by1 and G decreases it by1. Thus adding C invalidates the common rank,
exactly exposing neutral recycling instead of hiding it in nonincreasing
energy. The native H/G and R components also remain separated as proved
by the incoming modulo16/modulo9 argument; a shared rank creates no bridge.

**Terminal scope.** An episode may stop because no action is available.
This is not necessarily ROOT. Formula(5) is a uniform termination
certificate for the declared episode grammar, not for all Collatz orbits.
To use it for universal proof coverage, every stopped obligation must
receive a supplied terminal or a separately ranked authenticated bridge.

## A4. The source prior pays the actual episode rank

There is an unconditional averaging statement for this particular
controller, without assuming a universal ROOT rank. Sample the source
N=2m+1 from nu_r conditioned on m>=55, and start with zero credit.
The conditioning mass is exactly 56^{-r}. Sources with no permitted
first action have episode length zero. For j>=1 set

\[
 M_j=\left\lceil58\kappa^{-j}-3\right\rceil,
       B_N=b(N+5),\qquad\kappa=2349/2560.
\]

The definition of b gives the exact equivalence B_N>=j iff m>=M_j.
Telescoping the source prior, and using
M_j+1>=58kappa^{-j}-2>=56kappa^{-j}, gives

\[
 \Pr(B_N\ge j\mid m\ge55)
   =\left(\frac{56}{M_j+1}\right)^r\le\kappa^{rj}.       \tag{6}
\]

For any legal history-dependent selector, let L be the number of
completed episode actions. A3 gives L<=3B_N pointwise. Therefore

\[
 \Pr(L>t)\le\kappa^{r(\lfloor t/3\rfloor+1)},\qquad
 \mathbb E L\le\frac{3\kappa^r}{1-\kappa^r}.            \tag{7}
\]

Under nu_2 the last bound is exactly
16553403/1035799<16 completed actions. This is a bound under the declared
source distribution, including sources whose episode stops immediately.
Conditioning further on a guard or on survival changes the law and does
not retain this number automatically. The source measure is explicit;
no random-trajectory hypothesis or Haar-to-integer transfer enters (6).

The geometric tail here is consistent with A1's warning: it averages
over a specified source law, whereas no geometric bound uniform over
all fixed starting sources was asserted. The terminal remains the end of
an episode, which can be ungrounded. Thus (7) pays controller work, not
the missing ROOT bridge or the physical Collatz stopping time.

## Next obligation and reproduction

The next concrete target is a native bridge that strictly lowers (5), or
comes with a second well-founded account that pays any temporary increase.
The new source priors can integrate polynomial initial growth of such a
rank and reconstruct source scores with the cubic alias bound. They cannot
turn a neutral loop or an ungrounded component into a terminal.

[Program](../../04-computation/experiments/collatz_atomic_rank_interface_20261011.py),
[exact output](collatz_atomic_rank_interface_20261011.json). The controls
check the potential-to-rank conversion and exact moment boundary terms;
the underlying full receipt authentication is inherited from the linked
packets, not claimed to have been independently reimplemented here.
Normal and optimized runs agree byte for byte on 2,633 always-active
exact checks, including 768 guarded action examples, 320 finite
moment identities, and 96 source-rank tail thresholds. No independent
peer audit is claimed.

```text
python3 -B 04-computation/experiments/collatz_atomic_rank_interface_20261011.py
python3 -B -O 04-computation/experiments/collatz_atomic_rank_interface_20261011.py
```
