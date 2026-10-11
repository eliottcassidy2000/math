# Atomic actual-step bridges into the run family

**Status: PROVED** elementary template classification and guarded common-future
families; **FINITE-EXACT** replay controls; **OPEN** arbitrary terminal-child
grounding. This packet changes the action granularity and retains every carry.
It does not infer universal Collatz convergence or new coverage over the full
ordinary-descent algorithm.

## 1. Inheritance and the missing interface

Write $U_a(x)=(3x+1)/2^a$, with the letter $a$ legal only when it is the exact
valuation. All exported words are positive-odd, strict first-hit words.

The inherited [five-letter carry interface](collatz_carry_interfaces_20261004.md)
and [translation decoder](translation_phase_decoder_20261005.md) retain these
maps and their **native** source domains:

| Prefix | Map $x=F(z)$ | Source $z=c+Mu$ | Output $x=d+Lu$ |
|---|---|---|---|
| H | $(729z+669)/1024$ | $155+2048u$ | $111+1458u$ |
| G | $(9z+5)/8$ | $11+16u$ | $13+18u$ |
| A | $(81z+85)/128$ | $187+256u$ | $119+162u$ |
| B | $(81z+73)/128$ | $7+256u$ | $5+162u$ |
| L | $(9z-3)/16$ | $219+256u$ | $123+144u$ |
| LG | $(81z+53)/128$ | $219+256u$ | $139+162u$ |

Here and below every progression parameter is a nonnegative integer. `LG` is
the inherited composition L then G, not a new freely adjoined letter. Every
native L output is $11\pmod {16}$, so that G is automatically native there.

The [projective run receipt](collatz_projective_receipts_20261009.md) supplies

\[
R_r(y)=\frac{2^{r+1}(2y+1)}{3^r}-1,
\qquad y\equiv3\pmod4,\quad 3^r\mid2y+1,\quad r\ge4.
\]

It has $1<R_r(y)<y$ and the actual common future

\[
y\xrightarrow{(1)}J=\frac{3y+1}{2}
\xleftarrow{(1^r,3)}R_r(y).
\tag{1}
\]

Its union of native source domains is $y=283+324t$. The previous
[credit interface](collatz_run_credit_interface_20261011.md) proves a real
barrier: direct old-library output cannot feed R after an old prefix of
length at least two, and every R output is $15\pmod {16}$, blocking the five
old native actions. The neutral inverse-G/G cycle is the canonical hostile to
declaring extra credit merely from reversible receipts.

The new operation is instead the **atomic** batch

\[
z\xrightarrow{F}x\xrightarrow{U_a}y
\quad\text{followed by the paid dependency }y\leadsto R_r(y).
\tag{2}
\]

An actual step can change the output's ternary class. The preserved object is
the original source $z$ and a literal common-future receipt, not a controller
label or a residue alone. No neutral recycling is used.

## 2. A complete progression compiler for the declared template

Fix one row $F$ of the table, a positive integer $a$, and $r\ge4$. Put

\[
\tau_r=\begin{cases}7&r\text{ even},\\5&r\text{ odd},\end{cases}
\qquad y_r=(3^r\tau_r-1)/2,\qquad M_r=4\cdot3^r.
\]

These are exactly the positive native R$_r$ sources $y=y_r+M_rt$. For
$z=c+Mu$ and $x=d+Lu$, the batch is native exactly when

\[
3L u\equiv 2^a y_r-(3d+1)\pmod{2^aM_r}.
\tag{3}
\]

Let $g=\gcd(3L,2^aM_r)$. There is no source if $g$ does not divide the right
side. Otherwise, reducing (3) by $g$ gives one class $u=u_0+Tu'$, with
$0\le u_0<T$ and $T=2^aM_r/g$. Thus the complete positive source ray is

\[
z=(c+Mu_0)+MTu',\qquad u'\ge0.
\tag{4}
\]

There is no implicit source selection here: `receipt(packet,z)` checks the
supplied source against (4). The constructor merely describes its domain.
Because $0<c<M$, a negative $u$ would give a nonpositive source. Conversely,
every source in (4) has positive $x$ and positive odd $y$, and the oddness of
$y$ makes $a$ the **exact** valuation. Positivity forces the representative
$y=y_r+M_rt$ to have $t\ge0$ because $0<y_r<M_r$.

Every such batch pays its original source. A native run source is $1\pmod3$;
therefore $2^{-a}=1\pmod3$, so $a$ is even and at least two. H, A, B and LG
already decrease their native positive source, and $U_a(x)<x$ for $a\ge2$.
For the expanding G case,

\[
y=\frac{27z+23}{2^{a+3}}<z
\]

because $a\ge2$ and $z\ge11$. Together with (1), this gives
$1<R_r(y)<y<z$. No extra finite-height cut is needed.

This theorem is complete **only for (2)** with a declared prefix, one actual
valuation, and one specified R run. It is not a classification of all common
futures or all possible connectors.

### Exact all-height R4 valuation classes

Solving (3) at $r=4$ gives:

| Prefix | All positive eligible $a$ | Least source at its least $a$ |
|---|---|---:|
| H | $52\pmod {54}$ | 8795309731532069019 |
| G | $0\pmod {18}$ | 21981259 |
| A | $8\pmod {54}$ | 81851 |
| B | $32\pmod {54}$ | 2106266101767 |
| L | none | — |
| LG | $36\pmod {54}$ | 21972133598683 |

For H, A, B and LG the only ternary condition is respectively
$2^a=61,13,49,55\pmod {81}$; the binary gcd condition is automatic. The
element 2 has order 54 modulo 81, giving exactly the displayed classes. For G
the condition is $2^a=1\pmod {27}$, and the order is 18. These primitive
orders can be proved by LTE: $v_3(2^{2j}-1)=1+v_3(j)$, while odd powers are
$-1\pmod3$. This is an all-height argument, not extrapolation from the
finite controls.

For L, its output is $11\pmod {16}$, whose next actual valuation is 1. It
cannot reach a run source, which requires even $a$. Passing first through the
native G action changes this obstruction: LG has the last row above.

**Carry hostile.** Keeping only the old output class modulo 9 predicts H's
candidate $a=16$ and B's $a=14$; both fail (3). The full carry raises the first
solutions to 52 and 32. The strongest survivor of the coarse argument is a
necessary congruence, not a native bridge.

## 3. An infinite G family, including infinitely many grounded bases

For **every** $j\ge1$ and $t\ge0$, set $a=18j$ and

\[
N_j(t)=\frac{2264\,2^{18j}-23}{27}+96\,2^{18j}t.
\tag{5}
\]

The numerator is integral since $2^{18j}=1\pmod {27}$. Then

\[
N_j(t)\xrightarrow{G}
\frac{2^{18j}(283+324t)-1}{3}
\xrightarrow{U_{18j}}283+324t
\leadsto223+256t.
\]

The complete actual common-future words are

\[
N_j(t)\xrightarrow{(1,2,18j,1)}425+486t
\xleftarrow{(1,1,1,1,3)}223+256t.
\tag{6}
\]

Thus the source can first increase under G and nevertheless pay the unchanged
original source atomically. A primitive zero-credit G action would be rejected
by the inherited credit controller; that does not invalidate (6), whose
contract is at the batch endpoints.

For $j=1$, (5) is $21981259+25165824t$. At $t=0$, use the already frozen
strict ROOT word of 223 from the run packet. It begins with the right word in
(6). Replacing those five letters gives, for every $j\ge1$,

\[
\operatorname{ROOT}(N_j(0)):
\quad (1,2,18j,1)\;\Vert\;\operatorname{ROOT}(223)[5:].
\tag{7}
\]

This is an infinite fully grounded family with **23 odd edges** and total
valuation sum (Terras shortcut cost) **$18j+43$**. The first source has cost 61.
For the unaccelerated convention $C(n)=3n+1$ at an odd state and $C(n)=n/2$
at an even state, the step count instead adds the 23 odd operations, giving
**$18j+66$**. No orbit search is used
for any $N_j(0)$. In contrast, (6) for arbitrary $t$ exports an obligation at
$223+256t$ and does not assert that obligation is grounded.

### Exact comparison with declared old entry banks

Every source in (5) is $75\pmod {256}$, so it misses all native H, A, B and L
entry guards. Moreover

\[
N_j(t)\equiv2-j\pmod3.
\tag{8}
\]

For example, $2^{18j}=1+27j\pmod {81}$ proves (8). The parameter period has
exactly one factor of 3, so in a fixed $1\pmod3$ source class the parameters
are uniform modulo the remaining ternary precision. The old inverse bank
G1/G5/G17 has pairwise disjoint native classes $2\pmod3$, $4\pmod9$, and
one class $1\pmod9$ modulo $2187$; see
[the exact child normal form](collatz_child_normal_forms_20261007.md).
Consequently the proportions of the parameter ray **outside that bank** are

\[
\begin{array}{c|ccc}
j\pmod3&1&2&0\\\hline
\text{complement}&485/729&1&0.
\end{array}
\tag{9}
\]

The R domain is contained in the G5 native class. Hence the complements in
(9) also miss R, while already missing every old funded first action. The
first source $21981259$ is an explicit such entry. These are exact conditional
proportions in each declared progression, not a density statement about all
unresolved integers. The actual prefix $(1,2,18j)$ already descends, so there
is **no additional source coverage over an algorithm accepting all ordinary
first descents**.

## 4. A and L cascades with a reused grounded terminal

The A family admits one extra run letter without tightening the source guard:

\[
81851+131072t\xrightarrow{A}51797+82944t
\xrightarrow{U_8}607+972t\leadsto319+512t.
\tag{10}
\]

Indeed $2(607+972t)+1=243(5+8t)$, so R5 is native. The actual words to
$911+1458t$ are $(1,2,1,3,8,1)$ and $(1^5,3)$. The known ROOT223 path reaches
319 after its first six letters; its authenticated suffix therefore supplies
ROOT319. At $t=0$, (10) yields a strict ROOT81851 word of rank 18 and cost 45.
The old inverse-bank complement on this source ray is $1214/2187$; this is a
comparison with G1/G5/G17, **not** new mass over the old A entry itself.

The one-step L impossibility also has a precise surviving cascade:

\[
21972133598683+35184372088832t
\xrightarrow{LG}\cdots\xrightarrow{U_{36}}607+972t
\leadsto319+512t.
\tag{11}
\]

Its source word is $(1,2,1,1,38,1)$ and child word $(1^5,3)$, again meeting
at $911+1458t$. At $t=0$, the same supplied ROOT319 gives rank 18 and cost 73.
This crosses the L one-connector barrier using an explicitly retained extra
native action; it is not an assertion that the earlier barrier was false.

## 5. A native bridge for the incoming integer rank

The incoming checkpoint `9449eb146`,
[`collatz_atomic_rank_interface_20261011.md`](collatz_atomic_rank_interface_20261011.md),
converts the inherited H/G/R potential into an integer rank. Its relevant
definition is

\[
\kappa=2349/2560,\quad
b(E)=\max\{d\ge0:116\le E\kappa^d\},\quad
\mathcal R(x,k)=3b((x+5)(9/8)^k)+k.
\]

The new G batch preserves the current nonnegative credit $k$ and gives an
actual missing bridge for this rank. With $y=283+324t$ and $a=18j$,

\[
\frac{R_4(y)+5}{N_j(t)+5}
=\frac{16y+89}{6(2^a y+14)}
\le\alpha_j:=\frac{228}{N_j(0)+5}<\kappa.
\tag{12}
\]

The derivative has the sign of $224-89\cdot2^a<0$, so the bound is sharp at
$t=0$. It decreases further with $j$. Let

\[
d_j=\max\{d\ge0:\alpha_j\le\kappa^d\}.
\]

If $E'\le\kappa^{d_j}E$ and $E'\ge116$, then the definition of $b$ gives
$b(E')\le b(E)-d_j$. Every child here is at least 223, so the lower energy
domain is retained. Hence the atomic batch reduces $\mathcal R$ by at least
$3d_j$ at **unchanged credit**. Exact rational comparison gives
$d_1=133$, and thus every member of every G family pays at least **399 rank
units**. The uniform bound is sharp: at $j=1,t=0,k=1$, the exact ranks are
$427\to28$, a drop of 399. At the same source with $k=0$, they are
$423\to21$, a drop of 402; the bound is not asserted to be an equality
on every input.

The primitive intermediate G can be unfunded; only the authenticated atomic
edge is adjoined to the ranked action alphabet. The proof does not authorize
arbitrary subdivision into credit steps. Adding this edge preserves finite
episode termination by the rank, but supplies neither a next action at every
resulting child nor universal ROOT grounding. The rank measures completed
actions, not the valuation sum $18j+43$ or unaccelerated step count
$18j+66$ of the grounded base words.

## 6. Receipt boundaries, representation map and verification

For G, A, B the actual source word is their inherited literal word followed
by $(a,1)$. For H it is `HEAD` followed by $(a+2,1)$: `HEAD` ends at $4x+1$,
and the next valuation is $a+2$. For LG it is $(1,2,1,1,a+2,1)$ by the same
sibling identity. These identities, the native source guard and odd final
endpoint force every indicated valuation to be exact. Their common endpoint
is greater than 1, so a prior ROOT would have remained at 1 under exact
formal continuation and cannot occur. The exported words are therefore
strict pre-ROOT words. A supplied child ROOT certificate has the indicated
child word as its unique prefix, making the substitution first-hit as well.

| Source representation | Map / preserved predicate | Information lost without sidecar | Decisive control |
|---|---|---|---|
| Native affine prefix $F$ | Linear congruence (3) preserves actual next valuation and R guard | Output mod9 omits higher carry precision | H16 and B14 fail |
| Expanded action sequence | Atomic receipt preserves common future and original-source payment | Stepwise credit bookkeeping may reject an internally rising action | G18 increases first, pays after the batch |
| Paid child | Supplied ROOT prefix substitution | A smaller child is not itself a ROOT certificate | Wrong-child and ROOT-padding words rejected |
| L one-step failure | Retain the full LG prefix | One failed template does not exclude longer cascades | (11) succeeds |

The production APIs are `compile_bridge`, `receipt`, `at`, and `discharge`.
All cached packets are recompiled; exact integer types precede equality checks.
No function silently selects another source or discovers a terminal orbit.

Reproduce from the repository root:

```text
python -B -X utf8 04-computation/experiments/collatz_atomic_run_bridge_20261011b.py
python -B -O -X utf8 04-computation/experiments/collatz_atomic_run_bridge_20261011b.py
```

The finite universe is all six declared prefixes and $1\le a\le216$ for
the R4 classification; all six prefixes, $1\le a\le108$, $4\le r\le8$ for
the general compiler; 18 explicit G exponent families; complete finite
ternary residue comparisons; frozen-child proof splices; and malformed API
controls. This gives 112 classified R4 lifts and 100 general-run lifts.
The independent paths are literal valuation replay, the inherited native
controller evaluator, exact group order cycles, and the linear congruence
compiler. The incoming-rank delta adds exact power inequalities and the
sharp retained-credit equality control. Normal and optimized runs give
**4,238 checks**, with deterministic
output saved in [the output record](collatz_atomic_run_bridge_20261011b.out).

Correction lineage: the pre-freeze draft called the valuation sum an
"ordinary step cost." The two clocks are now explicitly separated above;
the actual words, valuation sums and rank claims were unchanged.
