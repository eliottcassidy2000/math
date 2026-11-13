# Minimum-budget escape at the unit-ready hole

Status: **PROVED** explicit strong indegree-one defect-two source and an
infinite root-preserving extension family. The minimum-budget assertion uses
the independently audited all-order unit-ready exclusion in the parent
packet. The final next-target unit-ready exclusion has a separate receipt;
it is not needed for the constructive theorem.

Let N=53084231=`405·2^17+71`. This number is attained; the obstruction concerns
the indegree-one unit-ready interface. Increasing the positive defect from1
to2 suffices at the *same exact target*.

## 1. Actual source and its complete rooted state

Take the four-vertex core

```
T = (6,4,8,3),  J = {1,2}.
```

Each vertex of J starts exactly one Hamilton path in T. Attach a transitive
tail, dominated by T, and a fresh vertex p which beats exactly J inside T.
For the tail word w, `1` means tail vertex→p and `0` means p→tail vertex.
Then append r beaten only by p. The final root is r and its sole predecessor
is the named vertex p.

Set `x=1133276`, and use its bits from least to most significant:

```
w = 001110110101001010001.
```

The source has27 vertices and rooted state `(H,s,q,F)`

\[
\boxed{(53084231,16819383,16819385,53084233)}.
\tag{1}
\]

Here s counts Hamilton paths starting at r, q=H(Total−r), and F counts ordered
spanning covers `(P,Q)` with P nonempty ending at r and Q possibly empty.
Thus `q−s=F−H=2`. The graph is strong because the core is strong, J is a
nonempty proper subset, and the final tail bit is1; the appended vertex
preserves strongness. The literal adjacency rows, p, and r are in
`defect_two_inverse_receipt.json`.

The parent packet excludes defect1 for every indegree-one rooted realization
of N. Therefore2 is the least positive defect in that specific rooted class.
No minimum over other interface types is asserted.

## 2. Universal word formula and proof

Define the shell energy

\[
g(0)=g(1)=1,\quad
 g(2a)=g(a)+a+(a\bmod2),\quad g(2a+1)=g(a).
\tag{2}
\]

For any x≥0, use its least-significant-first word, with x=0 denoting an
empty tail. This same fixed core and entry set gives

\[
\begin{aligned}
h(x)&=5g(x)+43x+3(x\bmod2)+18,\\
q(x)&=5g(x)+11x+3(x\bmod2)+4,\\
s(x)&=q(x)-2,\qquad F(x)=h(x)+2.
\end{aligned}
\tag{3}
\]

This is the inherited two-entry funnel mechanism, specialized to a new
prescribed target. Its exact definition and multi-entry boundary are in
[`two_entry_funnels_20261103/TWO_ENTRY.md`](../../two_entry_funnels_20261103/TWO_ENTRY.md)
and `research.py`; the polarized derivation is in
[`polarized_path_calculus_20261103/FUNNELS.md`](../../polarized_path_calculus_20261103/FUNNELS.md).

For clarity, the finite core response is retained, including the joint cut
term that cannot be dropped:

| core coordinate | value |
|---|---:|
| h=H(T) |5|
| P₂, two ordered bins allowing empties |32|
| e, full paths ending in J |3|
| F_J, first bin ends in J, second arbitrary |16|
| G_J, first bin starts in J, second arbitrary |11|
| b, full paths starting in J |2|
| C_J, slot1 ends in J and slot2 starts in J, both nonempty |4|

For the word associated with x, the inherited tail coefficients are
`L₂=x+1`, `B₁=g(x)`, `D₁=1−(x mod2)`. They follow directly by induction on
(2), or from the finite word automaton. The parent's Hamilton count and
ending-root two-bin response are

\[
q=L_2G_J+B_1h-D_1e-C_J,
\qquad E_p(2)=L_2P_2-F_J.
\]

The new root r then has

\[
H=q+E_p(2)-b,
\]

which gives (3). The predecessor starts exactly b=2 full paths: it must
enter the first core SCC through J and then traverse the tail. The
indegree-one birth identity gives `q−s=b` and `F−H=b`. At the specified
address, `g(1133276)=870669`, yielding (1). This proof applies to arbitrary
word length; the27-vertex source does not require a general subset census.

## 3. An infinite family retaining both physical ports

Prefixing one zero to the word inserts one vertex before the old transitive
tail. Its arcs are

\[
T\to z\to\text{old tail},\qquad p\to z,\qquad r\to z.
\]

Every old vertex and old arc is retained. In particular r and its unique
predecessor p stay fixed, and the new source remains defect-two.
After ℓ such steps the address is `2^ℓx`. Since x is even,

\[
g(2^\ell x)=g(x)+(2^\ell-1)x.
\]

Thus on27+ℓ vertices, for every ℓ≥0,

\[
\begin{aligned}
H_\ell&=54397248\,2^\ell-1313017,\\
q_\ell&=18132416\,2^\ell-1313031,\\
s_\ell&=q_\ell-2,\quad F_\ell=H_\ell+2.
\end{aligned}
\tag{4}
\]

All H_ℓ are odd. Apply the inherited two actual births NR, keeping the first
newborn as the active root. The balanced output is

\[
(H',s',q',F')=(3H_\ell+4,H_\ell,H_\ell,3H_\ell+4).
\]

It is primitive for every ℓ, since `gcd(3H_ℓ+4,H_ℓ)=gcd(4,H_ℓ)=1`.
Therefore the family reaches the established primitive balanced controller
and retained-chain interfaces in
[`h17_exact_repairs_20261112/constructors/CHAIN_REUSE.md`](../../h17_exact_repairs_20261112/constructors/CHAIN_REUSE.md).
Here H_ℓ≡2 modulo3 and3H_ℓ+4≡1, so the usual C₃ entry uses one balanced R
before C₃; the general connector supplies the remaining fixed finite bridge. This proves reusable access to those
interfaces, not a return to defect-two or unit-ready at every next target.

## 4. Exact target boundary and controls

At ℓ=0 the next NR count is159252697. The complete unit-ready inverse produced
an exclusion receipt `ready_after_defect_two.json` with16075330 retained
cores. An independent direct-outer-word audit in the sibling ports packet checks
the same16075330 cores and233412223 words, with no hit. This extra exclusion does not invalidate the balanced primitive
output: unit readiness and primitivity are different resources.

The neighboring number159252695=`3N+2` has an independently supplied
unit-ready source in the parent packet. Equation(3) does not claim to produce
that number from the defect-two source. Subtracting the two-unit debt would
require a further actual target-preserving operation.

Reproduction:

```
python3 05-knowledge/results/target_directed_cores_20261112/constructors/defect_two_inverse.py
```

The receipt contains the literal27-vertex graph, the exact core moments,
the target identity, the full root state, and independent subset/endpoint
DP controls for all16 words with integer addresses0 through15. These small
controls test the joint term and the empty-tail boundary; the symbolic
funnel law proves the unbounded formula. The search stops at the first
positive source and makes no completeness claim for its core menu.
