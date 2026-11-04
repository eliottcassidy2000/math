# All-tail funnel inversion and the actual component grammar

**Status: PROVED algebraic inverse; FINITE-EXACT core menus.** This packet
classifies the stated funnel families, not all marked components or all
tournaments. It does not prove a reusable ternary return or the `{7,21}`
classification. The successful component replacements in `../inverse/` leave
these smaller families.

## Inheritance and the question

The closest proved mechanism is the arbitrary-tail one-port funnel calculus,
`polarized_path_calculus_20261103/FUNNELS.md`, especially its F3 identities and
the complete unique-start grammar implemented in that folder's `recognize.py`.
The current object is the named strong suffix in
`triangle_return_170807_20261103/ZERO_READY.md`: a component $B$ at mark 31,
with $h=2195,F=42153$, preceded by TT2 and the zero-ready outer root 39.
Its total is $h+4F=170807$. The hostile example is the previously exhausted
fixed two-module continuation in `triangle_return_170807_20261103/TRIANGLE.md`.
The corrected near miss is that changing an outer module's scalar contribution
need not transport its marked component response. The less-used sidecar here
is the full polarized path-bin calculus, especially $S_w(2)$, not just $h,F$.

The live concepts are: the zero-ready clock; actual SCC grammar; polarized
path-bin fields; adaptive marks; native edge compilation; and the unit-F dial.
This packet asks whether the current component's recursive cores, with arbitrary
new funnel tails and marks, can supply the two immediate scalar targets. Every
operation below has an explicit actual-graph meaning; scalar hits would still
need named alignment and the joint attachment certificates handled elsewhere.

## A complete inverse on one fixed marked core

Let $U$ be a nonempty tournament with marked vertex $w$. Write

\[
h=H(U),\quad e=\#\{\text{Hamilton paths ending at }w\},\quad
P=P_U(2),\quad F=E_w(2),\quad S=S_w(2).
\]

Here $P_U(2)$ counts ordered covers by two directed paths, allowing empty
bins; $E_w(2)$ requires the first bin to be nonempty and end at $w$;
$S_w(2)$ requires it to start at $w$. A funnel appends a transitive tail
$v_1,\ldots,v_t$, places $U\to\{v_i\}$, and adds a fresh root $r$ with
$r\to w$, $U-w\to r$. Tail bit $ε_i=1$ means $v_i\to r$, and
bit 0 means $r\to v_i$. If $U$ is strong, a nonempty tail gives a strong
output exactly when its last bit is 1. An empty tail is allowed when
$|U|>1$; the singleton/empty-tail case is transitive 2 and is excluded from
strong-component menus.

For a nonempty strong tail, set

\[
z=\sum_{i=1}^t\epsilon_i2^{i-1}\ge1,\quad L=z+1,\quad
D=\mathbf1_{z\text{ even}},\quad B=B_{\epsilon}(1).
\]

The six-state before/after automaton of inherited F3 computes $B$. Equivalently,
assign every tail position to a before or after bin. The before bin must be
empty or end in bit 1; the after bin must be empty or start in bit 0. Thus
$0\le B\le2^t\le2(L-1)$. F3 at arguments 1 and 2 gives

\[
\boxed{h'=LS+Bh-De,\qquad F'=LP-F.}\tag{CI1}
\]

The separate empty-tail boundary is

\[
\boxed{h'=S+h-e,\qquad F'=P-F.}\tag{CI0}
\]

For positive integers $α,β$, the exact objective
$αh'+βF'=N$ becomes, with $K=αS+βP$ and $T=N+βF$,

\[
LK+\alpha Bh-\alpha De=T.
\]

Every possible nonempty tail must therefore satisfy

\[
\boxed{
\max\left(2,\left\lceil\frac{T+2\alpha h}{K+2\alpha h}\right\rceil\right)
\le L\le
\left\lfloor\frac{T+\alpha e}{K}\right\rfloor.}\tag{CI2}
\]

The left inequality uses the upper bound on $B$, and the right one uses
$B\ge0,D\le1$. Conversely each integer $L\ge2$ specifies exactly one
admissible tail: write $L-1$ in binary, least significant bit first. Evaluate
its $B$, then check CI1. Together with CI0 this is an exact terminating inverse
over **all tail lengths**, not a sampled length cutoff.

For a zero-ready transitive prefix of any size $s\ge1$, take
$α=1,β=2^s$. It suffices to inspect clocks satisfying

\[
h+2^s(P-F)\le N.\tag{CI3}
\]

Indeed the new graph contains $U$, so deterministic Hamilton-path insertion
gives $h'\ge h$, while CI0/CI1 give $F'\ge P-F\ge h\ge1$. The last
inequality follows because the $h$ covers whose first bin is empty are not
counted by $F$. The clock search therefore terminates. No positivity or
strongness assumption is silently inferred from a scalar solution: literal
funnels establish graph existence, and the boundary filter above establishes
the claimed strong universe.

The map is from an actual marked core and tail word to an actual marked
component, preserving the old induced core and its physical names. The scalar
inverse forgets the rest of the attachment response. Its necessary sidecar for
recursive use is the full polarized path-bin data and the literal named graph;
the five fields in CI1 only suffice for this one objective. The cheapest
decisive test is the interval CI2, followed by a literal graph reconstruction.

## Recovering the actual named source

The source $B$ has names
`[1,2,4,6,7,8,14,21,22,31,32,33,40,54]` and rows
`[11134,632,2682,112,512,80,528,9087,126,296,13311,10234,11263,894]`.
Its exact profile at 31 is $(h,s,q,F)=(2195,1,53,42153)$. Its unique
Hamilton path from31 is
`31,22,4,33,32,40,21,1,54,2,6,8,14,7`.
Its rank in the inherited order 14 unique-start rooted grammar is 8868.

The canonical successive cores are as follows. The receipt retains every
adjacency row and reconstruction permutation, including all named tail cuts.

| Current mark | Core names, in the recorded order | Entry | Tail names | Tail bits |
|---|---|---|---|---|
|31|1,21,32,33,4,22,40,54|22|2,6,8,14,7|10011|
|22|1,21,32,33,4,40,54|4|empty|empty|
|4|1,21,32,33,40|33|54|1|
|33|32|32|40,21,1|101|

The corresponding $(h,e,P,F,S)$ are

| Order / mark | $h$ | $e$ | $P$ | $F$ | $S$ |
|---|---:|---:|---:|---:|---:|
|14 /31|2195|53|672920|42153|2673|
|8 /22|53|22|1640|487|62|
|7 /4|31|9|644|157|31|
|5 /33|9|1|84|11|11|

In particular the core 8, tail 10011, and clock 2 recover the actual named $B$,
not just its total170807. The graph reconstruction is checked before using it
as the positive control.

There is also an all-mark whole-component obstruction. For every physical
mark $v\in B$, swapping the two ordered bins gives $2E_v(2)\le P_B(2)$:
the events that the first bin ends at $v$ and that the second ends at $v$
are disjoint. Every funnel wrapping all of $B$, at any mark, has
$F'\ge336460$. Every outer clock $s\ge1$ therefore gives total at least
$2195+2\cdot336460=675115$, exceeding both desired immediate targets.

## Exact menus and their limits

First, allow every physical entry mark on each actual nested core of orders
14,8,7,5; every strong funnel tail of any length; and every outer TT clock
$s\ge1$. These are new component candidates; the arbitrary tail length may
add or discard old tail names, so these menus are not assertions of an
all-old-vertices extension. Neither512423 nor512425 occurs.

| Fixed actual core order | Nonempty integer candidates for512423 | For512425 |
|---|---:|---:|
|14|0|0|
|8|47|47|
|7|152|152|
|5|1740|1741|

Every permitted empty tail is separately checked. For target170807 the sole
hit in these menus is the actual core 8 at 22, tail 10011, clock 2. This is an
essential positive control: the inverse does recover a known nontrivial target.

Second, change at most one old arc **inside the actual core 8**, then permit
every physical mark, every strong funnel tail, and every outer clock. There
are 29 cores before the strong filter and 25 after it. CI2 tests 1077 and 1078
integers for the respective targets, again with no hit. This is a complete
statement for that one-arc core menu, not a limitation on arbitrary changes
deeper in the grammar.

Third, leave the order 14 component size fixed and change at most one old arc
of $B$, with every mark and outer clock. This is 92 actual component graphs,
of which86 are strong; neither target occurs. The closest below-target output
in this menu is 512225: reverse the old arc between names1 and31, use mark54,
and clock 4, giving $(h,F)=(6785,31590)$. An independent Python endpoint DP
checks this C++ witness. This scan does not subsume the funnel menus, nor do
they subsume it.

Finally, the sibling recovery lane's three-$N$ unit-F dial would supply
both next outputs from a unique-start source $D$ satisfying

\[
9H(D)+11E_w(D,2)=512426.
\]

Apply CI2 to a single funnel over every rooted unique-start strong core through
order 8, plus the singleton:233 classes in orders 1,3,4,5,6,7,8. This excludes
all tail lengths, after 8944 candidate integers and all valid empty-tail checks.
There is no solution in this universe. By the complete inherited unique-start
grammar, **any such $D$, if it exists, must have a first strong core of
order at least9 after deleting its unique-start mark**. This does not exclude
larger cores, other unit-dial formulas, or the separately successful target
components. It identifies an actual structural requirement on this proposed
shared branch source.

The missing induction ingredient is thus a constructive way to change the
marked strong core and transport its joint attachment data, or a larger closed
component grammar. Arbitrarily lengthening the present tails cannot supply
the immediate branches. The outer clock $h+2^sF$ remains useful, but it is
not by itself a sufficient interface for repeated return words.

## Reproduction and audit

From the repository root:

```bash
python3 05-knowledge/results/zero_ready_component_returns_20261103/component/check.py
python3 -O 05-knowledge/results/zero_ready_component_returns_20261103/component/check.py
```

The producer compiles `menu.cpp` into a temporary directory, replays its receipt,
and removes the executable on exit. `check.json` and `check.optimized.json`
are byte-identical PASS receipts with inherited-source hashes. They include
255 independent before/after assignment controls;128 literal funnel, strongness,
and inverse controls; the hostile singleton/empty boundary; named source
reconstruction; all menu sizes; and every successful inverse representation.
No Python `assert` is used. Scope excludes larger old-arc rewrites, multiple
nested fresh funnel layers, arbitrary strong seeds, and full spectrum coverage.
