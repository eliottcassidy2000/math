# Selective core replacement beyond homogeneous inflation

Status: **PROVED** for the parameterized construction below; **FINITE-EXACT**
for the separately labeled menus. The Chan–Pak exact odd spectrum remains open.
This note does not assert that the entire `405·2^k+71` family is unit-ready.
The parent packet's independently audited complete inverse excludes its `k=17`
member from that interface.

## 1. Inheritance and the operation being changed

The inherited complete ready fibres and homogeneous-substitution obstruction
are in
[`multicore_metagraph_20261113/FAMILIES.md`](../../multicore_metagraph_20261113/FAMILIES.md).
In particular, no homogeneous substitution in either ready881 source produces
1691. The present operation deletes one nonport and inserts two selectively
attached vertices. It retains the physical ready root, its sole predecessor,
and every arc between retained vertices. Thus it changes precisely an
assumption of that obstruction.

The closest positive compiler is the nested one-start funnel in
[`session296_multicore_metagraph_20261113/cores/ENTRANCES.md`](../../session296_multicore_metagraph_20261113/cores/ENTRANCES.md),
with its exact six-coordinate implementation `verify_ready.py`.
The underused sidecar is its three-bin path count: the replacement couples
binary and ternary growth, rather than a single dyadic coefficient.
Incoming `changing_core_returns` independently supplies neutral repairs at881
and further finite ready addresses. The present881-to1691 retained-vertex
split is a different operation, changing the exact count.

Write `P_M(j)` for the number of ordered spanning `j`-tuples of directed paths
in M, with empty paths allowed. Put

\[
u=H(M),\quad v=P_M(2),\quad w=P_M(3).
\]

A rooted graph with indegree-one root r and predecessor p is called
*d-ready* here when `H(T-r)-#HP(T starting at r)=d`. For this root type its
ending-root two-path response F satisfies `F=H(T)+d`. This is the actual
rooted interface, not an unmarked count assertion.

## 2. Universal module construction

Start with the cyclic triangle of blocks

\[
\{x\}\longrightarrow\{a\}\longrightarrow M\longrightarrow\{x\},
\]

rooted at x. A funnel word replaces the active root by a fresh root y: every
old vertex dominates a new transitive tail, y beats only the old root in the
old graph, and its incidences with the tail are encoded by the word (`0`
means y beats that tail vertex). The old nonroot vertices beat y. A nonempty
word ends in1. Finally append r beaten only by the last active root.

Define:

- `Old(M)`: apply funnel words `1`, then `11`, then append r.
- `New(M)`: first replace M by the ordinal sum `M→{b}`; apply words `1`, then
  `01`, then append r.

Both graphs are strong for nonempty M. Their rooted states `(H,s,q,F)` are

\[
\begin{aligned}
\mathrm{Old}(M):\quad
&H=-7u+105v+52w,\\
&q=15u+15v,\quad s=14u+15v,\quad F=H+u;\\
\mathrm{New}(M):\quad
&H=u+166v+114w,\\
&q=17u+34v,\quad s=16u+34v,\quad F=H+u.
\end{aligned}
\tag{1}
\]

Thus both retain readiness d=u. In particular they are unit-ready exactly
when M is transitive. This does not confuse a primitive interface with a
unit-ready one.

**Proof of the count law.** Use the inherited closed funnel coordinates
`(H,E,S(2),P(2),E(2),S(3))`, where E and S are the ending/starting active-root
ordered-bin responses, with the distinguished path in the specified bin.
For the initial three-block triangle their values are

\[
(u+v,\ u,\ u+v,\ 2w+4v-2u,\ 2v-u,\ w+2v).
\tag{2}
\]

One direct verification of (2) is to enumerate the positions of the two
named exterior vertices among the ordered path bins. If `K_i(M)` counts
ordered *nonempty* i-path covers, then
`v=2u+K_2` and `w=3u+3K_2+K_3`. The six coefficient rows in the basis
`(K_1,K_2,K_3)` are respectively

\[
(3,1,0),\ (1,0,0),\ (3,1,0),\ (12,10,2),\ (3,2,0),\ (7,5,1).
\]

There can be no more than three runs of M because only the two exterior
vertices can separate them. These are exact word/bin counts for arbitrary M,
not interpolation from transitive examples. Apply the two inherited funnel
matrices. They give parent counts `15u+15v` and `17u+34v`; the ready-child
count is parent H plus parent E(2) minus u. For the second construction
`P_{M→1}(j)=jP_M(j)`, so substitute `(u,2v,3w)` before applying its matrices.
This gives (1). The predecessor's starting-path count remains u through all
funnels; deleting the newborn proves q-s=u and then F=H+u. ∎

## 3. Literal retained-vertex law and the exact881 escape

Take `M=TT_m`, m≥1, so `(u,v,w)=(1,2^m,3^m)`. The programs from a singleton
are `(1^m,1,11)` and `(1^(m+1),1,01)`. Then

\[
H_o(m)=105\,2^m+52\,3^m-7,
\qquad H_n(m)=166\,2^m+114\,3^m+1.
\tag{3}
\]

Label the old graph in increasing construction blocks as follows: the three
nested roots are0,1,2; the initial singleton is3; M has labels4 throughm+3;
the middle tail is m+4; the final tail is m+5,m+6; and r is m+7.
Delete m+5. In the new graph retain labels0 throughm+3 unchanged, and map

\[
m+4\mapsto m+5,\quad m+6\mapsto m+7,\quad r=m+7\mapsto m+8.
\]

The two inserted vertices are m+4 and m+6. Reading the funnel orientations
shows that every retained arc agrees. In particular predecessor0 and r are
the same physical vertices and the unique incoming root arc is unchanged.
The construction is a selective replacement, not a sequence of individually
count-neutral arc edits.

For m=2 the literal pair is

```
old: (514,388,457,499,487,455,387,257,1,510)
new: (1282,772,905,1011,999,967,903,771,512,1,1022)
```

and its rooted states are `(881,74,75,882)` and `(1691,152,153,1692)`.
This supplies the exact first dyadic transition `1691=2·881-71` while
retaining nine of ten old vertices.

The failure boundary is exact:

\[
H_n(m)-(2H_o(m)-71)=10\,3^m-44\,2^m+86.
\tag{4}
\]

It vanishes at m=2; m=3 already gives4. Accordingly (3) is an infinite
physical family, but it is not a uniform405 dyadic return.

## 4. Ternary completion gate and finite stopping results

A fixed transitive module giving `H_m=K·3^m-1` and constant readiness1 would
allow the literal module extension `m→m+1` to implement `H→3H+2` indefinitely.
The positive run expansion makes this a precise gate: at m=1 its first
three nonempty-run coefficients must be

\[
A_1=H_1,\quad A_2=H_1+1,\quad A_3=(H_1+1)/3,
\]

and every higher coefficient must vanish. A first-step equality alone does
not certify this gate.

The complete strong unique-start parent menu through order9 gives three
ready near-families:

| base H | funnel address | module label in append-root convention | exact H_m |
|---:|---|---:|---|
|431|`1,11,1`|0|`132·3^m+36·2^m-37`|
|551|`1,101,ε`|5|`174·3^m+30·2^m-31`|
|893|`1,11,11`|0|`272·3^m+78·2^m-79`|

Their first two module sizes obey `H_2=3H_1+2`; the third misses by72,60,156.
The missing all-module correction is respectively `12K_3`, `10K_3`,
`26K_3`: a trade invisible at module sizes1 and2, but visible at3. This
identifies a specific required response sector rather than a scalar fit.
The run coefficient beyond3 vanishes in each case; the inherited
positive-run no-gaps lemma turns this into an all-size polynomial identity.
The included menu script records all four run coefficients, with A₄=0,
and independently counts the parent-start response at sizes1 through4.
Therefore these three all-size identities and their constant readiness are
**PROVED** by the inherited positive-run certificate; the discovery menu
itself remains **FINITE-EXACT**. They are not dependencies of Sections2–3.

A focused second menu considered all strong unique-start parents through
order12 with ready counts431,551,893, and every single module avoiding the
ready root and predecessor. It tested63 eligible modules in eight fibres.
Only the same three first-step hits occurred; none met the third-size gate.
This is the stated bounded menu, not an all-order absence theorem.
The next desired counts3887,4967,8045 do have ready witnesses in the parent
packet, so the fixed-module failure is not a global target obstruction.

The separate two-/three-position synchronized inflation controls tested322
and3038 graphs from the inherited complete finite fibres and found no
proposed target. They are retained as stopping receipts, not a structural
obstruction to selective operations.

## 5. Reproduction and remaining target

From the repository root:

```
c++ -O3 -std=c++17 05-knowledge/results/target_directed_cores_20261112/constructors/count_profiles.cpp -o /tmp/target_directed_count_profiles
python3 05-knowledge/results/target_directed_cores_20261112/constructors/selective_split.py
python3 05-knowledge/results/target_directed_cores_20261112/constructors/pure_ternary_probe.py
python3 05-knowledge/results/target_directed_cores_20261112/constructors/ternary_fibres.py
```

`selective_split_receipt.json` includes symbolic matrix identities, five
literal transitive sizes, the retained maps, and nontransitive C3/Q4 hostile
controls. The arbitrary-module identities and physical map prove the
unbounded claim; finite controls only guard transcription.

Open: supply a target-preserving selective operation that changes the third
run sector while retaining readiness, or use a wider defect interface at the
all-order unit-ready hole53084231. The present result does not close the
Chan–Pak odd-count conjecture.

A positive wider-interface escape is proved in [DEFECT_TWO.md](DEFECT_TWO.md):
the exact unit-ready hole53084231 admits a defect-two source and an infinite
extension family preserving its physical root and predecessor.
