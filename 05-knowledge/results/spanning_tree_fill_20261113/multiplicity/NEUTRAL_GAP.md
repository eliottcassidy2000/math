# Neutral singleton relocation gives an exact marked response law

Status: **PROVED** for every nonempty tournament core, every selected subset,
and every finite binary suffix. The bounded receipt is a transcription
control. The operation controls an exact *difference*; it does not supply an
arbitrary prescribed absolute Hamilton total.

## 1. Physical operation and coordinates

Let C be any nonempty tournament on a fixed named vertex set, let u be a new
named singleton, and let v₀→…→vₜ₋₁ be a transitive suffix. Compare the parents

\[
T_-=C\to u\to TT_t,\qquad T_+=u\to C\to TT_t.
\tag{1}
\]

Their Hamilton counts both equal h=H(C). The preparation reverses exactly
the arcs between u and C; all vertices and other arcs are retained. It is a
simultaneous neutral preparation, not a sequence of individually neutral arc
flips.

In both parents adjoin the same named newborn w with the same incident arcs:

- w beats exactly J⊆V(C) inside C; put K=V(C)\J;
- ε=1 means w→u, and ε=0 means u→w;
- b_i=1 means v_i→w, with i starting at0 on the **left** of the suffix.

Put `x=Σ_(i=0)^(t−1)b_i2^i`, β=x mod2. Thus `0≤x<2^t`; the empty suffix has
t=0,x=0. Leading high zero bits are allowed, and length is retained even
though the formulas below depend only on x. Call the two children B₋ and B₊.
Their w-stars agree literally, and all old vertices survive.

For A⊆V(C), write:

- S_A,E_A: full Hamilton paths of C starting, respectively ending, in A;
- G_A,F_A: ordered spanning covers(P,Q) of C with P **nonempty** starting,
  respectively ending, in A, and Q possibly empty;
- P₂: all ordered two-path covers of C, allowing either bin to be empty.

Empty A gives zero. Exact complement identities are

\[
S_V=E_V=h,\qquad G_V=F_V=P_2-h.
\tag{2}
\]

The full rooted state at w is `(H,s,q,Froot)`: s counts full paths starting at
w, q counts paths after deleting w, and Froot counts ordered covers(P,Q)
with P nonempty ending at w and Q possibly empty.

## 2. Exact difference and the complete four-coordinate change

Let Δ=H(B₊)−H(B₋). Then

\[
\boxed{
\Delta=
\begin{cases}
F_K+S_K+xG_K-\beta E_K,&\varepsilon=1,\\
-(x+1)G_J-(1-\beta)E_J,&\varepsilon=0.
\end{cases}}
\tag{3}
\]

The entire rooted change is

\[
\boxed{
(H,s,q,Froot)_{B_+}-(H,s,q,Froot)_{B_-}=
\begin{cases}
(\Delta,S_K,0,F_K),&\varepsilon=1,\\
(\Delta,-S_J,0,-F_J),&\varepsilon=0.
\end{cases}}
\tag{4}
\]

In particular the deletion coordinate q=h is preserved. More explicitly,

\[
\begin{aligned}
s_-&=S_J,&s_+&=\varepsilon h,\\
Froot_-&=(2x+2-\varepsilon)P_2-F_J,
&Froot_+&=2(x+1)P_2-2F_J-\varepsilon h.
\end{aligned}
\tag{5}
\]

The indegree at w, also unchanged, is
`|K|+(1−ε)+Σb_i`. Thus preserving the newborn star does not assert readiness
preservation. If this indegree is1, its defect changes from `S_K` to
`(1−ε)h`. In particular the ε=1 preparation can consume a positive defect
and land at defect0. More sharply, if h>1 and the newborn has indegree1,
the prepared source is never unit-ready: its defect is0 or h. This is not
a universal ready repair.

The sign is exact. If ε=1, Δ>0 iff K is nonempty; if ε=0, Δ<0 iff J is
nonempty. Otherwise Δ=0. Thus equality holds exactly when w sees C∪{u}
uniformly. Indeed `G_A,F_A>0` for nonempty A: singleton P plus a Hamilton
path on C\P is already a cover. Also F_A≥E_A. These facts prove the strict
signs in(3), including x=0. Every Δ is even by Rédei parity of the two actual
children. No strongness assumption is used. The nonempty-core assumption is needed
for the absolute formulas: if C is empty the preparation is the identity,
and those first-SCC formulas are not applied.

A useful sufficient strongness guard is: C strong, ε=1, J and K nonempty,
t≥1, and the last suffix bit1. Then both children are strong. The general
formula also includes weak parents/children and empty J or K; these cases
must not be silently discarded.

## 3. Proof with the joint term retained until cancellation

We use the inherited selected-entry funnel response, also consistent with
the independently derived [newborn matrix kernel](../matrix/PROOF.md).
If a new apex beats J in C and sees suffix word b, then

\[
H=L_2G_J+B_1h-D_1E_J-C_J,
\tag{6}
\]

where C_J counts ordered two nonempty path covers whose first path ends in J
and second starts in J. The selected joint term cannot be dropped.
For the suffix of(1), `L₂=x+1=:L`, `D₁=1−β=:D`, and `B₁=g(x)`, with
`g(0)=g(1)=1`, `g(2a)=g(a)+a+(a mod2)`, `g(2a+1)=g(a)`.
These are the exact inherited tail coefficients from
[`two_entry_funnels_20261103/research.py`](../../two_entry_funnels_20261103/research.py).

For B₋, prepend the tail bit1−ε. Its coefficients become

\[
L'_2=2L-\varepsilon,\quad D'_1=\varepsilon,
\quad B'_1=g(x)+\varepsilon(L-D).
\]

For B₊, regard u→C as the core, with selected set J together with u when
ε=1. Its data are

\[
(h,\ E_J,\ G_J+\varepsilon P_2,
\ C_J+\varepsilon(F_J+S_J)).
\tag{7}
\]

To verify the last coordinate, place u in the first or second nonempty bin.
If the first bin contains core vertices, one obtains C_J. If the first bin
is u alone, its ending selection contributes εS_J. If u is in the second
bin, its starting selection contributes εF_J. The G-coordinate similarly
splits according to which bin contains u. Apply(6) on both sides; the joint
term cancels, and(2) simplifies the difference to(3).

For(4)–(5), a path starting at w must enter the first SCC of the deletion:
for B₋ this means a full C-path starting in J, whereas for B₊ it must start
with u, possible exactly when ε=1, then take any full path of C. This argument
still works when C is weak, because its full starting counts already encode
its own SCC order. Deleting w gives h on both sides.

Finally, for any parent T and cut X={v:w→v}, the ending-root cover count is
`Froot=P_T(2)−F_X(T)`. This follows by splitting off the last arc into w,
including the case P={w}. In(1), both parents have two-bin count
`2^(t+1)P₂`. Moving u changes the selected ending response by
`F_J+ε(h−P₂)`; contributions from selected suffix vertices are identical.
This proves ΔFroot in(4). Summing the suffix contributions
`Σ_(b_i=0)2^(i+1)P₂` gives(5). ∎

## 3a. Positive handoff to the retained singleton port

The newborn w is not the only available named port. For ε=1, the retained
vertex u has indegree1 in B₊, with sole predecessor w. Put

\[
q_u=(x+1)G_J+g(x)h-(1-\beta)E_J-C_J.
\]

Deleting u leaves precisely the selected-entry funnel with core C, apex w
and the original suffix. Its full paths starting at w are S_J. Since u beats
every other vertex except w, prefixing u to a deletion path is legal exactly
when that path does not start at w. Hence the **actual retained-u state** is

\[
\boxed{(H,s,q,F)_{u,B_+}
=(H(B_+),q_u-S_J,q_u,H(B_+)+S_J).}
\tag{8a}
\]

Therefore **S_J=1 gives a unit-ready port at u for every suffix**. Both u and
w are named from their creation and can be retained as a two-port interface;
selecting u does not resurrect a forgotten vertex. The preparation preserves
the parent count h, the newborn uses its fixed cut, and the resulting source
has an accepted unit-ready port for subsequent NR. Its next count is
`3H(B₊)+2`. This does not assert that the next count returns to this family.
The ready-source existence is the inherited unique-start funnel compiler;
the new contribution here is the explicit neutral relocation and marked
handoff connecting its two same-cut observations.

For C=C₃,J={0},ε=1 the formula becomes

\[
H(B_+)=15x+3g(x)+\beta+10,\qquad
q_u=3x+3g(x)+\beta+2.
\tag{8b}
\]

At x=1, the retained-u state is `(29,8,9,30)`, with u=3 and predecessor w=5
in the literal six-vertex rows. Thus the17→29 example lands at a genuine
unit-ready interface after choosing the already retained u, although its
newborn w has indegree3. For all-one suffixes x=2^a−1, a≥1, this specializes
to the familiar ready family `H=15·2^a−1`, `q=3·2^a+3`; that scalar family is
already present in the earlier H29 module analysis, and is not claimed new.

## 4. An exact two-phase difference inverse

For ε=1 and nonempty K, let `a=F_K+S_K`, `g=G_K>0`, `e=E_K`. Then

\[
\Delta(2j)=a+2jg,\qquad
\Delta(2j+1)=a+g-e+2jg\qquad(j\ge0).
\tag{8}
\]

A specified positive difference D is achievable in this fixed core/entry
chart iff one of the two arithmetic progressions in(8) contains D, with
`x=2j+β<2^t` if suffix length t is fixed. There are at most two candidate
addresses. Strongness may impose the additional final-bit and entry guards
above. For ε=0, the positive magnitude has progressions
`G_J+E_J+2jG_J` and `2G_J+2jG_J`, unless J is empty and the change is zero.

This is a complete inverse for the **difference**, not for an arbitrary
absolute target. For example the absolute unprepared count also uses C_J:

\[
H(B_-)=(2x+2-\varepsilon)G_J+
[g(x)+\varepsilon(x+\beta)]h-\varepsilon E_J-C_J.
\tag{9}
\]

To hit a prescribed N after preparation, one must solve(9)+Δ=N with the
same C,J,ε,x and physical guards. A solved difference cannot be added to an
unrelated witness's Hamilton count.

For C=C₃,J={0},ε=1, the complement moments are
`F_K=G_K=6`, `S_K=E_K=2`. Thus Δ(2j)=8+12j and Δ(2j+1)=12+12j.
At x=1 the states are

\[
(17,1,3,33)\longrightarrow(29,3,3,39),
\]

recovering the literal same-cut example in [BIRTH_FLOW.md](BIRTH_FLOW.md).
The newborn has indegree3, so these states are not unit-ready states.

## 5. Even the full named-subset H polynomial loses this preparation

For a tournament T on named vertices, define

\[
\mathcal P_T(z)=\sum_{S\subseteq V(T)}H(T[S])\prod_{v\in S}z_v,
\qquad H(\varnothing)=1.
\]

Both parents in(1) have exactly

\[
\mathcal P_{T_-}=\mathcal P_{T_+}
=\mathcal P_C(z_C)(1+z_u)\prod_i(1+z_{v_i}).
\tag{10}
\]

Indeed the induced graph on each named subset is an ordinal concatenation,
and its Hamilton count is the product over those factors. Therefore the
preparation preserves **every named induced-subset Hamilton count**, and in
particular every unmarked ordered-slot polynomial. Yet the same newborn cut
has different H by(3). Endpoint and selected-cut information cannot be
recovered from this entire induced-H polynomial alone. This is the precise
missing coordinate, not merely loss from a scalar H quotient.

Reproduction:

```
python3 05-knowledge/results/spanning_tree_fill_20261113/multiplicity/neutral_gap.py
```

The bounded replay uses five specified nonempty cores (including weak and
strong cases), representative empty/full/proper selected subsets, both ε,
and seven suffixes including the empty suffix and trailing high zero bits.
Its308 cases independently count both complete rooted states and checks the strict
signs, equality boundary, evenness and difference inverse. Its finite scope
is explicit in `neutral_gap_receipt.json`; the proof above supplies all cores
and all suffix lengths.
