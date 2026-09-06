# Joint overlap positions close an actual balanced decoder family

**Status: PROVED + FINITE-EXACT; [independent analytic/source audit and exact
replay passed](long_frontier2_sep06_lrc_audit.md).** The finite
scale certificate below closes a class of seven-tail translated grids and
81 actual, unitless, balanced decoder entries. The declared older individual
interval and scalar consumers fail throughout these entries. This is not a
universal LRC(14) claim or a claim that all other safety methods fail.

## 1. Inheritance and the retained coordinate

The closest proved mechanisms are the exact translated-grid count and
pointwise overlap credit in [the third grid note](third_20260906_grid.md),
and the individual open-component count in
[the continuing4 packet note](continuing4_20260906_lrc_packets.md).
The latter already closes its entire displayed scale class starting at
1334: although pooled overlap estimates can fail, the sum of the individual
component bounds succeeds. Reusing those controls would not give a new
residual family.

Here every individual open-component bound is zero. The missing coordinate
is the position of each component after multiplying by the grid size and
reducing modulo one. Their positions cannot all simultaneously realize
their separate worst counts. Two primitive pairs suffice to exploit this
fact on a declared 84-scale universe. The exact sweep itself is elementary
and is not presented as a new general algorithm or grid instrument.

The actual-entry inputs remain those of **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
and the full actual decoder normalization in
[the general decoder note](overnight14_20260906_lrc_general_decoder.md).
The signed-box and dual consumers are inherited from
[the exact box note](overnight12_20260906_lrc_gcd_semigroup.md),
[the native decoder note](open_frontier_sep06_decoder.md), and
[the dual-pair note](continuing1_20260906_lrc_dual_pair.md).
The automatic balanced closure in
[the third decoder note](third_20260906_decoder.md) applies when
`max U<=28`; our maximum is much larger. The third grid note closes
`t>=97097`; our scales are all below 791. Neither observation says these
older sufficient conditions are necessary.

The concept board has five entries: individual overlap lengths; their
joint projected positions; an adaptive choice between two pairs; actual
mixed-support exclusion; and the half-phase of an odd component. The
canonical hostile within this bundle is the primary pair's zero overlap
at the actual half-grid for five scales. The corrected near miss is to
mistake an unpooled length-bound failure for a joint-position failure. The
least-used sidecar is the common phase of the projected intervals, retained
exactly below. The map from physical tails to primitive pair ratios preserves
all-translates intersection minima because its common factor is coprime
to the grid size; it loses a prescribed phase, so the proof takes a minimum
over that phase. No unrelated divisor-bag ceiling is sharpened.

## 2. A complete finite-scale supplier for arbitrary remaining tails

Put

```
T = {t: 620<=t<=790, gcd(t,128*101*113)=1}.
```

This is a set of exactly **84 integers**, not a height bound on the tails.
For `t in T`, let `U` be seven distinct positive integers, each coprime to
`t`, containing, for some positive integer `h`, the triple

```
h*(11413,12928,14464) = h*(101*113,101*128,113*128).
```

**Theorem 1.** Every translated `t`-grid contains at least five points at
which all seven members of `U` have clearance at least `1/14`.

The other four tails and `h` may have arbitrary heights compatible with
the stated coprimality and distinctness. No decoder or physical-box
assumption is needed for this grid theorem. Its conclusion is weak safety.

Write `D_u={x mod 1: ||ux||<1/14}`. Since `gcd(t,u)=1`, multiplication by
`u` permutes a translated `t`-grid, so each of its seven danger counts is at
most `ceil(t/7)`. For any chosen pair, the union count is at most the sum of
these seven counts minus that pair's intersection count. Thus it is enough
to obtain an intersection count at least

```
E(t)+5,             E(t)=7*ceil(t/7)-t.
```

Let `kappa_(p,q)(t)` be the minimum danger-intersection count for the
primitive pair `(p,q)` over every translated `t`-grid, with the danger sets
open as defined above. The exact certificate establishes

```
max(kappa_(101,128)(t), kappa_(113,128)(t)) >= E(t)+5
```

for every `t in T`. The minimum surplus is exactly five. The first pair
alone fails to exceed `E(t)` exactly at

```
633, 667, 687, 741, 761,
```

and its minimum is zero at each. The second pair repairs every one.
The output records all 84 rows `[t,E,kappa_1,kappa_2]`.

The first and third members of the stated triple have common factor `113h`
and primitive ratio `(101,128)`; the first and second have common factor
`101h` and primitive ratio `(113,128)`. Each common factor is coprime to
`t`. Multiplication by it sends an arbitrary translated `t`-grid bijectively
to another translated `t`-grid. Consequently the displayed minimum applies
to each physical pair. We use the larger of two valid pair credits, not
their sum. This proves Theorem 1 from the finite certificate.

### Exact completeness of the phase certificate

For coprime `p<q`, express circle coordinates in units of `1/L`,
`L=14pq`. The literal open arcs for `D_p` have endpoints
`14kq-q,14kq+q`, and those of `D_q` have endpoints `14kp-p,14kp+p`.
Intersect their ordered interval lists and rejoin the artificial cut at
zero. This enumerates every open component of `D_p intersect D_q`.
An independent length formula checks their multiset:

```
2p once, and min(2p,p+q-14k) twice
for 1<=k<ceil((p+q)/14),
```

again in units `1/L`. There are `2ceil((p+q)/14)-1` components.

Every component length is at most `1/(7q)`. For both pairs and every scale
in the declared band, its length is strictly less than `1/t`. An open
component `(a/L,b/L)` therefore contains one point of the translated grid
`{(alpha+j)/t}` exactly when `alpha` lies in the open circle arc with
endpoints `ta mod L,tb mod L`, and otherwise contains no grid point.
Different original components contain different grid points. Hence the
intersection count is precisely the overlap multiplicity of these
projected arcs.

The program sweeps every projected endpoint in exact integer arithmetic.
At an endpoint it first removes arcs ending there, evaluates the count at
the endpoint with new starts still absent, and then adds arcs starting
there for the next open cell. It checks cyclic closure. This convention is
essential: weak-safe grid points include open-danger endpoints. Every
endpoint and every complementary open cell is represented, so the minimum
is over all real translates, not a sampled grid of translates.

As a separate path, the program literally enumerates all `t` grid points
at every critical phase and every open-cell midpoint for both pairs at
`t=621,633,645,687,761`. These counts agree with the sweep including all
endpoint counts. These ten full phase-cell controls supplement the complete
84-scale exact certificate; they are not its universal quantifier.

## 3. Eighty-one actual boxed equality entries

Fix `Q=91^6=567869252041` and define

```
q=(287,277,271,265,263), P=product(q),
V={P} union {69P/q_i},
r=(128,101,113,131,149,167,227), H=product(r),
U={H/r_i}.
```

In increasing order these are

```
V=(360993824985,374026093035,382307113545,
   390963123663,393936227265,1501525040155),
U=(4761938937472,6472815202432,7254766032256,
   8251604113024,8445001084423,9566018927488,10702575631744).
```

For every scale

```
620<=t<=790, gcd(t,H)=1,
```

let the physical row be `tV union U`. There are exactly **81 such rows**,
with `g=1`. They form a subset of Theorem 1's scale set. Here the required
triple has multiplier `h=H/(128*101*113)`.

**Theorem 2.** All 81 rows are primitive, positive, distinct physical rows
inside `sum speeds<Q^2`, have exactly two actual decoder components of
sizes six and seven, and satisfy `W_(Q,3)=V_dec`. Every inherited
seven-through-twelve-subset profile is retained. All 81 rows have a time
with clearance strictly greater than `1/14`.

The denominators within each displayed shape are pairwise coprime. Both
shapes are primitive and unitless, every `V` label is odd, and their prime
supports are disjoint. Since `gcd(t,H)=1`, every cross pair `(tv,u)` is
coprime. The whole physical box follows from
`790*sum(V)+sum(U)<Q^2`.

The five displayed star edges of `V` have primitive sums
`356,346,340,334,332`. A spanning tree of `U`, given by its denominator
labels, uses

```
128--167, 128--227, 227--113,
113--101, 113--149, 167--131,
```

with primitive sums `295,355,340,214,262,298`. These sums belong to the
actual strict atlas: each prime is `2 mod 3`, each exponent is at most two,
and the sum is at most 356. Every cross pair is coprime and its sum exceeds
356, so no cross edge exists. The source reconstructs the complete graphs,
not only these spanning edges.

There are 105 mixed supports with two `V` labels and one `U` label, and
126 with two `U` labels and one `V` label. Both orientations, including
support-two degenerations and every possible nonzero distinguished
coefficient, are excluded as follows. Exact shape extrema are

```
D_min=min gcd(v,w)=1303226805,
S_max=max (u+w)/gcd(u,w)=394,
max(floor(Q/D_min)+1,floor(394Q/min(V))+1)=620.
```

For a relation on `tv,tw,u`, put `D=gcd(v,w)`. The distinguished
coefficient of `u` must be a nonzero multiple of `tD`, because
`gcd(tD,u)=1`. Its magnitude exceeds `Q` for every selected scale.
For a relation on `u,w,tv`, put `D=gcd(u,w)`. Coprimality forces the
outside coefficient to be `jD`, where `j` is a nonzero integer. The
remaining equation requires

```
|j|tv <= Q*(u+w)/D <=394Q,
```

contrary to `t min(V)>394Q`. These arguments exclude higher distinguished
multiples as well as the minimal one. Every bounded support of size at
most three is therefore internal. The displayed internal spanning edges
have coefficient heights below `Q` and span the two internal weighted
kernels. They give the reverse containment and hence
`W_(Q,3)=V_dec`, of dimension `(6-1)+(7-1)=11`.

Every physical subset of size seven through twelve has gcd one. A mixed
subset contains a coprime cross pair; the only possible pure subset of
these sizes is all seven members of primitive `U`. Its complementary gcd
word is thus all ones. The source checks membership in the pinned complete
inherited profile bank, not only its numerical gcd ceilings. It additionally
enumerates all 4095 such subset profiles and all 231 mixed supports at
five named scales. The preceding arguments cover every one of the 81 rows.

Finally, `eta=1/2` is strictly safe for `V` because all its labels are odd.
The physical times `(eta+j)/t` preserve clearance `1/2` on `tV`. Theorem 1
supplies at least five weak-safe choices on `U`. Both `t` and the phase
numerator `2j+1` are odd. Six `U` labels are even, so their phase fractions
have odd reduced denominator; none can equal a weak boundary with reduced
denominator 14. Thus only the unique odd `U` label can be at a weak
boundary. If it is, move the phase slightly toward its good interior.
Every other label is already strict and stays strict for sufficiently
small motion. This proves the asserted strict safety. Independently, the
source's first literal safe grid point is already strict in all 81 rows;
it records the full-row rational phase and clearance.

## 4. A genuine failure and the older consumer boundary

At the five primary-pair failure scales the common factors of the selected
actual pairs are odd. The physical half-grid therefore maps to normalized
translate `alpha=1/2` for both pairs. Literal counts are

| t | primary intersection | secondary intersection |
|---:|---:|---:|
| 633 | 0 | 18 |
| 667 | 0 | 16 |
| 687 | 0 | 14 |
| 741 | 0 | 16 |
| 761 | 0 | 24 |

This is a hostile to using the first pair alone at an actual supplied
phase. The surviving theorem permits the second pair. It is not a hostile
to safety: every displayed row is strictly safe.

For every pair of `U`, its primitive larger coefficient `q` is at least
113, and each open intersection component has length at most `1/(7q)`.
Since `t<=790<791`, the individual interval bound `ceil(t*length)-1`
is zero for every component of every one of the 21 pairs. The common gcd
is coprime to `t`, so no lost orbit multiplicity modifies that conclusion.
Consequently all their summed individual bounds, weaker pooled bounds,
and forests made from those bounds fail to give a positive safe count.
The nested-origin credits vanish as well. This is stronger than the
continuing4 examples where the individual length bounds already succeeded.

The source also checks failure of these specified sufficient widths, in
the inherited balanced normalization `a=6,b=7,g=1,L=max U`:

* All 126 forward native widths `6R/(56Lv)<1`, with a pair `u<w` in `U`,
  `p=u/D,q=w/D`, `R=Q(p+q)-(p-1)(q-1)`, and outside `v in V`.
  Here every physical cross gcd `delta` is one.
* All 36 maximum-endpoint widths fail: `56Dv>6Q` for each `U` pair
  containing `L` and every outside `v`.
* All 105 dual native widths fail: `6Q/(56L gcd(v,w))<1`, for each
  pair in `V` and each distinguished `U` label.
* Both whole-safe-arc grid gates fail: `6t<56L` and `7<49 max(V)`.

The statements concern the width parts of those sufficient consumers,
so no target-containment test can repair them in these rows. We do not
assert failure of every possible packet, actual phase search, or future
joint-position theorem. The new supplier deliberately uses a sidecar those
listed quantities discard. In particular the family has now been closed;
it is a structural stopping object rather than a remaining unsafe family.

## 5. Exact artifact and limits

The standard-library source is
[long_frontier2_sep06_lrc.py](../../04-computation/long_frontier2_sep06_lrc.py),
with [the frozen transcript](long_frontier2_sep06_lrc.out). It imports no
producer. All checks remain active with `python -O`. Reproduce with

```
python3 -B 04-computation/long_frontier2_sep06_lrc.py
python3 -B -O 04-computation/long_frontier2_sep06_lrc.py
```

The only inherited data bank is
`lrc14_joint_shadow_empty_core_next_sep06.json`, SHA256
`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.
If a sparse checkout omits that tracked file, the source reads the same
bytes from `HEAD` and checks the pin. It does not mutate Git.

The declared universes are all 84 scales for both joint minima, all phase
cells for the ten independent literal controls, all 81 complete actual
graphs and physical half-grids, and the five full support/profile controls.
The family has 81 scale choices at fixed large coefficient heights; it is
not an infinite family inside the fixed physical box. The generic supplier
does allow arbitrary other tail heights but only the 84 stated grid sizes.
No claim about other scales, arbitrary balanced entries, or universal
LRC(14) follows from this certificate.

The normal and optimized outputs are byte-identical, with **49,880
always-active gates** and semantic SHA256
`54f257677f163e86e1f59dc5cfdaac6d65a4d00756474fb88b0481ceddcf76ed`.
Raw SHA256 pins are:

```
source 83ddf56931ae05e8d467ee130f11f862279433df3a649960ae1ae12019f8b614
output 90382eb5d6aad90b3f0988e9055764ca1238a7cb62107e79a0fe7b841d379dc5
```

The source and output are frozen. The independent audit accepts the full
proof, actual-entry types, open-endpoint handling and both reproduction
runs. Its only draft finding was the repaired prose identification of the
two physical pairs; no source or mathematical correction remained.
