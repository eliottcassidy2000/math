---
id: THM-3382
title: "Independent superunit affine-tail proof of canonical reflected-residue closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  This is an
  independent alternative to THM-3355's weighted-horn closure.  On every one
  of the 649 upper-median bodies, every physical pair with intrinsically high
  reduced level ratio has overlap strictly greater than one fifth of the
  universal singleton-debt maximum.  The resulting five-edge Hunter tree,
  THM-3352, and the same-level graph redundantly close every positive canonical
  reflected-residue level assignment on all 3,003 bodies.  THM-3360 later
  strengthens the pair floor to 1/294.  This is not k=0, arbitrary six-drift
  k=1, other residue packets, or projected k=2,3, and it does not prove LRC(14).
source: root/lrc-math-2026-08-12
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3350-connected-low-full-tree-atlas-dense-closure-and-uniform-tail
  - THM-3352-connected-low-all-head-universal-physical-forest-closure
related:
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
  - THM-3355-disconnected-low-affine-tail-and-reflected-branch-closure
  - THM-3360-uniform-reflected-high-pair-floor-by-admissible-affine-tails
verification:
  - 04-computation/lrc14_disconnected_low_geometry_verify_20260812.py
  - 04-computation/lrc14_disconnected_large_ruler_floor_20260812.py
  - 04-computation/lrc14_disconnected_35_small_ruler_symbolic_20260812.py
  - 04-computation/lrc14_disconnected_non35_g4_tail_20260812.py
  - 04-computation/lrc14_disconnected_head263_exact_scan_20260812.c
  - 04-computation/lrc14_generalized_dirichlet_reduction_20260812.py
  - 04-computation/lrc14_affine_turn_band_reduction_20260812.py
  - 04-computation/lrc14_affine_Tlt1_superunit_continuum_20260812.cpp
  - 04-computation/lrc14_affine_Tlt1_superunit_continuum_independent_audit_20260812.py
  - 04-computation/lrc14_affine_Tlt1_superunit_unified_K_20260812.py
  - 04-computation/lrc14_disconnected_bridge264_698_exact_scan_20260812.c
  - 04-computation/lrc14_disconnected_bridge264_698_independent_audit_20260812.py
  - 04-computation/lrc14_disconnected_bridge264_454_reference_audit_20260812.py
  - 04-computation/lrc14_disconnected_bridge455_678_all_argmins_audit_20260812.py
  - 04-computation/lrc14_disconnected_bridge455_678_full_context_audit_20260812.py
  - 04-computation/lrc14_disconnected_bridge679_698_all_argmins_audit_20260812.py
  - 04-computation/lrc14_disconnected_bridge679_698_full_context_audit_20260812.py
  - 04-computation/lrc14_disconnected_reflected_branch_synthesis_thm3355.py
turn_band_script_sha256: 8a4d59f9eaae04f6c836c1d8081a7dbaf74fa51ed9551784324c8e17f48bbf71
turn_band_output_sha256: 83785ea2e3d39952a7ab50e1f90e01db56a47f033de32d3b651e1a7f52bfabd8
continuum_script_sha256: 2db0af3e0e2b678ee6b13d79dbaf4f6b8f91c8e8a83ee0a1d201de344abaf06a
continuum_output_sha256: af1c6c944ea1e5997ed58809e4f8ddccde16326d4b90b39e3b777f4200498661
continuum_audit_sha256: 3d10014f093b2cbee89ba524fadd2e9a5b4b1dee92d28bc550da0869a18a2aa0
continuum_audit_output_sha256: 62a5eae88a76216fece5de8844f84a147b65392d8a77474e4cc88174096ab7c6
unified_K_script_sha256: db8575364318a66a8635ea9acebc80d1537520eaedfb880045a906017af576fb
unified_K_output_sha256: 427cec0fe83780fb3115fa0c4ccfc8d21d82c6d645d692122bc9abc97375bda7
bridge_script_sha256: f5981fb354c002d68eea6a589f015c76d124cb4d0fdf7e6fe6088abcf5423046
bridge_output_sha256: 532bedbb84d49e9dfb655bf66877483740a02259a1044a40eaada9fadf7e12ba
bridge_audit_sha256: d95423e6e7b710e155bbfefcfaec85a7168ef251824a5c9c162a1cc5400c6355
bridge_audit_output_sha256: 20afda5eb6592e931942334f3ab123e3f7e725ff88f3d8b5b387360cd10a98b4
bridge_264_454_reference_audit_sha256: 1438a1265b8a8ae412094b84757739e028d699cbfebdb339de99c75a4e94b7a7
bridge_264_454_reference_output_sha256: c5587d320666984856500960421f56ae0b808538254c5d18f9c7a2ab45dad741
bridge_455_678_argmin_audit_sha256: faacadca6b117f7a5cf0c88dca418a0870471037f5f0c8cfc5c1bfed3ed1346d
bridge_455_678_argmin_output_sha256: de672acb627da86bbef0e322915384742723daed408db32ea5f46debcf9031ca
bridge_455_678_control_audit_sha256: 976accfb57c0da5531cdbd23d2c6075b95611e30db7cdf4bf613d223a6a83216
bridge_455_678_control_output_sha256: 607a22b5d19e373d82ef0bb24540d745d27f2a481fd4c66598e7e7865297eb65
bridge_679_698_argmin_audit_sha256: c694bacf74c42790ab74b9157ba93ad3fc162f21b47b3de661b83a6972380196
bridge_679_698_argmin_output_sha256: bd717c2c2d69c8968043c873e8262a444390979810a8a415cac9fe38d00f74e5
bridge_679_698_control_audit_sha256: 16689e833f9894a51c5c560f50b1e70d4737d16c3080b13a1776bec7d2c2ddd3
bridge_679_698_control_output_sha256: 841aaaf4b0b0bbf37c6955ffb471bd865e76ab1a56c95ae6828c29e2fbfc917d
synthesis_sha256: 8bbc4817b75732e38935a3431902d6c57886d8a3a8bb0fea4c033cea4870dbb7
synthesis_output_sha256: ced9ee2e128290939fce706968631b20dc15e6cda966e327b22c1bd16e2fec89
hash_basis: LF-normalized bytes
---

# THM-3382 -- independent superunit affine tail and canonical reflected-residue closure

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The computation filenames and the frozen synthesis transcript retain the
provisional token `thm3355`: this proof and the independently developed
weighted-horn theorem were assigned that ID concurrently.  The weighted-horn
theorem keeps canonical ID THM-3355; this independent proof is THM-3382.

## 1. Statement

Let `E` be any of the `649` upper-median bodies behind THM-3350/3352,
let `(L,j)` be its fixed upper-median body-safe cell, and let `e!=f` be two
ordered endpoint labels in `E`.  For distinct positive raw levels `p,q`,
swap the two endpoints together if necessary so that `p<q`, and write

```text
g=gcd(p,q),       p=gP, q=gQ,       gcd(P,Q)=1.          (1)
```

If the reduced channel is intrinsically high,

```text
P+Q>=8,                                                (2)
```

then its exact reflected physical overlap satisfies

```text
I_(L,j,e,f)(p,q) > D_max/5,                            (3)

D_max=186636088362/11773143757375,
D_max/5=186636088362/58865718786875.                   (4)
```

The inequality is uniform over every feasible body and its fixed
upper-median safe cell, every endpoint ordering, and every positive pair of
levels.  It is a physical located statement, not merely a homogeneous or
projective limit.

Consequently:

1. every six-distinct-level assignment whose intrinsic low graph is
   disconnected closes by a five-edge high Hunter tree;
2. THM-3352 supplies the complementary connected-low closure;
3. every repeated-level assignment closes on the same `649` bodies by the
   same-level graph from THM-2941; and
4. every positive canonical reflected-residue assignment

   ```text
   z_e=q_e L-e,       q_e>=1,                            (5)
   ```

   closes on all `649` bodies.  Together with the `2,354` complementary
   bodies already closed in THM-2941, this exhausts all `3,003` bodies in
   the canonical reflected-residue branch.

## 2. Inherited geometry and the exact target

THM-3350 gives exactly `649` active upper-median bodies and `4,044` distinct
ordered physical contexts `(L,j,e,f)`.  The exact context split is

```text
2,530 contexts with L<4592,
1,514 contexts with L>=4592,                            (6)
```

with maximum small ruler `4,368` and minimum large ruler `4,620`.
Both ordered endpoint lanes are retained, so the orientation in `(1)` loses
no physical case.

For six distinct levels, join two vertices when their reduced ratio is low,
meaning `P+Q<=7`.  If this graph is disconnected with components
`V_1,...,V_k`, its high complement is complete multipartite.  Choose
`x in V_1` and `y in V_2`; join `x` to every vertex outside `V_1`, then
join `y` to every vertex of `V_1\{x}`.  These are

```text
(6-|V_1|)+(|V_1|-1)=5                                 (7)
```

cross-component edges and form a spanning tree.  The exact disconnected
component-shape atlas has `36,520` realizable unlabelled profiles; `(7)` is
structural and does not depend on enumerating them.

For any concrete six-level assignment, the singleton excess above `6/7` is
at most `D_max` by the rank reduction and exact `649*720` census of THM-3350.
Therefore the strict per-edge target `(3)` is exactly what is needed: five
tree edges contribute more than `D_max` in total.

## 3. Exhaustive pair-floor partition

After the simultaneous orientation `p<q`, the following cases are disjoint
and exhaustive.

### 3.1 Far ratios, all rulers

On a body-safe cell the first clause has exactly `p` full teeth.  The
mean-zero primitive of the `q`-clause indicator has oscillation
`6L/[49(qL-f)]`, so interval-by-interval integration gives

```text
I(p,q) >= pL/[49(pL-e)]-6pL/[49(qL-f)].                (8)
```

Whenever `q>=8p`, the right side is at least

```text
23/4655 > 1/294 > D_max/5.                             (9)
```

This case is taken first and applies to every ruler.

### 3.2 Moderate ratios on large rulers

Suppose `q<8p` and `L>=4592`.  The frozen large-ruler theorem proves

```text
I(p,q)>=1/294                                          (10)
```

for every primitive high channel and every common dilation.  Its analytic
envelope handles `P>=4`; the nineteen smaller primitive channels have
`28,766` exact physical controls, all positive.

### 3.3 Moderate ratios on the 29 small rulers

It remains to take `q<8p` and `L<4592`.  There are exactly the `2,530`
ordered contexts in the first row of `(6)`.

- On the primitive channel `(P,Q)=(3,5)`, an exact affine-residue compiler
  proves `I(3g,5g)>=1/294` for every `g>=1`.  It checks `56,191` residue
  classes; the largest stabilization point is `134`.
- Off `3:5`, the linked phase floor plus the all-channel midpoint estimate
  proves `I>=1/294` for every `g>=4`.  Only `118` channel-scale rows survive
  its envelope, and all `118*2530=298,540` literal exact masses pass.
- Thus only `g<=3` remains.  The exact integer-only head compiler closes
  every raw pair `p<264`: `201,377` channels and `509,483,810` physical
  masses, with zero failures and every reported argmin independently
  replayed.
- The finite bridge below closes every `264<=p<=698`.
- The affine argument below closes every `p>=699`.

This leaves no endpoint, orientation, common-scale, ruler, or ratio case.

## 4. Refined Dirichlet turn bands

For the remaining `g<=3`, small-ruler, moderate-ratio bank, Dirichlet gives

```text
1<=d<=8, 0<=a<=7d, D=d+a,
c=dq-Dp, |c|<=p/9,
A=Lc+eD-df,
r=p mod d, N=(p-r)/d,
z=Lp-e, w=Lq-f,
rho=|A|/z, T=N rho.                                    (11)
```

Grouping the first `dN` teeth into `d`-blocks gives a periodic piecewise
linear function `G_d` with

```text
mean(G_d)=dL/(49z),
osc(G_d)<=dL/(7w),
Var(G_d')<=4dL/w.                                      (12)
```

For `m<=T<m+1`, the composite trapezoid identity and `(12)` give the hostile
envelope

```text
B_m(p,L)=m(p-7)/[49(m+1)p]
          -8L/[14(L(p+1)-14)]
          -64(m+1)^2 L/[2(L(p+1)-14)(p-7)].            (13)
```

Exact positive-coefficient derivative certificates show that `(13)`
increases with `p` and `L` on the declared ranges.  Its exact hostile corners
prove:

```text
1<=T<4:  p>=264,
4<=T<5:  p>=277,
T>=5:    p>=264 by the inherited many-turn estimate.   (14)
```

If `T<1` and `c!=0`, the same inequalities sharpen the resonance to

```text
1<=|c|<=13.                                            (15)
```

For completeness, the inherited `T>=5` lower bound is explicit.  Put

```text
R(p,L)=(Lp/9+888)/(Lp-14),
B(p,L)=4(p-7)/(245p)
       -4L/[7(L(p+1)-14)]
       -L[pR(p,L)^2+8R(p,L)]/[2(L(p+1)-14)].           (15a)
```

The exact generalized-Dirichlet verifier proves that `B` increases in both
variables on `p>=264,L>=168`: after the hostile-corner shift, the two
derivative numerators have respectively `36` and `20` positive
coefficients.  At the corner,

```text
B(264,168)-D_max/5
=85330033783953387991/7166476998347435667648750>0.     (15b)
```

Thus the last line of `(14)` is a proved analytic inequality, not an unnamed
tail assumption.  The turn-band verifier proves every other rational corner
and both monotonicity claims without floating point.

## 5. The centered superunit tail

The limiting `d`-block convolution must be centered.  Put

```text
T_infinity=|A|/(Ld).                                   (16)
```

Finite `T<1` does **not** imply `T_infinity<=1`.  It implies only

```text
T_infinity < (p-e/L)/(p-r) < p/(p-r)
           <=679/672=97/96                             (17)
```

for `p>=679`, since `r<=d-1<=7`.

The exact centered continuum compiler evaluates every compatible tuple

```text
1<=d<=8, 0<=a<=7d, 1<=|c|<=13, gcd(a,d)|c,
all 2,530 contexts, |A|/(Ld)<=679/672,
a=0 => c>0,       a=7d => c<0.                        (18)
```

It performs `5,053,047` integer cross-product checks, including `362,926`
rows with `T_infinity>1`, `11,314` equality rows, and `29` zero-step rows.
The exact minimum is

```text
J_0=709/48048,                                         (19)
```

with two equality rows, one at

```text
(d,a,c;L,j,e,f)=(3,8,-1;168,90,12,1), A=-39.          (20)
```

An independent Python `PeriodicPL` implementation reconstructs the centered
periodic convolution from its breakpoints and checks `1,053` exact values:
`256` ordinary, `768` superunit, and all `29` zero-step rows.  It agrees with
the C++ compiler on every value.  Optimized, unoptimized, and sanitizer builds
of the compiler have byte-identical output.

For the finite-to-continuum comparison put

```text
zeta=3167/3168, omega=3155/3168, k_f=6864/22085,
delta_d=d-1+1/12,
tau_d=1+delta_d/(264 zeta),
KP_d=d k_f+(97/96)[d(d-1)/(2zeta)+d/zeta].             (21)
```

Both phase-drift terms in `KP_d` require the superunit factor `97/96` from
`(17)`.  The Peano term does not: the finite hypothesis `T<1` bounds it
directly by a smaller quantity.  Its physical error is at most

```text
|A|Ld/(2ZWp^2).                                        (21a)
```

Equivalently, after extracting the global `1/p`, its coefficient is at most
`|A|Ld/(2ZWp)<d^2/[2omega(p-r)]<=d^2/(1344omega)`, below
the conservative coefficient retained next.  A safe unified error constant is

```text
K_d^*=tau_d KP_d/d
      +2 tau_d delta_d/(7zeta omega)
      +tau_d d/(42omega)+d/(14omega)
      +d^2/(264zeta omega).                            (22)
```

It increases for `1<=d<=8`, and exact arithmetic gives

```text
K_8^*=1792138785426/221510098565,
K_8^*/(709/48048-D_max/5) in (698,699).                (23)
```

Thus every nonzero-step `T<1` point with `p>=699` satisfies

```text
I(p)>=J_0-K_8^*/p>D_max/5.                             (24)
```

The exact margin at `p=699` is

```text
46135211197112901571553/
4166631528336272997191190000>0.                        (25)
```

If `A=0`, complete `d`-blocks repeat and the separate zero-step constant is

```text
KP0_d=d k_f+d(d-1)/(2zeta)+d/zeta,
K0_d=KP0_d/d+(d-1)2/(7omega).                          (26)
```

Its maximum is `477044832/69943195`, so the `A=0` face already closes for
`p>=589`.  Finally, `c=0` would imply `P|d<=8`, and `g<=3` would force
`p=gP<=24`; hence it cannot occur in this tail.

## 6. Exact finite bridge and independent audits

The combined integer-only C scanner covers

```text
264<=p<=698, p<q<8p, gcd(p,q)<=3,
all 2,530 small-ruler ordered contexts.                 (27)
```

Its exact universe and result are

```text
1,211,966 raw channels,
3,066,273,980 physical mass comparisons,
0 failures.                                            (28)
```

The global weakest row is

```text
(p,q)=(698,2559), (L,j,e,f)=(168,90,12,1),
I=20682154/1400220127,
I-D_max/5=
956138253921819066776/82424964235704398433125>0.       (29)
```

The omitted `49.5 MB` ledger has SHA-256

```text
932057abf10f674e4bb31f334c1ea94f39e4e17627c6a950bd5d727f8e595186. (30)
```

It is byte-identical to the concatenation of three separately generated
segments:

```text
264..454: 397,502 rows, ledger 2865e79a...7053,
455..678: 734,566 rows, ledger f12b9cc8...908a,
679..698:  79,898 rows, ledger 0915c0c3...5421.         (31)
```

The repository preserves each segment scanner, frozen output, and reference
audit.  The slower THM-3352 reference engine re-evaluates every reported
argmin in all three segments: `397,502+734,566+79,898=1,211,966` exact
independent checks, with zero mismatches.  Full `2,530`-context minima are
also recomputed on `220`, `40`, and `20` deterministic hostile/spread
controls, respectively.  The combined assembly audit checks `(30)`, every
row of `(31)`, the exact task universe, strict target, and `(29)`.

The audit scripts use the canonical reference engine and context bank.  The
only omitted inputs are the regenerated ledgers themselves: the three segment
scanners write the exact filenames expected by the audits.  No promoted audit
depends on the superseded `455..823` full ledger.

Together, Sections 3--6 are an exhaustive proof of `(3)`: the analytic tail
starts at `699`, and the independently audited finite bridge ends at `698`.

## 7. Hunter closure and repeated levels

For six distinct levels with disconnected low graph, freeze the five-edge
cross-component tree from `(7)`.  By `(3)` its physical credit is strictly
greater than `D_max`.  Hunter's tree inequality therefore gives

```text
mu(union_e A_e) <= 6/7 + epsilon(E,q) - credit
                 < 6/7,                                (32)
```

because `epsilon(E,q)<=D_max`.  This closes the disconnected-low case.
THM-3352 closes the connected-low case, so every six-distinct-level
assignment closes on all `649` bodies.

For repeated levels, THM-2941's universal signed same-level graph is `K_6`
on `3,001` of the `3,003` bodies.  Its only two exceptions are

```text
(1,2,7,9,11,13), (2,4,7,9,11,13).                     (33)
```

Both have all `15` robust edges.  The active `649`-body universe is exactly
the robust-edge-`<=10` bank, so its intersection with `(33)` is empty.
Therefore any repeated level supplies a closing same-level pair on every
active body.  The synthesis companion independently reconstructs all `649`
bodies, `(33)`, the disjointness, the exact debt/tree inequalities, and
`2354+649=3003`.

## 8. Corrected near misses

Three tempting compressions were false or incomplete and are not used here.

1. The first continuum compiler integrated `[x,x+lambda]`; the physical
   convolution is centered at `x`, hence uses `x-lambda/2`.  The corrected
   exact minimum is `(19)`, not the earlier provisional `13/1022`.
2. The implication `T<1 => T_infinity<=1` is false.  A literal witness is

   ```text
   d=2,a=0,c=2,p=781,q=782,
   (L,j,e,f)=(784,420,2,1), A=1570,
   T=612300/612302<1, T_infinity=1570/1568>1.           (34)
   ```

   This is why `(18)` includes the superunit chamber.
3. In the first superunit repair, only one of two phase-drift losses carried
   the `97/96` factor.  Both do.  Correcting the second changes the rigorous
   integer tail start from the provisional `695` to `699`; the bridge was
   extended through `698` before promotion.

These errors affect only superseded scratch estimates.  The centered census,
the constants `(21)--(26)`, and the exact bridge are mutually dovetailed and
independently audited.

## 9. Scope and nonconsequences

This theorem closes the positive **canonical reflected-residue level branch**
`z_e=q_eL-e` on all `3,003` six-body carriers, including the former `561`
body residual.  It does not prove that an arbitrary six-drift `k=1` packet
has that residue form.  It does not close other residue packets, `k=0`, the
finite-but-uncensused projected `k=2,3` sectors, the seven-tail rung, physical
entry, or LRC(14).  LRC(14) remains open.

## 10. Reproduction

From the repository root, the short exact companions replay with

```bash
python3 04-computation/lrc14_affine_turn_band_reduction_20260812.py
python3 -O 04-computation/lrc14_affine_turn_band_reduction_20260812.py
python3 04-computation/lrc14_affine_Tlt1_superunit_unified_K_20260812.py
python3 -O 04-computation/lrc14_affine_Tlt1_superunit_unified_K_20260812.py
python3 04-computation/lrc14_disconnected_reflected_branch_synthesis_thm3355.py
python3 -O 04-computation/lrc14_disconnected_reflected_branch_synthesis_thm3355.py
```

The centered continuum census is

```bash
clang++ -std=c++17 -O3 -Wall -Wextra -Wconversion -Wshadow -Werror \
  04-computation/lrc14_affine_Tlt1_superunit_continuum_20260812.cpp \
  -o /tmp/lrc14-affine-superunit
git show HEAD:04-computation/lrc14_disconnected_head263_contexts_20260812.txt \
  | /tmp/lrc14-affine-superunit
python3 04-computation/lrc14_affine_Tlt1_superunit_continuum_independent_audit_20260812.py
python3 -O 04-computation/lrc14_affine_Tlt1_superunit_continuum_independent_audit_20260812.py
```

The full finite bridge is intentionally expensive.  Compile the combined
scanner and pass the canonical `2,530`-context bank, an output-ledger path,
and a thread count:

```bash
clang -std=c11 -O3 -pthread -Wall -Wextra -Wconversion -Wshadow -Werror \
  04-computation/lrc14_disconnected_bridge264_698_exact_scan_20260812.c \
  -o /tmp/lrc14-bridge264-698
/tmp/lrc14-bridge264-698 \
  04-computation/lrc14_disconnected_head263_contexts_20260812.txt \
  /tmp/lrc14-bridge264-698.ledger 8
```

The segment sources and audit companions provide the independent path
described in `(31)`.  Ordinary and optimized Python replays are byte-identical;
the stored outputs are the declared frozen transcripts.

The combined assembly audit is invoked after those four ledgers exist:

```bash
python3 04-computation/lrc14_disconnected_bridge264_698_independent_audit_20260812.py \
  /tmp/lrc14-bridge264-698.ledger \
  /tmp/disconnected_bridge264_454_exact_scan.ledger \
  /tmp/disconnected_bridge455_678_exact_scan.ledger \
  /tmp/disconnected_bridge679_698_exact_scan.ledger
```

For the path-pinned segment audits, first prepare their frozen `/tmp` names:

```bash
cp 04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py \
  /tmp/canonical_reference_engine_thm3352.py
cp 04-computation/lrc14_disconnected_head263_contexts_20260812.txt \
  /tmp/disconnected_head263_contexts.txt
cp 04-computation/lrc14_disconnected_bridge264_454_exact_scan_20260812.c \
  /tmp/disconnected_bridge264_454_exact_scan.c
cp 04-computation/lrc14_disconnected_bridge455_678_exact_scan_20260812.c \
  /tmp/disconnected_bridge455_678_exact_scan.c
cp 04-computation/lrc14_disconnected_bridge679_698_exact_scan_20260812.c \
  /tmp/disconnected_bridge679_698_exact_scan.c
cp 05-knowledge/results/lrc14_disconnected_bridge264_454_exact_scan_20260812.out \
  /tmp/disconnected_bridge264_454_exact_scan.out
cp 05-knowledge/results/lrc14_disconnected_bridge455_678_exact_scan_20260812.out \
  /tmp/disconnected_bridge455_678_exact_scan.out
```

Compile and run each segment source with the context path, its expected
`/tmp/...ledger` output path, and a thread count.  Then run the five
`bridge...audit_20260812.py` companions.  The all-argmin replays are
deliberately slower than the integer scanners.
