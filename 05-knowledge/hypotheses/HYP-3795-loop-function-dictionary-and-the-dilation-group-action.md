---
id: HYP-3795
title: THE LRC WITNESS SEARCH IS A GROUP ACTION -- a dictionary of ~23 loop-functions (maps of R/Z) that operate group-like (AGL(1) affine + PSL2(Z) modular + PSL2(R) Verblunsky/Blaschke + the Gamma_0(14)/PSL2(7) apex), and its clean-next-step synthesis: the LRC witness search over rational times a/q IS the dilation group G=(Z/q)* acting diagonally on the residue configuration; a time a/q is lonely iff a.c lands in the SAFE BOX B_r=[rq,q-rq]^13; so M(S)=max_q (1/q) max_{a in (Z/q)*} min_i ||a v_i||_q, and the covering-min lower bound = every covering config's dilation-orbit meets the safe box at some modulus q (<= MSS bound 91^12, S78). The orbit-count #{a: a.c in B_r} is a RAMANUJAN character sum (kps-S7 c_q(k) = the group Fourier): count = phi(q)(safe frac)^13 + sum_{k!=0}(prod ghat(k_i)) c_q(sum k_i v_i). UNIFICATION: the orbit-count at the binding modulus = #extremal atoms = flat-extension rank (S76) = #binders (S77): construction=2, AP=6=phi(14).
status: CONFIRMED (dictionary group relations verified; the group-action identity M(S)=max_q(1/q)(orbit-best) and the Ramanujan character-sum expansion verified exactly at q=183; orbit-count=#atoms cross-checked). A SYNTHESIS + reformulation (builds on opus-S13's circle-map dictionary and kps-S7's Verblunsky/Ramanujan atoms), NOT a proof: the covering-min lower bound becomes "every covering orbit meets the safe box", and the analytic engine is the resonant Ramanujan error (= the singular series THM-501/515) which the construction's AP residues make LARGE (main term overestimates count 15 vs actual 2) -- the crux is unchanged, but now group-framed and unified across S76/S77/S78/kps/opus.
source: mac-mini-2026-06-30-S79
related:
  - HYP-3792   # S77 safe-band frame (the box B_r); the witness a is the group element
  - HYP-3793   # kps-S7 flat-extension atoms = units, moments = Ramanujan sums (= this orbit's Fourier)
  - HYP-3789   # S76 moment/flat-extension (orbit-count = #atoms = rank)
  - HYP-3794   # S78 MSS finiteness (q and |G|=phi(q) bounded)
  - HYP-3768   # the Dedekind sum = the PSL2(Z) cocycle of these loop-maps (the margin)
  - THM-501    # the singular series = the resonant Ramanujan error (the crux)
results:
  - 04-computation/loop_function_dictionary_and_group_action_macmini_20260701.py
  - 05-knowledge/results/loop_function_dictionary_and_group_action_macmini_20260701.out
---

# HYP-3795 -- the loop-function dictionary and the dilation-group action

The owner's seed -- "pushing Verblunsky to the unit circle is a recursive metaphor for the LRC (points on
a circle = runners on a loop); define creative functions between points on a loop, a whole dictionary,
group-like; and work the clean next step, synthesizing." Two deliverables.

## (I) The dictionary (~23 loop-functions, group-like) -- builds on opus-S13
Maps of the loop `R/Z` (and configuration maps), organized by the group each generates:
- **Isometry `O(2)`**: rotation `R_c` (observer shift), reflection `iota` (antipode / `t<->1-t`).
- **Affine `AGL(1,R/Z)`**: dilation `M_v` (a RUNNER), **unit dilation `M_a`, `a in (Z/q)*`** (the WITNESS
  search generator), affine `x->vx+a`, the `times-n` endomorphism.
- **Mobius / Verblunsky (`PSL2(R)` = disk-aut boundary)**: Blaschke `B_A` (Verblunsky atom-pusher,
  `|A|->1` = atomic = lonely), Cayley, Mobius `(ax+b)/(cx+d)`, Gauss `{1/x}` (CF / `PSL2(Z)`), Farey
  mediant (builds `t*=[0;n-1,n]`).
- **Flow / dynamics**: runner flow `x->x+vt`, three-gap first-return (interval exchange), winding number.
- **Configuration maps**: scale-all `{v_i}->{cv_i}`, insert/delete runner (Mode A), convolution /
  Minkowski sum, Atkin-Lehner `W_N` (apex `X_0(14)`), **Ramanujan-hat `c->(c_q(k))`** (the config Fourier).
- **Cocycles / measures on the group**: **Dedekind `s(a,q)`** (the `PSL2(Z)` cocycle = the LRC margin),
  Verblunsky `alpha_j` (`|alpha_j|=1/(n-1-j)` for the AP, opus-S13), sawtooth `((x))` (iota-odd),
  distance `||x||` (iota-even, the loneliness).

Verified group relations: `M_a M_b = M_{ab}`; affine law `A_{v,a}A_{w,b}=A_{vw,vb+a}`; `iota M_v iota =
M_v`; semidirect `M_v R_a = R_{va} M_v`. So rotation+reflection+dilation generate `AGL(1)`; Gauss/mediant
generate `PSL2(Z)`; Blaschke gives `PSL2(R)`; the arithmetic apex is `Gamma_0(14)/PSL2(7)` (Klein quartic).
**Runners are the dilation subgroup.**

## (II) The clean next step -- the LRC witness search IS the dilation-group action
A rational time `a/q` gives runner positions `a v_i / q`; loneliness `>= r` means every `||a v_i||_q >= rq`,
i.e. the dilated config `a.c` (`c = (v_i mod q)`) lands in the **safe box** `B_r = [rq, q-rq]^13`. Hence

    M(S) = max_q (1/q) * max_{a in (Z/q)*} min_i ||a v_i||_q,

and the **covering-min lower bound = every covering config's `(Z/q)*` dilation-orbit meets `B_r` at some
modulus `q`** (`<= 91^12` by MSS, S78). The safe box has relative volume `(1-2r)^13 = (6/7)^13 = 0.135`
(LRC) -- the orbit must hit a box of this positive density.

**The orbit-count is a Ramanujan character sum** (expand the safe-interval indicator in additive
characters; the inner `sum_{a in (Z/q)*} e(a m/q) = c_q(m)`, kps-S7's object):

    #{a in (Z/q)* : a.c in B_r} = phi(q)(safe frac)^13  +  sum_{k != 0} (prod_i ghat(k_i)) c_q(sum_i k_i v_i).

Verified at `q=183`: main `phi(183)(156/183)^13 = 15.06`, but the **exact count = 2**. The gap is the
resonant Ramanujan error: the construction's residues are a maximally-resonant AP `{14,28,...,168,169}`,
so `c_q(sum k_i v_i)` is large on the AP sublattice (`c_183(61)=-60`, `c_183(3)=-2`) and pulls `15 -> 2`.
This is exactly why the singular series (THM-501/515) is hard -- the AP/construction is the worst resonant
case -- now framed as a group orbit-count.

## The unification (one number)
The **orbit-count at the binding modulus = #extremal atoms (S76 flat-extension rank) = #binders (S77)**:
construction `= 2` (the atoms `{14/183, 169/183}`), AP `= 6 = phi(14)` (the units `(Z/14)*`). The
group-action count, the moment-matrix rank, and the binding runners are literally one integer.

## Clean next step (honest)
Prove `#{a : a.c in B_r} > 0` for every covering `c` at some modulus -- i.e. bound the resonant Ramanujan
error below the `phi(q)(6/7)^13` main term. For generic (non-AP-resonant) covering configs the error is
small and `count ~ 15 > 0` easily; the crux is the near-construction resonant configs, where `count` drops
to `2` but stays positive (the binders). This is the S77 bulk/construction dichotomy in orbit-count
language, and the analytic tool is kps-S7's Ramanujan sums -- the group Fourier of the witness search.
NOT a proof; a unifying reformulation that puts S76/S77/S78, kps, and opus in one group-theoretic frame.
