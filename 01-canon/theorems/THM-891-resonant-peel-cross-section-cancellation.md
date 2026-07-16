---
id: THM-891
title: Exact cross-section cancellation for a resonant far peel
status: PROVED. The far-peel limit, exact miss-pattern formula, seven-residue reduction, consecutive-core constant, and pair-sector ray law are rigorous. THM-910 now closes the former negative residue-six gap, so every limiting residue satisfies the universal 0.097 propagation bound. The sharp 16/343 extremal and finite-t wall composition remain open.
source: codex-2026-07-16-S17
depends_on: [THM-727, THM-883, THM-884, THM-887-uniform-maxS-and-affine-witness-coordinate]
related: [THM-888, THM-889, THM-890, THM-892, THM-894, THM-896-level3-crossing, THM-908, THM-910, THM-911, THM-912, HYP-7021, HYP-7024, HYP-7027]
verification: 04-computation/lrc14_resonant_cross_section_cancellation_codex_S17.py -> 05-knowledge/results/lrc14_resonant_cross_section_cancellation_codex_S17.out
---

# THM-891 — exact cross-section cancellation for a resonant far peel

## 1. Statement

Let `E={0,e_1,...,e_5}` be a fixed six-offset core.  For almost every `x`, let
`M_E(x)` be the set of sectors in `Z/7Z` missed by the six points `{e x}`.  The
stationary offset forces `0 notin M_E(x)`.  Write

- `B_s = meas{x : M_E(x)={s}}`,
- `A_{s,c} = meas{x : M_E(x)={s,c}}` for `s<c`.

Append a far offset `t`, let `R_s(E,t)` be the set on which the seven-offset family
misses exactly sector `s`, and use the `THM-727` centered section function

`g_s(y)=1_[s/7,(s+1)/7)({y})-1/7`.

For a fixed positive integer `a`, define the two-scale error at the owner resonance

`Error_t(a)=sum_s integral_[R_s(E,t)] g_s(atx) dx`.

Then

`Error_t(a) = C_a(E) + O_E(1/t)`,

where

`C_a(E) = sum_s B_s sum_(c != s) J_a(c,s)
          + sum_(s<c) A_{s,c}(J_a(c,s)+J_a(s,c))`

and

`J_a(c,s)=integral_[c/7,(c+1)/7] g_s(ay) dy`.

Moreover, if `a=7q+r`, `0<=r<=6`, then

`a J_a(c,s) = (7 delta_r(c,s)-r)/49`,

where `delta_r(c,s)=1` exactly when `s` belongs to the cyclic interval
`{cr,cr+1,...,cr+r-1} mod 7`.  Consequently

> **`a C_a(E)` depends only on `a mod 7`, and it is zero when `7|a`.**

For the consecutive slow core `E={0,1,2,3,4,5}`, the six nonzero values of
`a C_a(E)` are

| `a mod 7` | exact `a C_a(E)` | decimal |
|---:|---:|---:|
| 1 | `-239/5145` | `-0.046452866861` |
| 2 | `209/20580` | `0.010155490768` |
| 3 | `-39/6860` | `-0.005685131195` |
| 4 | `47/20580` | `0.002283770651` |
| 5 | `-4/1715` | `-0.002332361516` |
| 6 | `83/6860` | `0.012099125364` |

Thus every owner resonance of this asymptotic family satisfies

`|C_a(E)| <= 239/(5145 a) <= 239/5145 < 0.097`.

This proves the signed cross-section mechanism behind the exact `THM-884` audit.  The
sectionwise absolute-value constant `0.8287` from the colliding `THM-887` is a valid
uniform envelope, but it discards precisely the signs retained above.

There is also an exact arithmetic law for any two distinct positive speeds.  Divide
`a,b` by their gcd and call the coprime pair `A,B`.  If `r=A mod 7`, then

`Pr(sector(ax)=sector(bx))
 = 1/7` when `A != B mod 7`, and

`Pr(sector(ax)=sector(bx))
 = 1/7 + r(7-r)/(7AB)` when `A = B mod 7`.

In particular every pair-collision mass is at least `1/7`, but it is not fixed.  For a
six-offset core this implies the universal low-miss estimate

> **`p_1+p_2 <= 45/49`.**

The same pair geometry also gives `p_1<=5/7`.  Consequently the exact residue-one
identity and the residue-five kernel norm imply

> **`|F_1(E)| <= 230/2401 < 0.097` and
> `|F_5(E)| <= 225/2401 < 0.097`.**

The full pair-sector law has only 21 arithmetic rays.  Exact quadratic certificates on
their endpoints additionally give

| residue | rigorous interval for `F_r=rC_r` |
|---:|---:|
| 2 | `-40/441 <= F_2 <= 230/2401` |
| 3 | `-230/2401 <= F_3 <= 19/196` |
| 4 | `-2/21 <= F_4 <= 232/2401` |

Every endpoint has absolute value below `0.097`.  Finally the raw residue-six kernel
gives `F_6<=230/2401<0.097`.  Therefore residues `1,...,5` and the positive side of
residue `6` close universally by the present theorem.  The negative side was the
remaining owner-resonant limiting problem at this checkpoint; `THM-910` subsequently
closes it at `32/343<0.097`.  The sharp `16/343` conjecture remains open for every
residue.

## 2. Conditioning identity

Let `c_t(x)=floor(7{tx})`.  Adding one point can delete at most one missed sector.
Therefore the seven-offset family misses exactly `s` in precisely two disjoint cases:

1. `M_E(x)={s}` and `c_t(x) != s`;
2. `M_E(x)={s,c}` and `c_t(x)=c`.

Ignoring the null boundary set gives the exact identity

`Error_t(a)
 = sum_s integral 1_[M={s}] 1_[c_t != s] g_s(atx) dx
 + sum_(s!=c) integral 1_[M={s,c}] 1_[c_t=c] g_s(atx) dx`.

For an unordered pair `{s,c}`, the two ordered terms are exactly the two summands
`J_a(c,s)` and `J_a(s,c)` in the theorem.

## 3. Averaging lemma and the limit

Each indicator `1_[M_E=T]` is a finite union of rational intervals, because its only
walls are `k/(7e_i)`.  If `F` is any such step function and `h` is a bounded
one-periodic step function, cell decomposition into the `t` intervals
`[j/t,(j+1)/t]` gives

`integral_0^1 F(x)h(tx) dx = (integral F)(integral h) + O_F(1/t)`.

Indeed, cells not meeting a wall of `F` contribute their constant `F` value times
`(integral h)/t` exactly; only `O_F(1)` wall cells remain, each of length `1/t`.

Apply this lemma with

- `h(y)=1_[floor(7{y}) != s] g_s(ay)`, whose mean is `sum_(c!=s)J_a(c,s)`;
- `h(y)=1_[floor(7{y}) = c] g_s(ay)`, whose mean is `J_a(c,s)`.

Substitution into the conditioning identity proves `Error_t(a)=C_a(E)+O_E(1/t)`.

## 4. Exact seven-microcell calculation

Under `z=ay`, the source interval `[c/7,(c+1)/7]` becomes `a` consecutive
microintervals `[k/7,(k+1)/7]`.  Let

`N_a(c,s)=#{k in {ca,...,(c+1)a-1}: k=s mod 7}`.

The target-sector part has measure `N_a(c,s)/(7a)`, while the centered `-1/7` part
integrates to `-1/49`.  Hence

`J_a(c,s)=N_a(c,s)/(7a)-1/49`.

If `a=7q+r`, every residue occurs `q` times and the cyclic block of `r` residues
starting at `cr` occurs once more.  Thus `N_a(c,s)=q+delta_r(c,s)` and

`aJ_a(c,s)=(q+delta_r)/7-(7q+r)/49=(7delta_r-r)/49`.

This proves the seven-residue law.

## 5. Exact consecutive-core calculation

The breakpoint sweep for `E={0,1,2,3,4,5}` gives, in sector order `s=1,...,6`,

`B=(1/35, 1/35, 5/84, 11/210, 1/28, 2/35)`.

The nonzero upper-triangular pair masses are

`A_12=1/35, A_14=13/420, A_16=1/28,
 A_24=2/105, A_25=1/21, A_26=11/210,
 A_35=3/140, A_36=5/84, A_46=1/35, A_56=1/35`.

All omitted `A_sc` are zero.  Substitution into the microcell formula produces the
six-entry table in the statement.  The verification script performs the same sweep on
the integer wall grid and checks every rational exactly.

As a direct finite-`t` referee at `a=1`, the exact errors are

| `t` | `Error_t(1)` |
|---:|---:|
| 25 | `-13/300` |
| 50 | `-689/14700` |
| 100 | `-169/3675` |
| 200 | `-449/9800` |
| 400 | `-391/8400` |

Their differences from `-239/5145` are `O(1/t)` exactly as the proof predicts.

## 6. Finite extremal evidence and the remaining crux

The same exact event sweep was run over all `15,246` primitive cores
`{0,e_1,...,e_5}` with `1<=e_1<...<e_5<=20`.  The unique largest value of
`max_r |r C_r(E)|` is

`16/343 = 0.046647230321`

at `E={0,1,2,3,4,6}`, in residue `r=1`.  This is below the propagation slack but is
not promoted to a theorem outside the finite box.  `HYP-7024` records the sharp
all-core target and the sole remaining limiting slack obligation, `F_6>-0.097`.

## 7. Exact pair-sector laws and quadratic certificates

### 7.1 Collision law

Let `I_s=1_[s/7,(s+1)/7)` and let `c_k` be the Fourier coefficient of `I_0`.
The same-sector indicator is `sum_s I_s(u)I_s(v)`.  Summing the translated Fourier
series over `s` kills every pair of frequencies except `k+l=0 mod 7`.  Integrating on
the line `(u,v)=(ax,bx)` and writing `a=gA`, `b=gB`, `(A,B)=1`, leaves

`1/7 + 7 sum_(h != 0, 7|(B-A)h) c_(Bh)c_(-Ah)`.

If `A != B mod 7`, the surviving `h` are multiples of seven and the sine factor in
`c_(Ah)` vanishes, giving exactly `1/7`.  If `A=B mod 7`, then
`Bh=Ah+7m`; the phase `(-1)^m` cancels the identical sign change of the sine, so every
remaining term is nonnegative and

`Pr(same)=1/7 + 7/(pi^2 AB) sum_(h!=0) sin^2(pi A h/7)/h^2`.

Using `sum_(h!=0) sin^2(pi h theta)/h^2=pi^2 theta(1-theta)` at
`theta=r/7` yields the claimed correction `r(7-r)/(7AB)`.  A stationary-versus-moving
pair has mass exactly `1/7` directly.

### 7.2 The 21 pair-sector rays

The preceding calculation refines without summing over equal sectors.  Let `P_ij(A,B)`
be the probability that the unordered sector pair of two positive runners is `{i,j}`;
thus the independent distribution `U` has diagonal masses `1/49` and off-diagonal
masses `2/49`.  Before symmetrizing, Fourier expansion gives

`P^lab_ij(A,B)=1/49
 + sum_(h!=0) c_i(Bh)c_j(-Ah)`.

After multiplication by `AB`, every nonconstant summand depends only on
`alpha=A mod 7` and `beta=B mod 7`: changing `A` by seven flips both its sine and
half-interval phase by the same sign.  If either residue is zero, the sine vanishes and
`P(A,B)=U`.  Otherwise, for a fixed unordered residue pair `{alpha,beta}`,

`P(A,B)=U+D_(alpha,beta)/(AB)`.

Let `(A_0,B_0)` be the coprime representative of that residue pair with least product
`m=A_0B_0`.  Then

`P(A,B)=(1-m/(AB))U + (m/(AB))P(A_0,B_0)`.

Thus every positive pair distribution lies in the convex hull of `U` and 21 endpoint
distributions.  The product-minimal representatives, with residue keys in triangular
order, are

`11:1,8; 12:1,2; 13:1,3; 14:1,4; 15:1,5; 16:1,6;
 22:2,9; 23:2,3; 24:2,11; 25:2,5; 26:2,13;
 33:3,10; 34:3,4; 35:3,5; 36:3,13;
 44:4,11; 45:4,5; 46:4,13; 55:5,12; 56:5,6; 66:6,13`.

This is the two-coordinate specialization of the concurrent `THM-890` relation-lattice
identity: the only primitive relation is `(B,-A)`, its residue class selects the ray,
and its Fourier weight is exactly proportional to `1/(AB)`.

### 7.3 Exact quadratic certificates for residues 2, 3, and 4

For a sector-count state `m=(m_0,...,m_6)` of the five moving runners, define

`n_ii=binom(m_i,2)` and `n_ij=m_i m_j` for `i<j`.

There are `462` such compositions.  Extend `K_r(M)` by zero when the inner miss set
`M` has size other than one or two, so `49F_r=E[K_r(M)]`.  Each row below specifies an
integer weight `w_ij`, a denominator `d`, and a sign.  Omitted weights are zero.

| `r,sign` | `d` | nonzero `ij:w_ij` | `L` |
|---|---:|---|---:|
| `2,+` | 2 | `05:1, 06:4, 13:4, 14:1, 15:3, 16:2, 24:1, 35:3, 46:1, 56:3` | `230/49` |
| `2,-` | 1 | `03:1, 06:1, 15:1, 23:1, 25:1, 26:1, 33:1, 44:2, 45:2, 66:1` | `40/9` |
| `3,+` | 10 | `01:2, 03:17, 05:1, 06:10, 11:3, 12:-1, 13:15, 14:-3, 15:3, 16:-2, 22:1, 23:29, 24:4, 25:10, 26:2, 33:5, 34:-4, 35:4, 36:4, 44:9, 45:-3, 46:18, 55:3, 56:-2, 66:3` | `19/4` |
| `3,-` | 1 | `02:1, 06:2, 15:1, 16:1, 25:1, 35:1, 36:1, 44:1, 45:3` | `230/49` |
| `4,+` | 5 | `04:5, 06:9, 13:3, 14:10, 15:9, 16:1, 23:3, 24:6, 26:1, 34:3, 36:8` | `232/49` |
| `4,-` | 1 | `05:1, 14:2, 15:1, 22:3, 25:-2, 26:1, 35:3, 45:1, 55:3` | `14/3` |

Writing `beta_ij=w_ij/d`, exact enumeration verifies for every composition

`sign K_r(M(m)) <= sum_(i<=j) beta_ij n_ij(m)`.

Exact evaluation on `U` and the 21 ray endpoints verifies

`10 E_P[beta] <= L`.

An actual core has ten positive runner pairs, each on one of those segments.  Averaging
the pointwise inequality and summing the ten pair expectations therefore proves

`sign 49F_r <= L`.

Dividing the six displayed `L` values by `49` gives exactly the interval table in the
statement.  This is a finite exact certificate, not a diameter scan.

### 7.4 Low-miss consequences and the remaining sign

Now let `C(x)` count colliding pairs among the six core points.  Summing the fifteen
pair laws gives `E[C]>=15/7`.  On the event `N<=2` (mass `p_1+p_2`), six points occupy
at least five sectors and hence `C<=1`; always `C<=15`.  Therefore

`15/7 <= E[C] <= (p_1+p_2)+15(1-p_1-p_2)`,

which rearranges to `p_1+p_2<=45/49`.

There is a second universal consequence.  Choose any two positive runners.  The event
counted by `p_1` forces these runners into distinct nonzero sectors.  Write `D` for
their same-sector event and `H_j` for the event that both lie in sector `j`.  Inclusion-
exclusion gives

`Pr(distinct and nonzero)=5/7-Pr(D)+2Pr(H_0)`.

Time reversal `x->1-x` gives `Pr(H_0)=Pr(H_6)`, while the two disjoint events `H_0`
and `H_6` lie inside `D`.  Hence `Pr(D)>=2Pr(H_0)`, proving `p_1<=5/7`.  Since

`F_1=-2(3p_1+p_2)/49`

and `3p_1+p_2=2p_1+(p_1+p_2)`, the two bounds yield

`|F_1| <= 2(10/7+45/49)/49 = 230/2401 < 0.097`.

For residue five, the exact integer kernel `49F_5` has absolute value at most `5` on
every singleton miss pattern and at most `4` on every pair miss pattern.  Therefore

`|F_5| <= (5p_1+4p_2)/49 <= 5(p_1+p_2)/49
         <= 225/2401 < 0.097`.

For residue six the largest singleton and pair kernels on the positive side are `6`
and `2`, so

`F_6 <= (6p_1+2p_2)/49
     = (4p_1+2(p_1+p_2))/49
     <= 230/2401 < 0.097`.

The negative side is concentrated: singleton kernels are `-1` except at sector `3`,
and the only negative pair kernels are `K_6({1,5})=K_6({2,4})=-12`.  Consequently the
sole limiting slack crux is a higher-order bound on `A_15+A_24`; the pair-distribution
relaxation alone does not close it.

This does not prove `HYP-7024`: the signed residue kernels need more than the scalar
collision count, and the sharp constant remains open.  A second tempting shortcut also
fails.  Inversion `x->1-x` does not
induce a symmetry of the six inner miss patterns because the stationary point lies on
the sector-zero boundary and remains fixed.  For the consecutive core, for example,
`B_1=1/35` but the naively paired `B_5=1/28`.

## 8. Synthesis with the concurrent routes

### 8.1 Exact finite owner-tooth reduction

For `F=E union {t}`, the finite correction is exactly

`Delta_F(w)=p_0(F union {w})-p_0(F)-p_1(F)/7`.

At the first owner tooth `w=t`, the duplicated phase cannot fill the unique missed
sector, so `Delta_F(t)=-p_1(F)/7<=0`.  It is therefore never an upper propagation
obstruction.  At `w=2t`, writing `q=floor(14{tx})` gives
`sec(tx)=floor(q/2)` and `sec(2tx)=q mod 7`; the whole finite remainder becomes the
signed fourteen-half-sector packet in `HYP-7083`.  This reduces the first bounded owner
band from an unspecified wall term to one explicit miss-label/parity incidence.

`HYP-7084` identifies the packet's class geometry exactly.  If `q=2c+eta` and
`u=2c mod 7`, doubling chooses the endpoint `u+eta` of the oriented edge
`u->u+1` on the `THM-913` parallel-class circle.  Symmetric mass and parity current
give an exact average-plus-divergence decomposition.  The `K_7` crossing Laplacian
proves `Delta_F(2t)^2<=(16/29)E_pc`, but the exact bank refutes that scalar energy as a
closure.  Restoring wall chronology gives the sharper identity

`Delta_F(2t)=C_2(E)+(1/t) sum_p(P_{m_-}({tp})-P_{m_+}({tp}))`,

where the `P_m` are a finite fourteen-cell miss-pattern palette.  The synchronized
walls `p=k/7` are residue-only and satisfy the sharp universal bound
`-23/784<=G_sync<=229/5488`.  The honest remainder is now the noncommon gcd-lattice
wall sum, the real-space counterpart of the general-cluster `THM-887` comb/CRT law.

- `THM-883` finds the owner comb in each section.  This theorem composes the sections
  before taking absolute values and identifies their exact cancellation.
- `THM-884` observes ratios `0.022--0.060` in full finite audits.  The value
  `239/5145` explains the consecutive-family resonant limit within that band.
- The other `THM-887` proves the coarse `0.8287` sectionwise envelope.  The gap to
  `0.0465` is now an explicit signed miss-pattern identity, not unexplained slack.
- `THM-888` diagonalizes the quadratic comb weights by residue.  `THM-891` is the
  linear, signed cross-section counterpart.
- `HYP-7021` retains short arithmetic relations for balanced owner combs.  Those
  relations are also necessary here: the attempted fixed collision-moment shortcut is
  false precisely on relation-rich pairs such as `(1,8)`.  Its resolution as the
  concurrent `THM-890` supplies the exact relation-spectrum language used by the pair
  rays and identifies residue-six's higher-order remainder.
- `THM-892` (the incoming invariant-mean theorem, renumbered after the `THM-891`
  collision) proves that the frame-averaged quadratic law is a divisor-lattice
  functional of the signed coincidence spectrum `N(h)`.  This is the exact arithmetic
  carrier missing from the failed relation-blind collision moment.  `THM-892` controls
  generic frames; `THM-891` computes the exceptional owner-resonant signed limit.
- `THM-892`'s S30 deviation addendum makes that complement sharper: its quadratic mean
  is a universal `pi^2/3` term plus a signed coincidence-class imbalance.  The present
  pair rays resolve every two-coordinate imbalance needed for residues `2,3,4`; the
  surviving `A_15+A_24` term is genuinely beyond pair marginals.
- `THM-894` independently names the level-three overlap tensor as the next LRC
  self-composition.  The negative residue-six kernel supplies its concrete target:
  certify the synchronized miss mass `A_15+A_24` using triple/higher relation data.
- `THM-896-level3-crossing` subsequently proves that order-three Bonferroni is coercive through its
  `m'<=11` wall and names the remaining triple-beat upper bound.  Its overlap observable
  is not identical to `A_15+A_24`, but the required relation-localized enhancement cap
  is the same analytic rung and should be developed once for both consumers.
- The resolved `HYP-7027` wall-movie probe proves that wall cycles are never silent and
  that recurrence tracks the relation spectrum.  Therefore the finite-`t` remainder
  should retain the expressive wall palette rather than quotienting cycles away.
- The constant-propagation ledger closes the coarse tail at `3.4 diam`.  A universal
  form of `HYP-7024`, plus a finite-`t` remainder bound, would attack the factor-17.6
  loss and could contract that residual band toward the natural `w>=diam` boundary.
  Pair rays now close limiting residues `1,...,5` and the positive side of `6` at the
  ledger's `0.097` slack.  The negative side of `6` and the uniform wall remainder are
  the honest remainder.

The `0.097` statement is only the current owner-resonant benchmark.  It does not replace
the ledger's all-`w` `0.2729 diam/w` envelope, certify every row margin, or make the
`O_E(1/t)` remainder uniform in the core.

## 9. Tournament and assumption challenge

The computation treats the six nonzero residues as vertices.  Its pairwise observable
is `|rC_r|-|qC_q|`; majority over primitive cores supplies the binary relation, and
lexicographic order is the tie gauge.  Both diameter-14 and diameter-20 tournaments are
transitive with Hamiltonian path `(1,2,3,4,6,5)`, score histogram one vertex at every
score `0,...,5`, singleton SCCs, no directed triangles, one Hamiltonian path, and zero
edge flips between the two boxes.

This tournament is useful telemetry, not a proof engine.  Runners, swing arcs, fixed
sections, section boundaries, wall events, residues, Fourier modes, miss patterns, and
proof obligations were all considered as vertices.  The miss patterns are the faithful
carrier: they preserve the exact two-scale coefficient but destroy finite-`t` aliasing
and wall chronology.  Residues preserve only the final error ordering.  Treating runners
or arcs as vertices would lose the signed quotient that creates the cancellation.
