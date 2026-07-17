# HYP-7085 — Farey divisor sheets and the Zarankiewicz guardrail

**Status:** EXACT DIVISOR-SHEET / COLORED-FIBER REDUCTION PROVED; ORDINARY
ZARANKIEWICZ AND CLASS-PARITY CLOSURES REFUTED; EXACT ADJACENT-PIN DENSITY
BOUNDS PROVED; THE DIVISOR-SUMMED FAR-ADDRESS BOUND REMAINS OPEN; THE
PROVISIONAL COEFFICIENT CAP IS REFUTED AND REPLACED BY AN EXACT LEAST-SPEED
ENVELOPE (codex-2026-07-16-S20).

Let `A={a_1,...,a_5}` be the positive slow core and `t>max(A)`.  This note
refines the noncommon owner-wall remainder in `HYP-7084`.

## 1. Farey/divisor sheets

At a slow wall, reduce

```text
7x=h/d,  gcd(h,d)=1.
```

A speed `a` crosses there exactly when `ah/d` is integral, hence exactly when
`d|a`.  Therefore the simultaneous owner set is

```text
S_d(A)={a in A:d|a}.                                         (1)
```

The synchronized deck is `d=1`; the noncommon remainder is the disjoint sum
over `d>=2`.  This is the precise Farey-layer version of the `THM-887`
comb/CRT law and the torsion-class dictionary in the parallel-class-circle
reflection.

Now fix a proper owner set `S`, put `g=gcd(S)`, and use columns `k mod 7g`.
The columns with exact owner set `S` are

```text
K_A(S)={k mod 7g : g|ak exactly for a in S}.
```

For `O=A\S`, inclusion-exclusion gives the proved count

```text
N_A(S)=7 sum_(T subset O) (-1)^|T| gcd(g,T),                 (2)
```

with `gcd(g,empty)=g`.  Indeed, imposing `g|bk` for every `b in T` leaves
exactly `7 gcd(g,T)` columns.  In particular every exact column count is a
multiple of seven.

Let `m_k^-` and `m_k^+` be the missed-sector sets immediately before and after
column `k`.  The exact owner-set contribution is the signed cyclic word

```text
G_S(A,t)=1/t sum_(k in K_A(S))
  [P_(m_k^-)({tk/(7g)})-P_(m_k^+)({tk/(7g)})].               (3)
```

Equations (1)--(3) were fraction-refereed at every wall of all `6,900` rows in
the `HYP-7084` bank.

## 2. Exact class-parity inclusion-exclusion

Define

```text
eta_t(k)=floor(2{tk/g}) in {0,1}.
```

For `h_T=gcd(g,T)`, `c_T=gcd(h_T,t)`, and `delta_T=h_T/c_T`, the exact parity
counts are

```text
N_eta(A,S;t)=7 sum_(T subset O) (-1)^|T| E_eta(h_T,t),        (4)
E_0(h,t)=c ceil(delta/2),  E_1(h,t)=c floor(delta/2).
```

To prove (4), the conditions indexed by `T` force `k=(g/h_T)r`; the fractional
parts `{tr/h_T}` run through the `delta_T`-grid with multiplicity `c_T`, seven
times.  Thus

```text
N_0-N_1=7 sum_T (-1)^|T| c_T 1_[delta_T odd].                (5)
```

This is the literal parity inclusion-exclusion suggested by the cyclic
Zarankiewicz book.  It is exact, but it does not control the signed word.

## 3. The affine colored `K_{8,7}` fiber

On a fixed divisor sheet, write

```text
h=a+jd,  a in U_d,  j in Z_7.
```

For every row speed `v`, its post-wall sector color on this seven-column fiber
is

```text
A_v^+(j)=floor(va/d)+vj mod 7,
A_v^-(j)=A_v^+(j)-1_[d|v].                                  (6)
```

Adding the rows `0,A,t,2t` makes an affine edge-colored `K_{8,7}`.  The two
distinguished rows obey

```text
A_(2t)(j)=2A_t(j)+eta_(d,a) mod 7,
eta_(d,a)=floor(2{ta/d}).                                   (7)
```

The waveform is not a crossing count.  It is a pointed rainbow-completion
score: among the first seven rows `F={0} union A union {t}`, if exactly one
sector `s` is absent, score `1_[A_(2t)=s]-1/7`; otherwise score zero.  The
primitive `P_m` in (3) integrates this colored score.  Thus the natural
Zarankiewicz object is a **colored** `K_{8,7}` with ordered columns and
distinguished rows, not its ordinary uncolored drawing.

For a unit `a in U_d`, call its seven-column sum the **pinned cyclic rainbow
discrepancy** `D_b(d,a;t)`.  It is equivalently a colored
Dedekind--Rademacher sum: it retains the root pin `b`, cyclic color order,
pre/post miss flag, and far source address.  The exact owner word is

```text
tG_S=sum_(d|g, d>=2, S_d(A)=S) sum_(a in U_d) D_0(d,a;t),
D_0(d,d-a;t)=D_6(d,a;t).                                   (7a)
```

Thus opposite unit fibers are precisely the adjacent-pin symmetrization of
the cut-open parallel-class circle.  Ordinary crossing number is the quotient
obtained after discarding every field in `D_b` that carries its sign.

## 4. What the Zarankiewicz values do and do not say

The ordinary crossing-number status in the incoming notes needed correction:

```text
cr(K_7,7)=81,  cr(K_7,8)=108,  cr(K_8,8)=144.                (8)
```

Woodall proved the `K_{7,7}` value in 1993.  Deleting one vertex in the
eight-vertex part gives

```text
(m-2) cr(K_(m,n)) >= m cr(K_(m-1,n)),
```

so `81` propagates to the lower bounds `108` and `144`; Zarankiewicz drawings
give equality.  These are known cases, not open instances.  The general
conjecture remains open.

For the unlabelled owner-column incidence graph, every owner meets every exact
column, giving `K_(|S|,N_A(S))`.  Proper owner sets have `|S|<=4`, so their
ordinary crossing numbers equal the Zarankiewicz formula unconditionally.  In
the exact bank:

- `34,791` proper owner words have only `Z in {0,9,18,42}`;
- `30,702/34,791 = 0.882470...` have `Z=0`, including both individual signed
  extrema `-53/1372` and `131/4116`;
- planar-word row sums already range from `-1441/41160` to `172/5145`;
- the positive-`Z` part ranges only from `-97/4116` to `15/686`;
- correlations of aggregate exact-`Z`, full-grid-`Z`, class cut, and crossing
  energy with `|G_noncommon|` are respectively about
  `-0.0116,-0.0113,-0.0459,-0.0449`.

The failures are exact, not just correlations:

1. The same owner set `(2,4,6,8)`, quotient arc `(1,2,3,4)`, parity split,
   `K_{4,7}` with `Z=18`, class cut, and crossing energy gives
   `-43/3430` at `A=(2,3,4,6,8),t=10` but `15/686` at
   `A=(1,2,4,6,8),t=12`.
2. With the core and owner set fixed,
   `A=(1,3,5,6,7), S=(7)`, the contribution is `-53/1372` at `t=8` and
   `11/882` at `t=9`; both have the same exact parity census `(21,21)`.
3. Stronger: `A=(1,2,3,4,5), S=(2,4)` has the identical full parity word
   `eta=(0,0,0,0,0,0,0)` at `t=6` and `t=8`, yet its contribution flips from
   `23/4116` to `-37/5488`.

Therefore ordinary crossing number, the entire binary parity word, cyclic cut,
and crossing energy all forget proof-bearing data: the full source address
`floor(tk/g) mod 7`, the pre/post miss word, the outside-divisibility mask, and
the cyclic column order.

## 5. Cut-open-circle reflection

Let `R(s)=-s-1 mod 7`.  Directly from the fourteen half-sectors,

```text
H_(Rm)(13-q)=H_m(q),
P_(Rm)(1-z)=-P_m(z).                                        (9)
```

If `W_S(k;b)` denotes a column contribution with the stationary pin at sector
`b`, then negating the column swaps pre/post states under `R`, so

```text
W_S(-k;b)=W_S(k;R(b)).                                      (10)
```

The actual pin is `b=0`, while `R(0)=6`.  Paired columns therefore give
**adjacent-pin symmetrization**, not cancellation.  This is the cut-open face
of the parallel-class circle and explains structurally why the adjacent class
has crossing profile `xi(1)=0` while the signed owner current remains nonzero.

There is also an exact owner-order gauge.  Let `c in (Z_7)^5` be the five slow
row colors and put `Phi_z(c)=P_(m(c))(z)`, with the stationary pin retained.
A wall with owner set `S` moves `c` to `c+1_S`.  For every ordering
`pi=(s_1,...,s_r)` of `S`,

```text
Phi_z(c)-Phi_z(c+1_S)
 =sum_(j=1)^r [Phi_z(c+1_{s_1,...,s_(j-1)})
                -Phi_z(c+1_{s_1,...,s_j})].                 (10a)
```

Thus a simultaneous owner wall is a path-independent exact gradient on the
five-row color cube; the ordering is a gauge that may be chosen to expose
runnerwise or reflection pairing.  Exact enumeration of the one-row color
edges on the fourteen-cell palette gives

```text
max single-edge |Delta Phi|=51/686.                         (10b)
```

This is the correct Hamiltonian-path role suggested by Tournament Analysis:
the path orders proof obligations inside one owner fiber, rather than ranking
runners globally.  It can reassemble multi-owner sheets into runner currents,
but it cannot by itself close singleton sheets, where there is no ordering
freedom.

The adjacent-pin gauge also gives rigorous local bounds.  Exact enumeration of
all five-row color states, owner masks, and the fourteen phase endpoints gives
the paired ranges

```text
|S|     min[W(k;0)+W(k;6)]     max[W(k;0)+W(k;6)]
 1              -1/7                    11/98
 2              -1/7                     1/7
 3             -13/98                   45/343
 4             -13/98                   15/98.              (10c)
```

The exact column set is invariant under `k -> -k`.  Its only possible proper
fixed point is the half-turn; there reflection makes its score exactly half a
paired score.  If `N=N_A(S)`, (10c) therefore implies in particular

```text
-N/(14t) <= G_S <= 11N/(196t)        for |S|=1,
|G_S| <= N/(14t)                     for |S|=2.              (10d)
```

Since `N<=7g`, the singleton absolute value is `<1/2` and its positive side is
`<11/28`.  Distinct pair owners force `max(S)>=2g`, hence the pair bound is
`<1/4`.  The higher-cardinality rows of (10c) similarly give absolute caps
`<13/84` and `<15/112`.  These are the first universal sheet estimates that
use the parallel-class reflection, but they remain too coarse when summed.

The correct no-double-counting charge is explicit.  Let `M_r(A)` be the number
of slow wall positions having exactly `r` owners and put

```text
I_k=7 sum_(T subset A, |T|=k) gcd(T).
```

Every `k`-fold wall intersection is counted once for each `k`-subset of its
owners, so binomial inversion gives

```text
M_r=sum_(k=r)^5 (-1)^(k-r) binom(k,r) I_k.                  (10e)
```

Consequently

```text
sum_(|S|=1) G_S <= 11M_1/(196t),
sum_(|S|=2) G_S <= M_2/(14t).                              (10f)
```

This is the honest divisor-simultaneity aggregate.  Its remaining loss is now
precise: it discards alignment of the far addresses and signs across sheets.

## 6. The sharpened remaining crux

The full bounded-bank frequency scan evaluates every residue frequency on all
`1,259` distinct proper owner words (`49,483` exact evaluations).  It gives

```text
-1019/2058 <= t G_S <= 89/196,
```

with both extrema on singleton, hence planar, sheets.  The lower value is
`-1/2+5/1029`.  This is a diameter-10 box artifact, not evidence for a
universal cap.  The first larger diameter already gives the exact refutation

```text
A=(7,8,9,10,11), S=(11), t=44:
t G_S=-216/343=-(6/7)^3,   G_S=-54/3773.                   (11)
```

There are 70 exact columns, 38 with nonzero local score.  Thus the periodic
coefficient can accumulate with the sheet denominator; it is the wrong
normalization.

This is not an isolated or merely favorable negative miss.  For every `g>=5`,
take the consecutive-tail singleton family

```text
A_g=(g-4,g-3,g-2,g-1,g),  S=(g),
C_g^(m)=mg G_S(A_g,mg),  m in {2,4}.
```

The far phases are locked to `mk/7`, the parallel-class circle coordinate.  For
an exact singleton column, the four outside rows occupy
`k-ceil(ck/g) mod 7`, `c=1,2,3,4`, while the owner moves from `k-1` to `k`.
The local word therefore depends only on `k mod 7` and the four ceilings.
Partitioning `0<=k/g<7` at the rational breakpoints of those ceilings gives a
degree-one quasipolynomial in `g` with period dividing
`lcm(1,2,3,4)*7=84`.  Exact evaluation of its 84 residue classes gives

```text
C_(g+84)^(2)=C_g^(2)+62/49,
C_(g+84)^(4)=C_g^(4)-149/49.                               (12)
```

Consequently `tG_S` is unbounded in both signs even on planar `K_(1,N)` owner
graphs.  The positive half-cap already fails at `g=33,t=66`, where

```text
tG_S=53/98,  G_S=53/6468.                                  (13)
```

The normalized currents instead have finite phase-locked limits

```text
G_S(A_g,2g) -> 31/4116,
G_S(A_g,4g) -> -149/16464.                                 (14)
```

Pair sheets have the same phenomenon.  On `g=11 mod 84`, set

```text
A=(g,2g-3,2g-2,2g-1,2g), S=(g,2g), t=3g.
```

Then the owner incidence is again planar and the exact phase ray satisfies

```text
C_(g+84)=C_g+797/343,  G_S -> 797/86436.                   (14a)
```

Thus neither singleton nor pair coefficients admit even a one-sided constant
cap; the normalized density bounds (10d)--(10f) are essential.

The number `84=12*7` is exactly the merge of the four-row Farey breakpoint
grid with the seven-vertex parallel-class circle.  Ordinary Zarankiewicz value
is identically zero on this family and cannot see the drift.

There is, however, an exact finite replacement.  Put `D=max(A)`, `p=7g`, and
write the numerator of (3) as the periodic coefficient

```text
C_(A,S)(r)=sum_(k in K_A(S))
  [P_(m_k^-)({rk/p})-P_(m_k^+)({rk/p})],   r mod p.
```

For `0<=r<p`, the least admissible far speed in that residue is

```text
tau_D(r)=r+p(floor((D-r)/p)+1).                             (15)
```

Every `t>D` with `t=r mod p` has the same numerator, so

```text
max_(t>D, t=r mod p) |G_S(A,t)|=|C_(A,S)(r)|/tau_D(r).       (16)
```

Consequently the infinite far-speed problem for a fixed word is exactly the
finite `p`-entry **least-speed envelope**.  This is the normalization that the
diameter-10 coefficient scan omitted.

The fourteen-cell palette also gives a universal, though presently too weak,
guardrail.  Since every primitive difference is piecewise linear on the common
half-sector grid, exact endpoint enumeration gives

```text
max_(m,n,z) |P_m(z)-P_n(z)|=135/1372.                       (17)
```

For a realized proper owner set, choosing one outside speed shows
`N_A(S)<=7(g-1)`.  Hence

```text
|G_S(A,t)| <= 135 N_A(S)/(1372 t)
             <= 945(g-1)/(1372 t) < 945/1372.               (18)
```

The bound (18) is not remotely summable enough to close the global remainder,
but it cleanly separates local palette size from divisor multiplicity.  The
next proof-bearing target is a **divisor-summed adjacent-pin discrepancy bound**
for the least-speed envelopes (16), beginning with singleton and pair sheets.
Those graphs are always planar, yet they carry the largest coefficients; the
outside mask, full far address, and simultaneous charge across divisors must
stay in the state.

## 7. The slice-Parseval splice

The concurrent `HYP-7103` result supplies the correct quadratic sidecar for
this target.  In the `THM-888` owner-frequency notation it proves, for every
owner `e` and residue class `r`,

```text
sum_(m=r mod e) |S_e(m)|^2=(P/e) C_hat_e(r),                 (19)
```

where `C_hat_e(r)` is the `r`-twisted same-class endpoint-coincidence sum.
This is exact coset Parseval, not an envelope: it removes the former
quadratic-in-owner-count floor and improved the tested two-owner tail by a
factor `636`.

The compatible linear splice to prove is a regrouping of the colored divisor
sheet into owner-frequency slices of the form

```text
T_e(r)=sum_(m=r mod e) K_hat_P(m) S_e(m).
```

Once that Fourier identity is matched exactly to (3), (19) gives immediately

```text
|T_e(r)|^2
 <= [sum_(m=r mod e)|K_hat_P(m)|^2] (P/e) C_hat_e(r).        (20)
```

Equation (20) would be the sharp norm input missing from (18), but it is not itself
the desired uniform signed estimate.  It controls one fixed far-address class
in `L2`; summing its square roots independently would again discard the
adjacent-pin signs and divisor simultaneity.  The remaining crux is now
precisely to prove the colored regrouping and combine (20) with (10e), retaining
the common twist `r=t` and the pre/post miss labels across all divisor charges.

This also challenges a tempting identification: the Parseval owner `e` is a
frequency-slice owner and should not be silently identified with the reduced
Farey denominator `d`.  The regrouping preserves the fixed-address LRC
functional and the twisted coincidence mass, while ordinary Zarankiewicz
crossing numbers discard both.  Thus `HYP-7103` strengthens the divisor-sheet
route without reviving the refuted crossing-energy quotient.

Tournament Analysis uses owner cardinalities `1,2,3,4` as vertices.  Positive
signed risk orders them `(1,2,4,3)`, whereas Zarankiewicz risk orders them
`(3,4,1,2)`; the switch flips five of six edges.  Both are transitive with
singleton SCCs and unique tie Hamiltonian paths, confirming that the scalar
ordering is telemetry, not the carrier.

Verification:
`04-computation/lrc14_farey_divisor_sheet_zarankiewicz_codex_S20.py` and the
matching result file.

External status source: D. R. Woodall, *J. Graph Theory* 17 (1993), 657--671,
doi:`10.1002/jgt.3190170602`.

-> `HYP-7084`, `THM-887`, `THM-913`, cyclic Zarankiewicz `THM-922`,
`HYP-7103`, and the parallel-class-circle reflection.
