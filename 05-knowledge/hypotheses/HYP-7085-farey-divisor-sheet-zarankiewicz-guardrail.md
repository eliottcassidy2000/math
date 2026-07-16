# HYP-7085 — Farey divisor sheets and the Zarankiewicz guardrail

**Status:** EXACT DIVISOR-SHEET / COLORED-FIBER REDUCTION PROVED; ORDINARY
ZARANKIEWICZ AND CLASS-PARITY CLOSURES REFUTED; THE ADJACENT-PIN SIGNED-WORD
BOUND REMAINS OPEN (codex-2026-07-16-S20).

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

## 6. The sharpened remaining crux

The full bounded-bank frequency scan evaluates every residue frequency on all
`1,259` distinct proper owner words (`49,483` exact evaluations).  It gives

```text
-1019/2058 <= t G_S <= 89/196,
```

with both extrema on singleton, hence planar, sheets.  The lower value is
`-1/2+5/1029`.  This supports, but does not prove, the local conjecture

```text
|t G_S(A,t)| <= 1/2.                                       (11)
```

Even (11) alone does not sum to the required global propagation bound without
using divisor-sheet simultaneity.  The next proof-bearing target is an
**adjacent-pin discrepancy bound** for (3), beginning with singleton and pair
sheets.  Those graphs are always planar, yet they carry the largest signed
coefficients; the outside mask and full far address must stay in the state.

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

-> `HYP-7084`, `THM-887`, `THM-913`, cyclic Zarankiewicz `THM-922`, and the
parallel-class-circle reflection.
