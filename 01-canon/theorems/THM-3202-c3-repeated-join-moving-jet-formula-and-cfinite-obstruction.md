---
id: THM-3202
title: "C3 repeated-join moving-jet formula and C-finite obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  nonempty tournament A, the fixed-depth path-cover coordinates of
  C3[A^{join r},A^{join r},A^{join r}] have an exact quotient-walk formula
  obtained by applying the C3 kernel to the full finite-difference jet of
  Q_A(t)^r.  The factor jet costs quadratic arithmetic in the block order and
  the fixed-depth C3 transfer is a bounded-width diagonal band.  Nevertheless
  every fixed positive output coordinate has a universal factorial-cube lower
  bound, so its sequence in r has no fixed finite polynomial-exponential
  form, no eventual constant-coefficient recurrence, and no rational ordinary
  generating function.  Scalar Hamiltonian data is insufficient; the
  complete seed path-cover polynomial and quotient-walk series are the
  established restoring sidecars.  No arbitrary-tournament complexity
  collapse, variable-coefficient recurrence obstruction, or
  information-theoretic minimality claim is made.
audit: >
  The primary moving-jet engine checks A=K1 through r=8, A=C3 through r=5,
  fixed output depths 1,2,3, quotient-kernel versus Hamiltonian-diagonal
  agreement, singleton endpoints, factorial-cube floors, integrality, the
  Stirling boundary, and the same-H mask-40/mask-76 hostile.  An independent
  vertex engine imports neither the quotient kernel nor finite differences;
  Held--Karp endpoint counts and canonical set partitions verify four
  low-profile controls, K1 lifts through 12 vertices, and C3 repeated-join
  lifts through 18 vertices.  Normal, optimized, and stored transcripts match
  byte for byte.
source: root/frontier-synthesis/tournament-moving-jet/2026-08-02
depends_on:
  - THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence
  - THM-3121-path-cover-walk-content-substitution-kernel
  - THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform
related:
  - THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant
script: 04-computation/tournament_c3_repeated_join_moving_jet_thm3202.py
output: 05-knowledge/results/tournament_c3_repeated_join_moving_jet_thm3202.out
script_sha256: e17421e4093d8f2d4ea27db604191e832a555817bbf0171971aac3e4fda7d1db
output_sha256: 05447a87519f536f685866020cbc66624608926839f55a8d4bb4caf8febd3c39
independent_script: 04-computation/tournament_c3_repeated_join_moving_jet_independent_audit_thm3202.py
independent_output: 05-knowledge/results/tournament_c3_repeated_join_moving_jet_independent_audit_thm3202.out
independent_script_sha256: 7bda0449f9a3e0ddc05649e1bd10ade62dbb58f976889e4519fb43e4bf2a6ac9
independent_output_sha256: 8634e87568dab6e85e6beb4a8e5dde8e1eeb6dde95924b606d00136f7823531f
hash_basis: LF-normalized bytes
---

# THM-3202 -- `C3` repeated-join moving-jet formula and C-finite obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3181 makes every fixed path depth of a repeated order join a finite
exponential sum in the join level.  Cyclic substitution changes that behavior
at a precise boundary: its quotient walks revisit a factor and sample the
factor jet all the way to its moving endpoint.  The resulting coordinates
remain exactly and efficiently computable in arithmetic operations, but a
universal factorial-cube contribution rules out every fixed C-finite form.

## 1. Setup and exact factor jet

Let `A` be any nonempty tournament of order `a>=1`.  For every integer
`r>=1`, put

\[
 B_r=A^{\triangleright r},\qquad N=ar,\qquad
 U_r=C_3[B_r,B_r,B_r].                                    \tag{1}
\]

For a tournament `S`, let

\[
 p_S(c)=\#\{\text{spanning covers by }c\text{ unordered directed paths}\},
 \qquad F_S(c)=c!p_S(c),                                  \tag{2}
\]

and let

\[
 Q_S(t)=\sum_c p_S(c)(t)_{\underline c}.                  \tag{3}
\]

Write `F_r(c)=F_(B_r)(c)`, extended by zero outside `1<=c<=N`.  THM-3181
gives `Q_(B_r)(t)=Q_A(t)^r`.  Falling-factorial inversion therefore gives,
for every `0<=c<=N`,

\[
 \boxed{F_r(c)=\Delta_j^c(Q_A(j)^r)|_{j=0}
 =\sum_{j=0}^c(-1)^{c-j}\binom cj Q_A(j)^r.}              \tag{4}
\]

For each fixed factor depth `c`, (4) is a finite exponential sum in `r`.
The cyclic operation below, however, requires every depth through the moving
endpoint `N=ar`.

## 2. Exact fixed-output moving-jet formula

Orient `C_3` as `1->2->3->1`.  Its THM-3121 quotient-walk series is

\[
 W(x,y,z)=\frac{x+y+z+xy+xz+yz+3xyz}{1-xyz}.              \tag{5}
\]

Put

\[
 w_d(c_1,c_2,c_3)=[x^{c_1}y^{c_2}z^{c_3}]W(x,y,z)^d.      \tag{6}
\]

The complete fixed-depth substitution law is, for every `r>=1` and every
`1<=d<=3N`,

\[
 \boxed{F_{U_r}(d)=
 \sum_{1\le c_1,c_2,c_3\le N}
 w_d(c_1,c_2,c_3)F_r(c_1)F_r(c_2)F_r(c_3).}               \tag{7}
\]

In particular, the balanced diagonal and its two nearest shells give

\[
 \boxed{H(U_r)=3\sum_{c=1}^{N}\left(
 F_r(c)^3+F_r(c+1)F_r(c)^2+F_r(c+1)^2F_r(c)\right).}      \tag{8}
\]

Equations (4), (7), and (8) are explicit exact formulas.  If (4) is expanded
inside (7), its exponential bases have the form

\[
 Q_A(j_1)Q_A(j_2)Q_A(j_3),\qquad 0\le j_i\le N.           \tag{9}
\]

Their index range and spectrum grow with `r`; this is a moving finite sum, not
a fixed finite polynomial-exponential sequence.

## 3. Complexity and the diagonal band

The complete row `(F_r(1),...,F_r(N))` is obtained from
`Q_A(0)^r,...,Q_A(N)^r` by an in-place forward-difference table in

\[
 O(N^2)\text{ arithmetic and }O(N)\text{ storage}.        \tag{10}
\]

For fixed `d`, the transfer in (7) takes only `O_d(N)` additional arithmetic.
Indeed, with

\[
 P=x+y+z+xy+xz+yz+3xyz,
\]

the expansion `W^d=P^d(1-xyz)^(-d)` gives

\[
 w_d(c_1,c_2,c_3)=\sum_{q\ge0}\binom{d+q-1}{q}
 [x^{c_1-q}y^{c_2-q}z^{c_3-q}]P^d.                       \tag{11}
\]

Every exponent of `P^d` lies in `[0,d]^3`, so the nonzero kernel cells form a
bounded-width diagonal band.  For `d=1`, (8) is an `O(N)` transfer after the
factor jet.  Thus any fixed output depth has an exact `O(N^2)` arithmetic,
`O(N)` storage evaluator.  These are arithmetic-operation bounds; no
unit-cost bit-complexity or universal algorithmic lower bound is claimed.

## 4. Universal factorial-cube lower bound

Every `N`-vertex tournament has exactly one cover by `N` singleton paths, so

\[
 F_r(N)=N!.                                                \tag{12}
\]

Fix an integer `d>=1` and suppose `N>=d`.  In `P^d`, choosing `3xyz` from
every factor contributes `3^d(xyz)^d`; in `(1-xyz)^(-d)`, the coefficient of
`(xyz)^(N-d)` is `binom(N-1,d-1)`.  Hence

\[
 w_d(N,N,N)\ge3^d\binom{N-1}{d-1}.                       \tag{13}
\]

Every term in (7) is nonnegative.  Retaining only the cell
`(c_1,c_2,c_3)=(N,N,N)` proves

\[
 \boxed{F_{U_r}(d)\ge
 3^d\binom{N-1}{d-1}(N!)^3,}                              \tag{14}
\]

\[
 \boxed{p_{U_r}(d)\ge
 \frac{3^d}{d!}\binom{ar-1}{d-1}((ar)!)^3.}              \tag{15}
\]

In particular, for every `r>=1`,

\[
 \boxed{H(U_r)\ge3((ar)!)^3.}                            \tag{16}
\]

At every fixed positive integer `m`, positivity in the falling-factorial
basis gives

\[
 Q_{U_r}(m)=\sum_{d=1}^m\binom mdF_{U_r}(d)
            \ge mH(U_r),
\]

and therefore

\[
 \boxed{Q_{U_r}(m)\ge3m((ar)!)^3.}                       \tag{17}
\]

## 5. C-finite, fixed-exponential, and rational-OGF obstruction

Every eventually constant-coefficient recurrent complex sequence is a finite
sum of terms `P_i(r)lambda_i^r` and is bounded in magnitude by
`C r^M R^r`.  The factorial-cube bounds (14)--(17) exceed every such bound.

Consequently, for every fixed nonempty seed `A`, every fixed integer `d>=1`,
and every fixed integer `m>=1`, each sequence

\[
 \{p_{U_r}(d)\}_{r\ge\lceil d/a\rceil},\quad
 \{F_{U_r}(d)\}_{r\ge\lceil d/a\rceil},\quad
 \{H(U_r)\}_{r\ge1},\quad
 \{Q_{U_r}(m)\}_{r\ge1}                                 \tag{18}
\]

has all three properties:

1. it is not a fixed finite polynomial-exponential sum;
2. it satisfies no nonzero constant-coefficient linear recurrence valid for
   all sufficiently large `r`; and
3. its ordinary generating function in `r` is not rational.

The third conclusion follows because coefficients of a rational formal
ordinary generating function satisfy an eventual constant-coefficient
recurrence.  The word **fixed** applies to the seed, output coordinate or
colour, proposed bases, and recurrence coefficients.  Variable-coefficient,
P-recursive, holonomic, differential-algebraic, and factorial-normalized
descriptions are not excluded.

This positive growth obstruction is stronger than a failed recurrence fit to
a finite prefix.

## 6. Lost data and restoring sidecar

The cyclic quotient can revisit a factor block.  Passing from the complete
factor profile to the scalar `H(A)` forgets the number of maximal block runs
and the allocation of vertices among those runs—the component index `c` in
(4)--(8).

The loss is explicit.  The THM-3134 labelled order-five tournaments with edge
masks `40` and `76` both have Hamiltonian count `15`, but factorial profiles

\[
 (15,78,198,240,120),\qquad(15,90,210,240,120),            \tag{19}
\]

and equal `C_3` lifts with

\[
 H(C_3[T_{40},T_{40},T_{40}])=178036299,
\]

\[
 H(C_3[T_{76},T_{76},T_{76}])=193215375.                  \tag{20}
\]

Thus scalar `H` is refuted as a substitution sidecar.  The **minimal
established functorial sidecar** is the complete seed path-cover profile,
equivalently `Q_A`, together with the quotient walk series `W`.  At level `r`
it generates the moving endpoint jet (4).  Information-theoretic minimality
among arbitrary nonlinear compressions is not asserted.

The source is the fixed seed profile, the target is a fixed path-cover
coordinate of the cyclic lift, and the map is order-join multiplication,
finite-difference inversion, then maximal-block-run threading through `W`.
The preserved predicate is the exact output coordinate; scalarization
destroys run/component depth; `Q_A`, the moving jet, and quotient walk content
restore it.

## 7. Seed controls and scope boundary

For `A=K_1`, `Q_A(j)=j`, `B_r` is transitive, and

\[
 F_r(c)=c!S(r,c).
\]

The first Hamiltonian values are

\[
 3,45,2721,421425,132518793,\ldots.                       \tag{21}
\]

For `A=C_3`, `Q_A(j)=j^3+2j`, and the first values are

\[
 3159,82069875945,195780715669875147291,\ldots.           \tag{22}
\]

This theorem supplies an exact polynomial-arithmetic evaluator in the join
level and a definitive obstruction to C-finiteness.  It does not compute an
arbitrary seed profile, make arbitrary tournament path counting
polynomial-time, produce a bounded endpoint jet for cyclic substitution,
exclude variable-coefficient descriptions, or establish
information-theoretic sidecar minimality.

The live questions are whether a natural factorial normalization admits a
useful variable-coefficient recurrence, and whether (8) can be evaluated
subquadratically without explicitly forming the full finite-difference row.

## 8. Exact verification

Run

```bash
python3 04-computation/tournament_c3_repeated_join_moving_jet_thm3202.py
python3 -O 04-computation/tournament_c3_repeated_join_moving_jet_thm3202.py
python3 04-computation/tournament_c3_repeated_join_moving_jet_independent_audit_thm3202.py
python3 -O 04-computation/tournament_c3_repeated_join_moving_jet_independent_audit_thm3202.py
```

The primary moving-jet engine checks:

- `A=K_1` through `r=8` and `A=C_3` through `r=5`;
- every displayed output coordinate of depths `1,2,3`;
- the general quotient-walk coefficient engine against the Hamiltonian
  diagonal formula;
- singleton endpoints, factorial-cube lower bounds, and ordered/unordered
  integrality;
- the Stirling boundary for the transitive seed; and
- the same-`H` hostile (19)--(20).

The independent engine imports neither the transfer kernel nor the
finite-difference formula.  A flat Held--Karp endpoint DP and a canonical
least-vertex set-partition recurrence verify four complete low-profile
controls, `A=K_1` through 12 vertices, and `A=C_3` through 18 vertices.  Both
normal and optimized transcripts byte-match their stored outputs, and all
four artifact hashes match the YAML header.

**End of proof.**
