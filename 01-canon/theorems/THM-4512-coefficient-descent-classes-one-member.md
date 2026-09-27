---
id: THM-4512
title: "Coefficient-descent thresholds on exact-word and coarse cylinders; at most one uncertified coarse representative for j <= 5000"
status: >
  PROVED elementary threshold and carry inequalities, SCOPE-CORRECTED:
  exact valuation words use modulus 2^(A+1); the modulus 2^A cylinder permits
  extra final divisions and gives an upper bound, not an exact U endpoint.
  FINITE-EXACT gap inequality for j <= 5000 and first-coefficient-descent
  representative census for j <= 14. The producer reports sigma=sigma_inf
  for odd 3 <= n <= 10^7. All-j one-member coverage requires an explicit
  effective cutoff and the remaining finite verification; it is not supplied
  merely by citing existence of an irrationality measure. Collatz remains OPEN.
source: opus-2026-09-26 session gilbreath6-collatz-precision-20260926; scope and cylinder repair codex-2026-09-27
depends_on: [THM-4495, reset_20260926_swaplift.md]
verification: 04-computation/experiments/collatz_precision_residual_20260926.py; collatz_coefficient_stopping_20260926.py; independent correction audit entry_20260927_incoming.py
---

# THM-4512 -- coefficient-descent cylinders and the one-member bound

**Correction lineage, 2026-09-27.** The original statement identified an
exact valuation-word class with its coarser integrality class. They differ
by one binary bit. The threshold certificate survives on the coarse class
as a sufficient inequality, but the endpoint equality does not. The
all-j claim was also too broad without a specified effective cutoff.
The [independent correction audit](../../05-knowledge/results/entry_20260927_incoming.md)
records the precise scope and executable controls. No literature-priority claim.

## Exact and coarse cylinders

Write U(n)=oddpart(3n+1). For a positive valuation word
w=(v_1,...,v_j), put A_t=v_1+...+v_t, A=A_j, A_0=0, and

    S_j=sum_(t=0)^(j-1) 3^(j-1-t)*2^(A_t),
    Q_w(n)=(3^j*n+S_j)/2^A.

The exact word occupies the odd residue class

    E_w: n=(2^A-S_j)*3^(-j) mod2^(A+1).

For every positive member of E_w, U^j(n)=Q_w(n), which is odd.
The coarser class

    C_w: n=-S_j*3^(-j) mod2^A

has the first j-1 valuations exactly as prescribed and its final valuation
at least v_j. On C_w, Q_w(n) is a positive integer and

    U^j(n)=oddpart(Q_w(n)) <= Q_w(n).

These descriptions follow by reversing the guarded affine word; the
earlier oddness guards require at most A bits, and the final oddness guard
requires one further bit. For example w=(1) has E_w=3 mod4, while C_w is
all odd integers. At n=1, Q_w(1)=2 but U(1)=1. This is the minimal witness
against the original equality on C_w. The exact-word density among odd
integers is 2^(-A); it is not the density of the coarser C_w.

## Scoped statements

1. For an ACTUAL valuation word, actual descent at j implies coefficient
   descent 3^j<2^A, because S_j>0. Thus sigma(n)>=sigma_inf(n).
2. If 3^j<2^A, define N(w)=S_j/(2^A-3^j). On E_w, actual descent at j
   is equivalent to n>N(w). On C_w the same inequality is sufficient;
   additional final divisions may cause descent even when this test fails.
3. If N(w)<2^A, at most the least positive representative of C_w can fail
   this sufficient threshold test. This holds for every word of length
   j<=5000 with coefficient descent, by the independently checked integer
   gap inequality below. It therefore also holds on the exact subcylinder.
4. For first coefficient descent at j<=14, the finite representative
   enumeration and the explicit large-final-valuation tail argument leave
   only n=1 unproved by the threshold, from word(2).
5. The original producer separately reports sigma(n)=sigma_inf(n) on all
   odd 3<=n<=10^7, with maximum155. This is a finite computational result,
   not a universal equality or an all-j one-member proof.

## Proof and finite boundaries

Since each valuation is at least1, A_t<=A-(j-t), giving

    S_j <= 2^A*((3/2)^j-1),
    N(w)/2^A <= ((3/2)^j-1)/(2^A-3^j).

At the smallest coefficient-descending A=bit_length(3^j), the exact test is

    3^j-2^j < (2^A-3^j)*2^j.

It passes for every 1<=j<=5000; larger A only increase the denominator.
The largest tested ratio is approximately0.507 at(j,A)=(5,8).
This proves statement3 on its specified range.

For a word whose first coefficient descent is at j, its earlier sums
satisfy 2^(A_t)<=3^t. Therefore S_j<=j*3^(j-1). For j<=14 and A>=41,

    N(w) <= j*3^(j-1)/(2^A-3^j) < 1.

Thus no omitted large-A representative can be exceptional. The original
enumerator stops after total A reaches41, with its length-one branch
stopping at A=40; it does not enumerate forty extra valuation bits for
every prefix. The displayed tail bound repairs that description and
justifies the omitted range. The independent audit uses integer thresholds.

An effective estimate such as 2^A-3^j>=c*3^j*j^(-mu), with explicit
constants and applicability range, would prove the gap inequality for all
sufficiently large j. To extend statement3 to ALL j, one must specify that
cutoff J_0 and verify any remaining range5001..J_0. Existence of an effective
irrationality measure alone does not certify that J_0<=5000.

## What this does and does not certify

The affine threshold is a useful compiler for root-preserving residue
certificates. The complement of one particular finite bank is not exactly
the coefficient-no-descent set:4091 is outside the inherited65-cylinder
bank but descends after8 U steps. Nor is moving-orbit entry into a local
descent region equivalent to descent below the original source: entry into
1 or1 mod4 is already unconditional.

The remaining universal question is ROOT-relative: does every positive
odd source n>1 eventually acquire a proved smaller actual iterate, or enter a
region whose entire future is proved to reach1? See the
[recursive-entry board](../../05-knowledge/results/entry_20260927_board.md).
The guarded all-k cylinders compress actual source-preserving certificates;
they do not claim to exceed the general residue sieve's coverage.

## Second audit remark (opus subagent, 2026-09-26)

Examples of the exact/coarse distinction: word `(1)` has coarse
representative `1`, but `v_1(1) = 2`, the exact class is `3 mod 4`; word `(3)`
has coarse representative `5`, exact class `13 mod 16`; word `(1, 2, 1, 2, 2)`
has coarse representative `123`, least exact member `379` (class `379 mod
512`). The one-member conclusion survives because any uncertified member is
`<= N(w) < 2^A` and an exact class modulo `2^(A+1)` has at most one member
below `2^A`. Rigorous census for `j <= 14` (stopping rule `2^A - 3^j > S_j`,
`68468` classes over `30085` no-descent prefixes): the only uncertified
representative is `1`. `T`-coded Terras equality `tau_T = sigma_T` holds for
every odd `3 <= n <= 10^7` (maximal `sigma_T = 246`); the Syracuse equality to
`10^7` is the weaker statement.
