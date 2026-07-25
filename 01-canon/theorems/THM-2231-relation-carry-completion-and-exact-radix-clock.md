---
id: THM-2231
title: "Relation-carry completion, primitive depth pumping, and the exact radix clock"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every radix q>=2 and integer relation m,
  infinite digit streams admitting integral carries are exactly the kernel
  of m on the q-adic completion. The carry alphabet is uniformly bounded;
  for mixed-sign m it has at most ||m||_1-1 states. Ordinary nonnegative
  relation rows are precisely the eventually-zero paths, so the finite
  carry system recognizes formal compatibility but not termination. Sharply,
  the primitive thirteen-speed rows {1,...,11,13^N,13^N+1} have one fixed
  support-three, coefficient-height-one relation and repeat the same
  carry/owner/current-digit state for N-1 levels. The scalar radix length of
  the quotient decreases by exactly one per level and is an exact
  termination clock. This is a carrier theorem and no-go for a uniform
  state-count horizon, not a proof of LRC(14).
source: codex-2026-07-25-relation-carry-clock
depends_on: []
related:
  - THM-2163-radix-relation-carry-descent
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
script: 04-computation/lrc14_relation_carry_completion_radix_clock_thm2231.py
output: 05-knowledge/results/lrc14_relation_carry_completion_radix_clock_thm2231.out
script_sha256: 244232a057dc078734da9006c4452f22b3a847d2ec7a84339efc47a347e5f876
output_sha256: 2d31ed3a58869093149301773bbb5a8c1042353c9eeafcb9f1f5949bf3544435
hash_basis: working-tree bytes (LF)
---

# THM-2231 -- relation-carry completion and the exact radix clock

THM-2163 proves that every finite integer relation row gives a bounded
zero-to-zero radix carry path.  The converse question has two different
answers depending on the ambient object.  The recurrence alone describes a
relation in the radix completion.  An LRC row needs the stronger condition
that its digit word actually ends.

This theorem identifies that distinction exactly and gives a sharp primitive
family showing that no bound on the number of digit levels follows from the
carry, owner mask, and current digit.

## 1. The completed relation-carry theorem

Fix

```text
q>=2,                     0!=m=(m_1,...,m_d) in Z^d.       (1)
```

Write

```text
Zhat_q=lim_(J) Z/q^J Z                                      (2)
```

for the `q`-adic completion.  This notation is used for every integer radix,
including composite `q`; no field structure is asserted.  Every
`V in Zhat_q^d` has a unique coordinatewise digit expansion

```text
V=sum_(j>=0) q^j D_j,
D_j in {0,...,q-1}^d.                                    (3)
```

Call the digit stream **carry-admissible** when there are integers
`kappa_j` satisfying

```text
kappa_0=0,
q kappa_(j+1)=kappa_j+m.D_j             for every j>=0.   (4)
```

> **Completed carry theorem.**  The digit expansion (3) gives a bijection
>
> ```text
> {carry-admissible digit streams}
>        <--> ker(m:Zhat_q^d -> Zhat_q).                    (5)
> ```
>
> The carries are unique.  If `C=||m||_1`, every carry obeys
>
> ```text
> |kappa_J|<C.                                             (6)
> ```
>
> More sharply, if
>
> ```text
> P=sum_(m_i>0)m_i,             N=sum_(m_i<0)(-m_i)        (7)
> ```
>
> are both positive, then
>
> ```text
> -N<kappa_J<P.                                            (8)
> ```
>
> Thus the mixed-sign carry alphabet contains at most
>
> ```text
> P+N-1=||m||_1-1                                         (9)
> ```
>
> integer states, independently of digit depth.

### Proof

For a digit stream put

```text
V_<(J)=sum_(j=0)^(J-1) q^j D_j.                          (10)
```

Multiplying (4) by `q^j` and summing from `0` through `J-1`
telescopes:

```text
m.V_<(J)
 =sum_(j<J)q^j(q kappa_(j+1)-kappa_j)
 =q^J kappa_J.                                           (11)
```

If the carries are integral, (11) says

```text
m.V_<(J)=0 mod q^J                                      (12)
```

for every `J`.  The compatible residues of `V_<(J)` are exactly `V`, so
`m.V=0` in the inverse limit.  This proves the forward direction of (5).

Conversely, suppose `m.V=0` in `Zhat_q` and take the unique digits (3).
For every `J`, the kernel condition gives

```text
m.V_<(J)=0 mod q^J.
```

Define

```text
kappa_J=(m.V_<(J))/q^J in Z.                           (13)
```

Since

```text
V_<(J+1)=V_<(J)+q^J D_J,
```

equation (13) gives exactly (4).  It also proves uniqueness: (11) determines
every carry from the digit prefix.

Each coordinate of `V_<(J)/q^J` lies in `[0,1)`.  Therefore

```text
|kappa_J|
 =|m.(V_<(J)/q^J)|
 <sum_i |m_i|=C,                                      (14)
```

which proves (6).  If both signs occur, separating positive and negative
coefficients gives (8).  The integers strictly between `-N` and `P` are

```text
-N+1,-N+2,...,P-1,
```

exactly `P+N-1` values.  This proves (9) and the theorem. QED.

## 2. Ordinary rows are the eventually-zero paths

Embed `Z_(>=0)` in `Zhat_q` by its ordinary base-`q` expansion.  Then

> **Termination classification.**
>
> ```text
> V in Z_(>=0)^d and m.V=0
> ```
>
> corresponds under (5) exactly to a carry-admissible path for which
>
> ```text
> D_j=0              for every sufficiently large j.      (15)
> ```

Indeed, a nonnegative integer has an eventually-zero base-`q` expansion, and
every eventually-zero expansion is a nonnegative integer.  If (15) begins at
`J_0`, recurrence (4) gives

```text
kappa_(J_0)=q^r kappa_(J_0+r)             for every r>=0. (16)
```

The fixed integer on the left is divisible by every power of `q`, so it is
zero.  Thus the path is eventually the zero carry as well.  Conversely, (5)
puts the finite digit vector in the completed kernel.  If an ordinary integer
is divisible by every `q^J`, it is zero, so the same vector lies in the
ordinary integer kernel.

There is an exact owner formulation.  For an arbitrary digit stream define
the tail-activity mask

```text
O_J(D)={i:there is some r>=J with D_(r,i)!=0}.             (17)
```

On an ordinary row, with

```text
Z_J=floor(V/q^J)              coordinatewise,             (18)
```

this is precisely

```text
O_J={i:Z_(J,i)>0}={i:V_i>=q^J}.                            (19)
```

Consequently

```text
ordinary nonnegative relation row
 iff O_J is eventually empty.                              (20)
```

The finite carry alphabet in (6)--(9) therefore recognizes the completed
relation kernel.  Eventual emptiness is a separate termination condition.

## 3. Primitive fixed-height pumping inside the LRC arity

The distinction in Section 2 persists under every elementary normalization
relevant to the support-three LRC relation branch.  For `N>=2`, put

```text
V^(N)=(1,2,...,11,13^N,13^N+1) in Z_(>0)^13             (21)
```

and use the fixed relation

```text
m=(1,0,...,0,1,-1) in Z^13.                              (22)
```

Thus the three nonzero entries of `m` occur at coordinates `1,12,13`.
Every row (21) has thirteen distinct positive coordinates, is primitive
because it contains `1`, and satisfies

```text
m.V^(N)=1+13^N-(13^N+1)=0.                              (23)
```

The relation has

```text
support(m)=3,          ||m||_infinity=1,       ||m||_1=3. (24)
```

For the base-`13` quotient at every level `1<=j<=N`,

```text
Z_j=(0,...,0,13^(N-j),13^(N-j)),
R_j=(1,2,...,11,0,1).                                  (25)
```

Hence

```text
kappa_j=(m.R_j)/13^j=0,
O_j={12,13}.                                            (26)
```

At levels `1<=j<N`, the current digit is

```text
D_j=0^13,                                               (27)
```

while at level `N` its last two entries are `1,1`.  Finally,

```text
Z_(N+1)=0,             O_(N+1)=empty.                  (28)
```

Therefore:

```text
(kappa_j,O_j) repeats for N consecutive levels;
(kappa_j,O_j,D_j) repeats for N-1 consecutive levels.   (29)
```

At the fixed level `j=1`, every `V^(N)` has exactly the same full local state
in (27), but its owner mask needs `N` more quotient steps to become empty.
As `N` is arbitrary, no finite-valued function of

```text
(carry, owner mask, current digit)                       (30)
```

can uniformly bound the remaining radix depth.  This remains true at fixed
dimension thirteen, after primitive normalization and distinctness, and for
one fixed support-three coefficient-height-one relation.

The family (21) is a carrier hostile control.  No assertion is made that its
members are LRC counterexamples.

## 4. The exact scalar repair

The full quotient `Z_J` in (18) restores termination, but a single unbounded
integer already suffices.  Put

```text
M_J=||Z_J||_infinity,

h_J=0                              if M_J=0,
h_J=1+floor(log_q M_J)             if M_J>0.             (31)
```

Equivalently, for `M_J>0`, `h_J` is the unique positive integer such that

```text
q^(h_J-1)<=M_J<q^h_J.                                  (32)
```

Since Euclidean quotienting gives

```text
Z_(J+1)=floor(Z_J/q)             coordinatewise,
M_(J+1)=floor(M_J/q),                              (33)
```

equations (31)--(32) imply the exact recursion

```text
h_(J+1)=max(h_J-1,0).                                  (34)
```

Moreover,

```text
h_J=0 iff O_J=empty.                                   (35)
```

Thus `h_J` is exactly the number of further radix divisions required to
empty the owner mask.  It is a scalar termination clock: it records only the
remaining digit length, not the quotient values, labels, or relation
geometry.

The clock is deliberately unbounded.  Equation (29) proves why replacing it
by another finite state count cannot supply a uniform all-height horizon.

## 5. Cross-frontier interpretation

The theorem separates three logically different coordinates:

```text
carry recurrence:      compatibility in the completed relation kernel;
owner mask:            whether a coordinate has any future nonzero digit;
radix clock h_J:        how far ordinary termination remains.             (36)
```

This is the relation-carry analogue of two established mechanisms.

1. THM-2225's cyclic checksum decides the fair output only after the initial
   critical-run length selects an adaptive dyadic block.  The response image
   controls balance; the run length controls stopping.  Here the carry
   controls the relation response and `h_J` controls stopping.
2. THM-2228 separates a compatible `2`-adic orbit from stabilization to an
   ordinary integer.  Equations (5), (15), and (34) give the corresponding
   exact split for arbitrary radix relation paths.

The source-to-target contract is:

```text
source:       a finite or infinite coordinatewise radix word;
quotient:     its bounded affine relation carry;
preserved:    m.V=0 in the q-adic completion;
destroyed:    eventual digit support and remaining magnitude;
sidecars:     O_J for termination truth, h_J for exact remaining depth;
hostile test: the primitive family V^(N).                         (37)
```

## 6. Boundaries

This theorem proves no finite height bound and does not make LRC(14) finite.
It also does not say that a finite automaton cannot recognize an
eventually-zero infinite word when supplied an acceptance condition.  The
proved no-go is narrower and exact: the current finite state (30) supplies no
uniform numerical horizon.

The family (21) does not obstruct:

- a search performed inside an independently proved height box;
- a descent carrying the full quotient `Z_J`;
- a proof with a separate decreasing magnitude potential; or
- a target whose natural objects already live in the `q`-adic completion.

It does obstruct deriving an all-height termination bound merely by counting
carry, owner, and current-digit states.

## 7. Exact referee

Run

```text
python 04-computation/lrc14_relation_carry_completion_radix_clock_thm2231.py
python -O 04-computation/lrc14_relation_carry_completion_radix_clock_thm2231.py
```

The companion uses only exact integer arithmetic and explicit raising checks.
It exhausts `826,353` finite cylinder vectors across radices `2,3,4`, with
`q=4` as a composite-radix control, and checks equivalence between terminal
divisibility and every integral carry transition.  It checks all accepted
cylinders against (8) and (11).

Independently, it verifies (21)--(29) for every `2<=N<=128`, including
`8,255` carry/owner plateau states and `8,128` full
carry/owner/current-digit plateau states.  It checks `8,509` clock
transitions on those rows and another `26,751` transitions for radices
`2..17`.  Ordinary and optimized runs are byte-identical and match the stored
output. QED.
