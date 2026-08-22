---
id: THM-3159
title: "Right-neighbor exceptional-prime reciprocal-root zero-face separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  At the
  right-neighbor prime q=249727, the two resonant quadratic factorial-moment
  polynomials have distinct positive q-adic slopes.  Their height-zero faces
  are coprime by an exact reciprocal-root descent to a degree-249727 extended
  gcd over F_q.  Consequently the window r=q-3=249724 is null-free.  The
  theorem does not claim the analogous finite-field gcd for every prime q.
audit: >
  The first candidate audit caught an odd partial-fraction sign error on the
  reflected cubic pole.  MISTAKE-355 records the false quartic transform.
  A fresh independent audit rederived the coefficient valuations and Newton
  slopes, Wilson sign, corrected partial fractions and determinant, the
  reciprocal polynomial P(x)=2xU(x)+1, both singular charts, and the exact
  degree-249727 gcd and Bezout identity.  Normal and optimized replays match
  the stored LF-normalized transcript and declared hashes.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
  - THM-3153-four-step-prime-resonance-second-euclidean-newton-separation
script: 04-computation/factorial_right_neighbor_exception_reciprocal_root_thm3159.py
output: 05-knowledge/results/factorial_right_neighbor_exception_reciprocal_root_thm3159.out
script_sha256: 7b57f5df8969a8f78cdb54f8b6fe05865a2b817553c83bc386e8aa82740110bc
output_sha256: cbce9ea29197c5e40aff38f8f7a707bf415650db285dcdb406173816ff08c77d
hash_basis: LF-normalized bytes
---

# THM-3159 -- right-neighbor exceptional-prime reciprocal-root separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(t^j)=j!,                         Q(t)=a+bt+ct^2            (1)
```

with `abc!=0`, and put

```text
q=249727,                          r=q-3=249724.             (2)
```

Then the three consecutive moments

```text
L(Q^r), L(Q^(r+1)), L(Q^(r+2))                              (3)
```

cannot all vanish.

The proof has two genuinely different parts.  The `q`-adic Newton polygon
separates every nonzero slope.  A reciprocal-root chart then turns the two
large height-zero faces into one linear-size reverse truncated exponential,
where an exact extended gcd disposes of the remaining unit roots.

## 1. Resonant integral pair

By THM-3124, `(3)` would force `b/a=-1/d` with

```text
d=r+2=q-1.                                                   (4)
```

Divide by `a`, put `v=d(c/a)`, and define

```text
A_n(v)=L((d-t+v t^2)^n) in Z[v].                            (5)
```

It is enough to prove that

```text
A=A_(q-3),                         B=A_(q-2)                 (6)
```

are coprime.  Their coefficients are

```text
[v^j]A_n
 =binom(n,j) sum_(ell=0)^(n-j)
  binom(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!.               (7)
```

## 2. The two Newton polygons

Write `m=(q-1)/2`.  For `n=q-c`, with `c=3` or `2`, every coefficient
strictly above index `m` is divisible by `q`: every factorial in `(7)` then
has argument at least `q+1`.  The coefficient at `m` is a unit because only
the `ell=0` summand survives modulo `q`, and the leading coefficient

```text
[v^n]A_n=(2n)!                                                (8)
```

has valuation exactly one.  The exact constant residues and midpoint
residues are

```text
A_(q-3)(0) == -112210,       [v^m]A_(q-3)==156079,
A_(q-2)(0) ==   25306,       [v^m]A_(q-2)==124864   (mod q). (9)
```

All four are nonzero.  Integral coefficients cannot lie below height zero,
and every coefficient after `m` lies strictly above the open endpoint
segment.  Therefore

```text
NP_q(A):(0,0)--(m,0)--(q-3,1),   slopes 0, 2/(q-5),
NP_q(B):(0,0)--(m,0)--(q-2,1),   slopes 0, 2/(q-3).          (10)
```

The positive slopes are distinct.  A common factor over `Q`, after base
change to `Q_q`, could therefore have only slope zero.  It remains to exclude
a common nonzero root of the reductions of `A` and `B`.

## 3. Frobenius endpoint faces

In `F_q[t]`, set

```text
E(t)=sum_(j=0)^(q-1) t^j/j!,       H(t)=1-t-vt^2,
F_c(v)=[t^(q-1)] E(t)/H(t)^c.                              (11)
```

Wilson's theorem gives, for every polynomial `f`,

```text
L(f) == -[t^(q-1)]E(t)f(-t)                    (mod q).     (12)
```

Terms of degree at least `q` vanish on the left.  Since

```text
(d-t+vt^2)|_(t -> -t) == -H(t),
H(t)^q == 1                                      (mod t^q), (13)
```

equations `(11)--(13)` give the exact faces

```text
A_(q-c)(v) == (-1)^c F_c(v)                     (mod q),
c in {2,3}.                                                  (14)
```

Thus a common unit root in `(10)` would be a common root of `F_2,F_3`.

## 4. Reciprocal-root descent

Choose `alpha` in the algebraic closure with

```text
v=alpha(alpha-1),
H(t)=(1-alpha t)(1-(1-alpha)t).                              (15)
```

First suppose `alpha` is not `0,1,1/2`.  Define the reverse truncated
exponential

```text
U(alpha)=alpha^(q-1)E(alpha^(-1))
        =-sum_(k=0)^(q-1)(-1)^k k! alpha^k.                 (16)
```

For `k=1,2,3`, put

```text
S_k(alpha)=[t^(q-1)]E(t)/(1-alpha t)^k.                     (17)
```

The identities `E'=E+t^(q-1)` and

```text
tE''+(1-t)E'-E=0                                             (18)
```

give

```text
S_1(alpha)=U(alpha),
S_2(alpha)=-(U(alpha)+1)/alpha,
S_3(alpha)=(U(alpha)-alpha+1)/(2alpha^2).                   (19)
```

Apply the ordinary partial-fraction decompositions of `H^-2` and `H^-3`
and substitute `(19)`.  The odd pole on the reflected branch carries the
minus sign

```text
(beta*t-1)^(-3)=-(1-beta*t)^(-3).                           (20)
```

With that sign retained, `F_2=F_3=0` is a linear system for
`U(alpha),U(1-alpha)`.  After clearing the displayed denominators, its
determinant is `-2alpha(alpha-1)(2alpha-1)^3`, and solving gives

```text
U(alpha)=-1/(2alpha),
U(1-alpha)=-1/[2(1-alpha)].                                 (21)
```

Equivalently, every common root away from the three singular charts would
be a common root of

```text
P(x)=2xU(x)+1,                     P*(x)=P(1-x).             (22)
```

These are explicit polynomials in `F_q[x]` of degree `q=249727`.

## 5. Exact certificate and the singular charts

The companion constructs `(16),(22)` coefficient by coefficient in
`F_q[x]` and performs both a gcd and an extended gcd.  It obtains

```text
gcd(P,P*)=1,
S(x)P(x)+T(x)P*(x)=1.                                      (23)
```

and verifies `(23)` by exact polynomial multiplication.  No floating point,
random evaluation, or truncated subresultant is used.

The divisions in `(19)--(21)` excluded two values of `v`.  Direct exact
coefficient recurrences give

```text
(F_2(0),F_3(0))                 =(25306,112210),
(F_2(-1/4),F_3(-1/4))           =(191906,38381)   in F_q^2. (24)
```

Both charts are therefore disjoint from the common zero locus.  Equations
`(20)--(24)` prove `gcd(F_2,F_3)=1`.  Together with `(10),(14)`, this proves
that `A,B` are coprime over `Q`, contradicting `(3)`.  QED.

## 6. Exact companion

Run

```text
python 04-computation/factorial_right_neighbor_exception_reciprocal_root_thm3159.py
python -O 04-computation/factorial_right_neighbor_exception_reciprocal_root_thm3159.py
```

and compare byte-for-byte with the declared output.  Besides the
degree-`249727` gcd and verified Bezout identity, the companion checks
primality of `q`, all four Newton endpoints, both singular charts, and the
reciprocal-root transform at every odd prime from `5` through `47`.  Every
one of those corrected small transformed pairs is coprime.

## 7. Connection contract and scope

The source is THM-3124's exact quadratic resonant pair.  The first map is
the `q`-adic Newton functor, which preserves common factors but collapses all
unit roots to the height-zero face.  The second map is the two-sheeted
reciprocal-root chart `(15)`; it replaces a large symmetric polynomial in
`v` by the ordered pair `alpha,1-alpha`.  The preserved predicate is common
face vanishing.  The destroyed datum is the ordering of the two sheets; the
reflection `P*(x)=P(1-x)` is precisely the required sidecar.

The number `249721` is the largest eligible prime in THM-3148's offset-five
endpoint resultant.  Its window has

```text
d=249721+5=249726=q-1,       r=d-2=q-3,                     (25)
```

so this theorem closes that exact exceptional window.  It does not by itself
prove the generic positive-slope analysis for the entire offset-five family;
that Euclidean-Newton layer must be stated and audited separately.

Most importantly, the exact computation `(23)` is for `q=249727`.  Small
controls suggest the analogous `gcd(P_q,P_q*)=1` for every odd prime, but no
all-prime symbolic identity is claimed.  Nor does the theorem settle an
arbitrary prime gap, arbitrary-support `SFC(3)`, `GMC(2)`, `NC(2)`, or
`LRC(14)`.

## 8. Candidate-audit correction

The first pushed candidate used a quartic polynomial in place of `(22)`.
Its derivation converted the reflected cubic pole with the wrong sign.  The
independent audit caught that error before promotion; MISTAKE-355 records the
minimal failed implication.  Equations `(20)--(23)`, the companion, stored
output, and hashes were all replaced rather than silently retaining the old
certificate.

The repaired candidate then received a fresh independent hostile audit.  It
rederived `(7)--(10)`, including all four endpoint residues; checked the
Wilson sign in `(14)`; reconstructed the corrected partial-fraction system,
its determinant and solution `(21)`; and verified that `(22)` is exactly the
necessary off-chart condition.  It independently replayed the two singular
charts, the full gcd and extended-gcd certificate, and both normal and
optimized executions.  The repaired statement and evidence are therefore
promotion-safe; the superseded quartic candidate is not a dependency.
