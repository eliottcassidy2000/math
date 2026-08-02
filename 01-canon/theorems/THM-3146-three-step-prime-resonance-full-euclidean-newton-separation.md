---
id: THM-3146
title: "Three-step prime resonance full Euclidean-Newton separation"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT AUDIT.  For every
  odd prime p, no exact quadratic has three zero factorial moments beginning
  at r=p+1.  At resonance d=p+3, subtracting the full linear Euclidean
  quotient gives a remainder with generic slopes 2/(p+1),2/(p-1), disjoint
  from the first polynomial's slopes 0,2/p.  The exact endpoint residues 87
  and 37 isolate p=29 and p=37; p=37 remains separated by a four-vertex
  polygon, while p=29 and p=3 lie in THM-3124's proved finite range.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3138-left-factorial-resonance-newton-slope-separation
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
script: 04-computation/factorial_three_step_prime_resonance_euclidean_newton_thm3146.py
output: 05-knowledge/results/factorial_three_step_prime_resonance_euclidean_newton_thm3146.out
script_sha256: 90482ac88ee20831a5384d189282531cb98a63e864f6c5250920a78166f46d07
output_sha256: 61ab8000639af13475a2025a1fc3c457c9bed49d51dc7716dfc01c33faa74cd4
hash_basis: LF-normalized bytes
---

# THM-3146 -- three-step prime resonance full Euclidean-Newton separation

**PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT AUDIT.**

Let

```text
L(t^k)=k!,                    q(t)=a+bt+ct^2,                 (1)
```

with `abc!=0`.  For every odd prime `p`, put

```text
r=p+1,                       d=r+2=p+3.                      (2)
```

Then the three consecutive moments

```text
L(q^r), L(q^(r+1)), L(q^(r+2))                               (3)
```

cannot all vanish.

This is the next Euclidean-Newton layer after THM-3143.  Cancelling only the
leading term no longer separates the polygons; the full linear quotient does.

## 1. Resonant pair and its first polygon

THM-3124 says that `(3)` forces `b/a=-1/d`.  After division by `a`, write
`u=c/a`, set `v=d u`, and define

```text
A_n(v)=d^n L((1-t/d+u t^2)^n)
      =L((d-t+v t^2)^n) in Z[v].                             (4)
```

It is enough to prove coprimality of

```text
A=A_(p+1),                         B=A_(p+2).                 (5)
```

As before,

```text
a_(n,j)=[v^j]A_n
 =binom(n,j) sum_(ell=0)^(n-j)
  binom(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!.                (6)
```

Now `d==3 (mod p)`.  For `A`, the constant and linear coefficients are both
`6 (mod p)`.  Coefficients `2<=j<p` have at least one `p` from their binomial
factor; starting at `j=(p+1)/2`, their factorial gives a second.  The last
two indices have the same required lower bounds, and

```text
[v^(p+1)]A=(2p+2)!,                   v_p((2p+2)!)=2.        (7)
```

For every `p>=5`, this gives

```text
NP_p(A):(0,0)--(1,0)--(p+1,2),        slopes 0 and 2/p.      (8)
```

The small prime `p=3` is handled in Section 5.

## 2. The full linear quotient

The leading quotient coefficient is

```text
C=(2p+4)(2p+3).                                             (9)
```

Cancelling the next coefficient gives

```text
q=-2(p+2)(p+3)/(2p+1).                                     (10)
```

The denominator `D=2p+1` is a `p`-adic unit.  Clear it and define the
integral remainder

```text
R(v)=D B(v)-[D C v-2(p+2)(p+3)]A(v).                        (11)
```

Its terms of degrees `p+2` and `p+1` vanish exactly, so `deg R<=p`.
Multiplication by the unit `D` does not change Newton slopes or the
common-factor predicate.

The quotient formula follows directly from the four top coefficients:

```text
[v^(p+1)]A=(2p+2)!,
[v^p]A=(p+1)(2-p)(2p)!,
[v^(p+2)]B=(2p+4)!,
[v^(p+1)]B=-p(p+2)(2p+2)!.                                 (12)
```

## 3. Low coefficients and the midpoint unit

Lucas reduction of `(6)` gives the only nonzero low residues needed:

```text
A_0==A_1==6,
B_0==15,       B_1==0,       B_2==72                  (mod p), (13)
```

with all remaining relevant residues zero.  Since `D==1`, `C==12`, and
`2(p+2)(p+3)==12 (mod p)`, equation `(11)` gives

```text
R_0==15+12*6=87,               R_j==0 for every j>0  (mod p). (14)
```

Thus all positive-degree coefficients are divisible by `p`, while the
constant is a unit unless `p=3` or `p=29`.

Put

```text
m=(p-1)/2,                         k=m+1=(p+1)/2.             (15)
```

For `p>=5`, the only term surviving in `a_(p+1,m)/p` is `ell=0`, and

```text
a_(p+1,m)/p
 ==(-1)^(m-1)3^(m+2)/[m(m-1)] !=0                  (mod p).  (16)
```

At index `k`, the `B_k` and `A_k` contributions to `(11)` have valuation at
least two.  Hence

```text
R_k/p
 ==(-1)^m 12*3^(m+2)/[m(m-1)] !=0                 (mod p).  (17)
```

The midpoint coefficient has valuation exactly one.  For every `k<j<p`,
all three terms of `R_j` have valuation at least two.

## 4. High coefficients and the two arithmetic chambers

Direct simplification of the high coefficients gives

```text
R_p=-2p(p+1)(p+2)(24p+37)(2p-2)!,                          (18)

R_(p-1)=2p(p-1)(p+1)(p+2)
         (24p^2-115p-246)(2p-4)!.                           (19)
```

For `p!=37`, `(18)` has valuation exactly two.  If also `p!=29`, `(14)`,
`(17)`, and the intervening bounds force

```text
NP_p(R):(0,0)--(k,1)--(p,2),
slopes             2/(p+1),  2/(p-1).                       (20)
```

The prime `p=37` raises the last point to height three.  Formula `(19)` has
valuation exactly two there, and the exact replacement polygon is

```text
NP_37(R):(0,0)--(19,1)--(36,2)--(37,3),
slopes                1/19,    1/17,    1.                  (21)
```

These slopes are disjoint from `{0,2/37}`.

At `p=29`, the constant residue `87` vanishes exactly once.  The complete
integer polygon is

```text
NP_29(R):(0,1)--(15,1)--(29,2),       slopes 0 and 1/14.     (22)
```

Its zero slope overlaps `(8)`, so this observer no longer proves
coprimality.  This is an exact failure of the Euclidean-Newton certificate,
not a factorial counterexample.

## 5. Coprimality and finite exceptional closure

For every `p>=5` except `29`, the slope sets of `A` and `R` in
`(8),(20),(21)` are disjoint.  A common factor of `A` and `B` would divide
the exact combination `R`; after base change to `Q_p` it would contribute a
common lower Newton slope.  Therefore `A` and `B` are coprime for all those
primes.

The two remaining primes are already in THM-3124's independently audited
exact gcd range:

```text
p=3:   r=p+1=4,
p=29:  r=p+1=30.                                             (23)
```

THM-3124 excludes both resonant windows.  Hence `(3)` is impossible for
every odd prime.  QED.

## 6. Exact controls

The companion reconstructs the full integer polynomials and cleared
remainder for

```text
p=5,7,11,29,37,41,101,211.                                  (24)
```

It verifies the quotient, identities `(18),(19)`, all regular hulls, the
overlapping `p=29` hull, and the separated four-vertex `p=37` hull.  The
`p=41` control is deliberate: its penultimate factor in `(19)` acquires an
extra `p`, but that point lies above the unchanged endpoint segment and
creates no chamber.  A repeated-product bivariate constructor independently
checks the original polynomials for `p<=7`.  Finally, the closed residues
`87`, `37`, and `(17)` are scanned over every odd prime below `10,000`; their
only exceptional primes are exactly `29` and `37`.

Run

```text
python3 04-computation/factorial_three_step_prime_resonance_euclidean_newton_thm3146.py
python3 -O 04-computation/factorial_three_step_prime_resonance_euclidean_newton_thm3146.py
```

and compare byte-for-byte with the declared output.

## 7. Connection contract and scope

The source object is THM-3124's resonant pair `(5)`.  The map is full linear
Euclidean reduction `(11)` followed by the neighboring-prime Newton polygon.
It preserves common factors and destroys root residue/phase.  The required
sidecars are the constant residue `87`, midpoint unit `(17)`, and high
identities `(18),(19)`.  The `p=29` hostile identifies the first lost
coordinate exactly; THM-3124 supplies the finite fallback rather than an
invented slope argument.

Combining THM-3131, THM-3138, THM-3143, this theorem, and THM-3124, any
still-open bad exact-quadratic window must satisfy

```text
r>=201,
d=r+2, d-1, d-2, and d-3 are all composite.                  (25)
```

The resonance therefore lies after a prime gap of length at least four.
This suggests an iterated subresultant/Newton resolution indexed by the
prime gap, but no arbitrary-gap induction is asserted here.

This is an exact `{0,1,2}` / `SFC(1)` factorial-moment theorem.  It does not
settle arbitrary-support `SFC(3)`, `GMC(2)`, `NC(2)`, or `LRC(14)`.

**End of proof candidate.**
