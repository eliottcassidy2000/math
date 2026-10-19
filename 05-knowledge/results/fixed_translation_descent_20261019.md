# Fixed-translation isolated collisions: a subdyadic height bound and smaller-source receipts

**Status, 2026-10-19:** PROVED elementary extensions, self-checked by exact
source controls. No assertion of universal positive-integer convergence or
universal equal-clock merging is made.

## 1. Inheritance and the question now answered

The [integer stabilizer packet](integer_stabilizer_20261019.md) proves that
isolated positive collisions T^N(x+d)=T^N(x), with unequal odd counts, satisfy
x<d/2+2^N/12. Its sharpness family allows d to grow proportionally to 2^N;
it leaves the **fixed translation d**, especially d=1, open.

This addendum proves two genuinely stronger source-wise statements:

1. For each fixed positive integer d, every isolated integer root at depth N
   is bounded by d/2+2^{N−c sqrt(N)+O_d(1)}, with the explicit constant
   c=1/sqrt(2 log_2 3). Thus the isolated integer head is o(2^N).
2. Outside the linear core x<(N+d)/2, every such isolated merge forces an
   actual smaller positive source z<x and clocked receipts from **both**
   initial sources to z. An odd-count-aware version is stronger still.

The first proof uses the nonvanishing integer carry absent from generic
real or Haar-only arguments. The second uses the actual positive-orbit
multiplicative correction, retaining the branch counts and clocks. We
inherit the affine pair chain and its measure boundary from
[THM-4569](../../01-canon/theorems/THM-4569-the-terras-clock-recurrence-of-two-collatz-orbits-is-unconditional.md)
and [THM-4581](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md),
and the receipt interpretation from the
[clock holonomy calculus](collatz_clock_holonomy_20261007.md).
The canonical hostile is (1,2), which never merges at equal clocks. The
near miss is applying a nonzero-carry argument after a whole-cylinder
identity hit. The least-used sidecar is the exact 3-adic valuation of the
translation, not another odd-modulus census.

Throughout T is the Terras map on positive integers. Fix x,d,N∈Z_{>0}, put
u_t=T^t(x+d), v_t=T^t(x), and suppose u_N=v_N with unequal final odd counts
a≠b. Let a_t,b_t be the counts through time t, and put K=max(a,b).

## 2. Isolated integer histories never have zero translation

Write the actual branch numerators as

    u_t=[3^{a_t}(x+d)+C_u(t)]/2^t,
    v_t=[3^{b_t}x+C_v(t)]/2^t.

If a branch has at least one odd step, its C is a 3-adic unit: its last
odd update is C↦3C+2^j, and subsequent even steps do not change C. Let

    k_t=a_t−b_t,   e_t=u_t−3^{k_t}v_t.

Then, whenever k_t≠0,

    v_3(e_t)=min(0,k_t).                            (2.1)

**Proof.** If k_t>0, then

    2^t e_t=3^{a_t}d+C_u(t)−3^{k_t}C_v(t)

is a 3-adic unit, because a_t≥1 and C_u(t) is the unique unit term. If
k_t<0, multiply this identity by 3^{−k_t}; the term −C_v(t) is the unique
unit. ∎

If k_t=0 and e_t=0, the two literal endpoints have already merged with
identity state; all their subsequent parities agree, so k_N=0. This
contradicts a≠b. Hence **e_t never vanishes at any time 0≤t≤N**.

Define the integer normalized carry

    E_t=3^{max(0,b_t−a_t)}u_t
          −3^{max(0,a_t−b_t)}v_t.                   (2.2)

It is nonzero throughout the isolated history. When k_t≠0 it is even a
3-adic unit. This latter unit statement is not needed for the height bound,
but explains why a pure nonidentity dilation cannot appear along an integer
translation history.

## 3. A depth bound from odd counts and common-even blocks

**Lemma 3.1.** At every time t≤N,

    0<|E_t|<(d+1)3^K.                              (3.1)

**Proof.** Put H_u=C_u/3^{a_t} and H_v=C_v/3^{b_t}. The previous packet's
sharp offset bound gives 0≤H_u,H_v≤2^t/6 for t≥1. Directly from (2.2),

    E_t=3^{max(a_t,b_t)}[d+H_u−H_v]/2^t.

Therefore |E_t|≤3^K(d/2^t+1/6)<(d+1)3^K. At t=0, E_0=d gives the same
bound. Nonvanishing is §2. ∎

Let L=floor(log_2((d+1)3^K)). Every consecutive block during which both
paths are even has length at most L: during that block E is halved at each
step, and a nonzero integer cannot be halved more than L times from a
starting absolute value less than (d+1)3^K.

There are at most a+b≤2K−1 steps at which at least one path is odd, since
a≠b. They separate at most 2K common-even blocks. Consequently

    N+1≤2K [1+floor(log_2((d+1)3^K))].              (3.2)

This is a deterministic bound for actual integer sources. It does not use
independence, a random parity driver, or a statement about typical paths.

## 4. The fixed-d isolated head is subdyadic

A refinement of the previous height proof is

    x<d/2+2^{N−K}.                                 (4.1)

Indeed, writing m=min(a,b), D=|a−b| and q=3^D, the branch collision
formula bounds x by d/2+2^{N−m}/(q−1), and q−1≥2^D. If a>b the d/2
term can be omitted. The strict inequality follows from the positive d or
the positive lower-branch affine offset, as in the previous packet.

Define the exact integer threshold

    K_d(N)=min {k≥1:
        N+1≤2k[1+floor(log_2((d+1)3^k))]}.          (4.2)

Combining (3.2) and (4.1) proves the fully explicit bound

    x<d/2+2^{N−K_d(N)}.                             (4.3)

The right side can be evaluated using integer bit lengths; no floating
logarithm is required in a certificate checker.

For a closed-form estimate, put A=log_2 3 and B_d=1+log_2(d+1). Equation
(3.2) implies 2AK²+2B_dK≥N+1, hence

    K ≥ [sqrt(B_d²+2A(N+1))−B_d]/(2A)
      ≥ sqrt((N+1)/(2A))−B_d/(2A).                 (4.4)

It follows that

    x < d/2 + 2^{N−sqrt((N+1)/(2 log_2 3))
                        +(1+log_2(d+1))/(2 log_2 3)}.   (4.5)

For every fixed d, the maximum isolated positive integer root at depth N,
with an empty maximum taken as zero, is therefore o(2^N). The gain is an
explicit factor 2^{−c sqrt(N)+O_d(1)}, c=1/sqrt(2 log_2 3).

For d=1, the exact thresholds K_1(N) at N=64,256,1024,4096,16384 are
5,9,18,36,72. These are theorem-derived bounds, not enumerations of roots
at those depths. They are not claimed to be optimal.

**Why the earlier sharpness family is consistent.** In that family d is
of order 2^N and K=2. Then B_d is of order N, so (3.2) permits arbitrary
N with bounded K. Fixed d is an essential quantifier. No bound o(2^N)
uniform over d≤2^N is asserted.

## 5. Exact positive correction products force a smaller source

For the actual positive trajectories, define

    P_u=product over odd u_t, t<N, of (1+1/(3u_t)),
    P_v=product over odd v_t, t<N, of (1+1/(3v_t)).

The branch formulas in multiplicative form give

    u_N=(3^a/2^N)(x+d)P_u,
    v_N=(3^b/2^N)x P_v,
    P_v/P_u=3^{a−b}(x+d)/x =: R.                  (5.1)

All quantities are positive exact rationals except the optional logarithms
below. This is the classical telescoping product applied to a pair, not a
claim that residue data determine the real correction factors.

**Lemma 5.1 (rational product threshold).** If m positive factors have the
form 1+1/(3z_i) and their product is at least r>1, then

    min_i z_i < mr/[3(r−1)].                        (5.2)

**Proof.** Otherwise each factor is at most 1+t/m, with t=1−1/r∈(0,1).
The binomial theorem gives

    (1+t/m)^m ≤ Σ_{j=0}^m t^j/j! < Σ_{j=0}^∞ t^j=1/(1−t)=r,

a contradiction. ∎

When R>1, apply the lemma to P_v≥R; when R<1, apply it to P_u≥1/R.
When R=1, this comparison supplies no forced small value. The exact
unified thresholds are

    some odd lower value < bR/[3(R−1)]      if R>1;
    some odd upper value < a/[3(1−R)]       if R<1.  (5.3)

The ratio R=1 can occur, for example for (45,15) at its isolated merge.
It is a genuine boundary of this sufficient test, not a failure of literal
coalescence. Additional directly observed descents may still exist there.

### A simple count-aware descent criterion

Suppose x>d/2. If a<b, then 3^{b−a}≥3 and (5.1) gives P_u≥3x/(x+d)>1.
If a>b, it gives P_v≥3(x+d)/x>1. Thus

    some odd upper value z < a x/(2x−d),       a<b;
    some odd lower value z < b(x+d)/(2x+3d),  a>b.  (5.4)

Let m=min(a,b). In either case

    m≤2x−d  ⇒  an actual premerge value z satisfies 0<z<x.   (5.5)

For the second case, the verification is
(2x−d)(x+d)<x(2x+3d); the first case is immediate. In particular,

    N≤2x−d  ⇒  both initial sources have smaller-source receipts.   (5.6)

To make the receipt explicit, choose the **first actual value below x**
on the path whose correction product supplies the small value. If this
value is z at time j, then 0<j<N, and

    T^N(x+d)=T^N(x)=T^{N−j}(z).                     (5.7)

One path also reaches z directly at time j. These are actual clocked
receipts to the same smaller positive source. If z is already certified,
both original sources are certified; no prior convergence assumption for
x itself is needed. Both clocks are retained, so the argument does not
silently discard the phase obstruction.

The logarithmic refinement, when desired, follows from log(1+t)<t:
replace (5.2) by min z_i<m/(3 log r). The rational form above is weaker
but exact and convenient for an integer certificate checker.

## 6. A sharper division of work for the integer atlas

At any fixed depth N and translation d:

- every isolated positive integer root lies in the subdyadic head (4.3);
- every such root with x≥(N+d)/2 has a proved receipt to a source below both
  initial integers;
- only the linear core x<(N+d)/2 escapes that depth-only sufficient descent
  rule; actual odd counts or the exact ratio (5.3) often certify more.

This is a classification of **existing merge receipts**, not proof that
all integer sources merge. It reduces a genuinely new obstruction: a
coalescence edge between two unresolved sources can be replaced by two
edges to a smaller source whenever (5.5) holds. It supplies the rank decrease
that a bare common endpoint lacks, without treating a Haar-null atom as
Haar mass.

The finite atlas remains incomplete as a global convergence mechanism.
There may be sources for which no receipt has been proved, and equal-clock
coverage of d=1 is itself false at x=1. The results give no bound on when
an arbitrary integer first enters either the cylinder or isolated atlas.
A phase-aware covering theorem remains OPEN.

## 7. Focused exact controls and reproduction

The [checker](../experiments/fixed_translation_descent_20261019.py) reuses the frozen atlas
from [the preceding packet](integer_stabilizer_20261019.json), rather than
running another residue census. It performs **12,865 explicit checks**;
the [published normal receipt](fixed_translation_descent_20261019.json) matches its optimized-mode replay.

It replays 183 isolated-source/depth rows, containing 34 distinct first
births. Every normalized carry, 3-adic valuation, even-block bound and exact
product identity is checked along the actual integer paths. Among those 34
births, the linear criterion certifies 15, the odd-count criterion 16, and
the exact rational product criterion 27. These are finite sufficiency
counts, not convergence densities. The remaining rows are not declared
counterexamples: some merely lie outside the sufficient test, including
the already certified lower source 1 and the R=1 boundary.

The checker also verifies all ten changing-d sharpness examples at
N=6,12,...,60 against the fixed-d formula, and checks the identity-cylinder
and (1,2) hostile cases. Large-depth rows are evaluations of the integer
threshold (4.2), not new source experiments.

```sh
python3 05-knowledge/experiments/fixed_translation_descent_20261019.py > /tmp/fixed-d.json
python3 -O 05-knowledge/experiments/fixed_translation_descent_20261019.py > /tmp/fixed-d-O.json
cmp /tmp/fixed-d.json /tmp/fixed-d-O.json
cmp /tmp/fixed-d.json 05-knowledge/results/fixed_translation_descent_20261019.json
```
