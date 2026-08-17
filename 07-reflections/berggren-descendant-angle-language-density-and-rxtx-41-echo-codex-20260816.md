# Berggren angle languages: the exact `16/41,9/41,16/41` split and the RXTX 41 echo

**Status: PROOF-COMPLETE ELEMENTARY CANDIDATE + VERIFIED-EXACT;
independent audit pending.**  The exact companion is
`04-computation/berggren_descendant_angle_language_density_20260816.py`.
This note advances THM-3334's pointwise angle classification to an all-level
tree language and density theorem.  It proves no tournament current, LRC
statement, Jacobian statement, or RXTX tensor identity.

## 1. Inheritance and the two walls

THM-2596 gives the three Berggren maps on the reduced Euclid slope
`x=m/n in (0,1)`:

```text
A(x)=1/(2-x),       B(x)=1/(2+x),       C(x)=x/(2x+1).   (1)
```

Their images are respectively `(1/2,1)`, `(1/3,1/2)`, and `(0,1/3)`.
They form a `PGL_2(Z)` reduction cross-section, not a cyclic order-three
action.  On coprime opposite-parity pairs they generate the primitive
Pythagorean tree rooted at `x=1/2`.

THM-3334 proves that the triangle formed by the three descendants of the
parent `(a,b,c)` is U-obtuse, acute, or D-obtuse according as

```text
b/a < 8/9,       8/9 < b/a < 9/8,       b/a > 9/8.      (2)
```

Since

```text
a=n^2-m^2,       b=2mn,       b/a=2x/(1-x^2),           (3)
```

and the last function is strictly increasing, the two slope walls are

```text
alpha=(sqrt(145)-9)/8,
beta =(sqrt(145)-8)/9.                                  (4)
```

They satisfy

```text
4alpha^2+9alpha-4=0,       9beta^2+16beta-9=0,
0<alpha<beta<1.                                         (5)
```

Every tree slope is rational, so neither equality wall is ever attained.

## 2. The unexpected period-four mechanism

Let `T` be the inverse branch selected by the three intervals in (1):

```text
T(x)=x/(1-2x) on (0,1/3),
T(x)=1/x-2     on (1/3,1/2),
T(x)=2-1/x     on (1/2,1).                              (6)
```

Exact arithmetic in `Q(sqrt(145))` gives two pure four-cycles:

```text
alpha --B--> (sqrt(145)-7)/8
      --A--> (17-sqrt(145))/12
      --B--> (sqrt(145)-7)/12
      --B--> alpha,                                     (7)

beta  --B--> (sqrt(145)-10)/9
      --C--> sqrt(145)/29
      --B--> -2+sqrt(145)/5
      --B--> beta.                                      (8)
```

Thus the two boundary words are the periodic words `BABB` and `BCBB`.
This is the structural source of every recurrence below.

## 3. Exact level counts

For a rational starting slope `x_0`, let

```text
F_n(t;x_0)=3^(-n) #{w in {A,B,C}^n : w(x_0)<t}.         (9)
```

The disjoint branch intervals give, for every distribution supported on
rational slopes,

```text
t in (0,1/3): F_(n+1)(t)=      F_n(Tt)/3,
t in (1/3,1/2): F_(n+1)(t)=2/3-F_n(Tt)/3,
t in (1/2,1): F_(n+1)(t)=2/3+F_n(Tt)/3.                 (10)
```

The minus sign in the middle row records that `B` reverses orientation.
Composing (10) around (7)--(8) yields

```text
F_(n+4)(alpha;x_0)=32/81-F_n(alpha;x_0)/81,
F_(n+4)(beta ;x_0)=50/81-F_n(beta ;x_0)/81.             (11)
```

Crucially, (11) holds for **every rational starting slope**, not only the
root.  This is the uniform cylinder estimate needed later.

At the root `x_0=1/2`, let `U_n,A_n,D_n` count U-obtuse, acute, and D-obtuse
parents at depth `n`.  Then

```text
U_(n+4)=32*3^n-U_n,       U_0..U_3=(0,1,4,10),
A_(n+4)=18*3^n-A_n,       A_0..A_3=(0,1,2,6),
D_(n+4)=32*3^n-D_n,       D_0..D_3=(1,1,3,11).          (12)
```

Solving the affine recurrences gives

```text
U_n/3^n -> 16/41,
A_n/3^n ->  9/41,
D_n/3^n -> 16/41.                                      (13)
```

The convergence is much stronger than an asymptotic estimate.  The exact
residual triples

```text
(41U_n-16*3^n, 41A_n-9*3^n, 41D_n-16*3^n)              (14)
```

are periodic with period eight:

```text
(-16,-9,25), (-7,14,-7), (20,1,-21), (-22,3,19),
(16,9,-25), (7,-14,7), (-20,-1,21), (22,-3,-19).       (15)
```

So the count error is bounded at every level.  The first direct rows are

```text
n       0   1    2      3       4        5          6
U_n     0   1    4     10      32       95        284
A_n     0   1    2      6      18       53        160
D_n     1   1    3     11      31       95        285.  (16)
```

The companion directly enumerates all `88,573` nodes through depth ten and
then checks (12)--(15) through depth forty.

## 4. Regular languages and harmonic subseries

Each threshold language is regular.  Read an ancestry word from leaf to
root.  The undecided state consists of one of the four points in (7) or (8)
and an inequality orientation.  A branch whose image lies wholly on one side
of the current wall enters an accepting or rejecting sink.  The unique branch
containing the wall advances the phase, and `B` flips the inequality.  This
is a finite automaton with at most eight undecided states and two sinks per
wall.  Regular languages are closed under reversal, and the product of the
two wall automata gives the three angle languages in root-to-leaf order.

Now enumerate the ternary tree in shortlex order, using any fixed order of
`A,B,C`, and let `S_U,S_A,S_D` be the three subsets of positive indices.
Equation (11), applied after an arbitrary fixed prefix, says that a suffix
cylinder of size `3^r` contains

```text
delta_T*3^r+O(1)                                       (17)
```

members of type `T`, uniformly in the prefix.  A lexicographic initial
segment of one level is a disjoint union of at most two prefix cylinders at
each depth.  Combining (17) over those cylinders and over the completed
earlier levels gives

```text
#{m<=N:m in S_T}=delta_T*N+O(log N),                    (18)
```

where

```text
(delta_U,delta_A,delta_D)=(16/41,9/41,16/41).            (19)
```

Partial summation therefore gives constants `C_T` such that

```text
sum_(m<=N,m in S_T) 1/m
  =delta_T log N+C_T+O((log N)/N).                      (20)
```

This makes the user's harmonic-series observation nontrivial: the labelled
subseries retains the exact angle language, while its logarithmic leading
coefficient remembers only one of the three rational densities in (19).
The coefficient is independent of the chosen fixed letter order; the
constant `C_T` need not be.

The raw consecutive-Fibonacci slope `F_k/F_(k+1)` is at least `1/2` for
`k>=2`, but this is not always the canonical primitive Euclid parameter.
When `k=1 mod 3`, both entries are odd and THM-3339 requires

```text
T(m,n)=((n-m)/2,(n+m)/2),       x'=(n-m)/(n+m)<=1/3.
```

Thus the normalized Fibonacci three-ray locus has periodic chamber word

```text
D,D,U,D,D,U,... .                                      (20a)
```

It is a sparse path in both obtuse languages, not the source of either
global density.  MISTAKE-418 records the repaired normalization boundary.

## 5. The lawful size-four object is a signed partial tournament

The two four-periodic wall orbits have the same orientation-sign word:

```text
alpha: B A B B,       beta: B C B B,
sign:   - + - -,      product=-1.                       (21)
```

Here the sign is the slope of the affine CDF update in (10): `B` reverses
the order and `A,C` preserve it.  The four phases and their successor edges
therefore form an oriented `C4`.  Of the six unordered pairs on four
vertices, four are observed transition edges and the two antipodal pairs are
missing.  If inverse transitions are retained as well, the same carrier is a
four-edge bidirected cycle with two missing pairs.  This is precisely a
directed graph with missing/both-way edges; it is not yet a tournament.

There are four ways to orient the two missing diagonals and complete the
oriented cycle to a tournament.  Exact exhaustion shows that none is
invariant under phase rotation.  Indeed the half-turn fixes each antipodal
pair setwise and swaps its endpoints, so an invariant arrow on such a pair
is impossible.  A complete `T4` would therefore add a gauge absent from the
recurrence.

Modulo two, the negative-edge indicators in (21) are a cocycle on `C4` with
odd cycle sum.  Vertex sign changes add coboundaries, while the odd sum is
the nonzero class in

```text
H^1(|C4|;F2)=F2.                                         (22)
```

This sign holonomy explains both the plus sign in `3^4+1` and the residual
period eight.  If `R_n` is any of the three centered count residuals in
(14), then

```text
R_(n+4)=-R_n,             R_(n+8)=R_n.                  (23)
```

Equivalently the ordinary generating functions are

```text
U(z)=(z+z^2-2z^3+2z^4)/((1-3z)(1+z^4)),
A(z)=(z-z^2)/((1-3z)(1+z^4)),
D(z)=(1-2z+2z^3-z^4)/((1-3z)(1+z^4)).                  (24)
```

The pole at `z=1/3` carries the densities; the roots of `1+z^4` carry the
eight-periodic bounded discrepancy.  This is an explicit word-to-`H^1`
map, but only for the wall automaton.  It does not furnish the missing LRC
same-copy edge or a characteristic-zero JC flux class.

## 6. The exact RXTX denominator echo is not a transfer

Rybin--Zhang--Luo's RXTX preprint has the exact cost recurrence

```text
R(n)=8R(n/4)+26M(n/4),       M(n)=n^(log_2 7).           (25)
```

Thus its asymptotic coefficient solves

```text
r=(8/49)r+26/49,       r=26/(49-8)=26/41.               (26)
```

The Berggren U-density solves

```text
u=32/81-u/81,          u=32/(81+1)=16/41.               (27)
```

The common denominator is the exact integer identity

```text
3^4+1=82=2(7^2-2^3).                                   (28)
```

This is a useful **resolvent-denominator echo**: both constants come from
closing a recursive equation after four-way structure is exploited.  It is
not a common tensor, tournament, or branch action.  RXTX's `41` is the cost
gap `49-8`; the angle `41` is half the orientation-reversing return
denominator `81+1`.  There is no map from RXTX's 26 general products to the
three angle languages, and none is needed for (12)--(20).

The paper's genuine size-four content is its symmetric `4 x 4` block output,
with four diagonal and six off-diagonal positions.  The Berggren language's
period four instead belongs to a quadratic wall orbit.  Equal numerals do not
identify those carriers.

## 7. Boundaries

- The three angle classes color tree nodes; they do not orient the six edges
  of a `T4` or `T6`.
- The rational densities are ancestry-shortlex densities, not densities by
  hypotenuse height or by Farey denominator.
- The shortlex harmonic series uses word indices.  Reciprocal triple
  coordinates along a Fibonacci or parabolic ray are different series.
- Nothing here supplies an LRC current, a D5 class, a Jacobian flux, or a
  matrix-multiplication lower bound.

## Reproduction

```text
python -B 04-computation/berggren_descendant_angle_language_density_20260816.py
python -B -O 04-computation/berggren_descendant_angle_language_density_20260816.py
```

Normal and optimized transcripts agree.  Current LF-normalized hashes are
recorded in the matching result file.
