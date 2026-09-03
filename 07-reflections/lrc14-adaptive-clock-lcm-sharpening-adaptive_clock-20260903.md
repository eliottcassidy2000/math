# LRC(14): one minimal wall LCM is already a complete adaptive owner clock

Status: **PROVED AS THM-4371 + VERIFIED-EXACT ARITHMETIC; INTEGRATED.**  This
sharpens only the explicit clock in THM-4349.  It does
not prove that the owner relation is empty for every body and does not prove
LRC(14).

## 1. Exact candidate

Let `H` be an eleven-element set of distinct positive integers.  Put

```text
h=max H,
O_h={a: 1<=a<12h and a is odd},
D_H=lcm({2c:c in H} union O_h),
Q_H=7D_H.                                                (1)
```

There is a useful exact normal form.  If

```text
L_h^odd=lcm(O_h),               nu_H=max_{c in H} v_2(c),
```

then every odd part of every `c in H` already occurs in `L_h^odd`, and hence

```text
D_H=2^(1+nu_H) L_h^odd.                                 (1a)
```

Thus, at fixed height, the body dependence of the clock is only its largest
`2`-adic core scale; all odd prime powers are forced already by the finite
tail box.

Then the four equivalent conditions in THM-4349 are also equivalent to the
single explicit condition

```text
R_(Q_H)(H)=empty.                                        (2)
```

Thus THM-4349's displayed clock

```text
28(42h+1)^2 lcm(1,...,12h)                              (3)
```

may be replaced by the body-sensitive clock `(1)`.  This removes the
quadratic padding, cancels the factor two in every doubled wall endpoint,
and also removes the formerly retained final factor two between witness and
certificate clocks.

There are two useful height-only forms.  Define

```text
D_h=lcm({2,4,...,2h} union O_h).                         (4)
```

Equivalently,

```text
D_h=2^(1+floor(log_2 h)) L_h^odd.
```

For every eleven-body of height `h`, the clock `7D_h` is complete.  The
simpler, slightly larger statement

```text
Q_tilde_h=7 lcm(1,...,12h-1)                            (5)
```

is also complete.  Indeed `D_H|D_h|lcm(1,...,12h-1)`, and THM-2066 divisor
transport makes emptiness propagate from a clock to all its multiples.

These clocks are still astronomical.  They are uniform completeness clocks,
not practical first-certificate estimates.

## 2. Inheritance and connection contract

The inherited mechanisms are:

- THM-2061's strict metric box: an unsafe odd pair has `a,b<12h`;
- THM-2066's exact labelled relation and divisor transport;
- THM-4070's continuum spoiled-set equivalence; and
- THM-4349's closed safe arc, rational endpoint witness, and long-block
  extraction.

The corrected near miss is `(3)`: all of its extra padding was sufficient,
not necessary.  The canonical hostile remains the max-`72` body that escapes
clocks `15,...,43` but closes at `47`; adaptivity remains real even after
improving this enormous universal bound.

The live concept board was

```text
doubled wall | physical speed species | minimal wall LCM | closed safe arc
long block | centered residue | owner top bit | adjacent owner rigidity.   (6)
```

THM-4363/4365/4367 are an analogy, not a dependency: a quotient can discard
the scale or address needed by its next consumer.  Here centering modulo
`Q_H` appears to lose one top bit of a tail class modulo `2Q_H`.  Unlike the
ray-scale collisions in those theorems, this loss is repaired internally:
on two adjacent eligible packet points the centered owner is constant, while
the missing top bit alternates.  Complementarity therefore forces the two
tails to have the same top bit, which cancels.

The connection contract is

```text
source:    rational safe endpoints for all bounded actual tail pairs
target:    one finite labelled owner relation R_(Q_H)(H)
map:       y=2x, common wall denominator, centering, adjacent owner rigidity
preserves: closed safety, strict danger, both lift labels, odd parity
apparently loses: the top bit modulo 2Q_H
repair:    two consecutive core-safe residues force equal top bits
test:      isolated top-bit flip followed by its adjacent-point failure.   (7)
```

## 3. Doubling gives the minimal common wall clock

Assume every distinct positive odd pair `a,b<12h` gives a safe row

```text
S(H;a,b)=2H union {a,b}.                                 (8)
```

For a fixed pair, its full safe-time set is a nonempty closed finite union of
intervals and points.  It is not the whole circle because time zero is
unsafe.  Hence some component endpoint lies on the safety wall of a physical
speed `m in S(H;a,b)`:

```text
x=(14k+/-1)/(14m).                                      (9)
```

The body phase used by the owner packet is `y=2x`.  Therefore

```text
y=(14k+/-1)/(7m),                                      (10)
```

whose reduced denominator divides `7m`, rather than `14m`.  The possible
physical speed species are exactly

```text
m=2c with c in H,              or              m in O_h. (11)
```

Thus every bounded pair has a safe lift over an unspoiled body phase on the
`Q_H`-grid, because

```text
Q_H=lcm({14c:c in H} union {7a:a in O_h})=7D_H.         (12)
```

This is the least common wall multiple furnished by the endpoint argument.
If the safe set is a singleton, its wall is still a valid witness: safety is
closed, whereas tail danger and the spoiled set are strict.  No positive
width assumption is used.

## 4. The wall LCM itself clears the long-block threshold

Cited lower-dimensional LRC gives THM-4349's closed body-safe arc

```text
I={y:dist(y,y_0)<=1/(84h)},              |I|=1/(42h).   (13)
```

At a clock `Q`, its grid points contain a cyclic block of `K` consecutive
residues with

```text
K>=Q/(42h)-1.                                           (14)
```

Because `I subset G(H)`, every residue in this block belongs to the labelled
core-safe packet `A_Q(H)`.

If a tail class `z mod 2Q` is eligible, subtracting its strict `Q/7` bounds
at the first grid point from the next `K-1` points gives

```text
|jz|_Q<2Q/7,                       1<=j<=K-1.           (15)
```

The elementary no-wrap lemma from THM-4349 says that, for the centered
representative `q` of `z modulo Q`, `(15)` implies

```text
|q|<2Q/[7(K-1)]
   <=2Q/[7(Q/(42h)-2)].                                (16)
```

The final expression is below `12h+1` whenever

```text
Q>1008h^2+84h.                                         (17)
```

No additional factor is needed to meet `(17)`.  The two largest members of
`O_h` are coprime:

```text
gcd(12h-1,12h-3)=1.
```

Moreover `D_H` is even because it contains every doubled body speed.
Consequently

```text
D_H is divisible by 2(12h-1)(12h-3),
Q_H=7D_H
 >=14(12h-1)(12h-3)
 =(1008h^2+84h)+42(24h^2-18h+1)
 >1008h^2+84h.                                         (18)
```

The last inequality is valid for every integer `h>=1`; an eleven-body has
`h>=11`.  Since `Q_H` is even, reducing an odd class modulo `Q_H` preserves
odd parity.  Equations `(16)--(18)` therefore give

```text
q odd and |q|<12h+1,              hence |q|<=12h-1.     (19)
```

Also `(14)` contains many consecutive points: `(17)` gives
`K>24h+1`.  In particular the top-bit argument below has two adjacent
core-safe packet residues.

## 5. Adjacent owner rigidity cancels the missing top bit

Assume the bounded rows in `(8)` are all safe, and suppose for contradiction
that

```text
(u,v) in R_Q(H),                    Q=Q_H.               (20)
```

Both odd classes are individually eligible by the definition of `R_Q`.
Let `a,b` be their centered representatives modulo `Q`.  By `(19)`, they are
nonzero odd integers with

```text
|a|,|b|<12h.                                            (21)
```

Because `u,v` live modulo `2Q`, there are unique top bits
`epsilon,delta in {0,1}` such that

```text
u=a+epsilon Q mod 2Q,
v=b+delta Q mod 2Q.                                    (22)
```

We now type the top-bit step exactly.  Extend a packet label cyclically to an
integer `s`; replacing its canonical representative by `s+Q` flips both odd
owner bits (each nearest integer increases by an odd number), so complementary
ownership is unchanged.  Choose two consecutive lifted labels `s,s+1` in the
block `(14)`.  For `q=a` or `q=b`, write

```text
n_q(t)=nint(qt/Q),
e_q(t)=qt/Q-n_q(t).                                    (23)
```

Eligibility at both labels gives `|e_q(s)|,|e_q(s+1)|<1/7`.  Hence

```text
n_q(s+1)-n_q(s)=q/Q-e_q(s+1)+e_q(s),
|n_q(s+1)-n_q(s)|<1/2+2/7=11/14<1.                    (24)
```

The left side is an integer, so it is zero.  Thus both centered owner bits
are constant between these adjacent labels.  On the other hand `(22)` gives

```text
omega_u(t)=n_a(t)+epsilon t mod 2,
omega_v(t)=n_b(t)+delta t mod 2.                        (25)
```

Complementarity at `s` and `s+1`, subtracted modulo two, now yields

```text
epsilon+delta=0 mod 2,
so                         epsilon=delta.               (26)
```

The common top-bit term in `(25)` cancels at every packet label.  Eligibility
also depends only on reduction modulo `Q`.  Therefore

```text
(a,b) mod 2Q is itself in R_Q(H).                       (27)
```

If `|a|!=|b|`, the `Q_H`-grid witness from Section 3 for the distinct actual
tails `|a|,|b|` leaves one lift unspoiled.  Signs preserve literal danger
masks, contradicting `(27)`.  If `|a|=|b|`, the two signed tails have
identical masks and cannot cover both lifts over any nonempty packet point;
the block `(14)` supplies such points.  This also contradicts `(27)`.

Thus `R_(Q_H)(H)` is empty.  Conversely, emptiness of any owner relation
certifies safety for all actual odd tails by THM-2066.  Together with the
finite cutoff and continuum equivalence already proved in THM-2061/4070/4349,
this proves the candidate in Section 1.  THM-2066 divisor transport then
gives the height-only clocks `(4)--(5)`.

Repeated residue classes cause no gap.  They remain allowed in `R_Q`, as
required because distinct actual tails can coincide modulo `2Q`.  After
centering, equal absolute values have identical literal masks and are
excluded using the nonempty safe block, not by silently assuming distinct
residue classes.

## 6. Hostile: the top bit flips at one point but fails on the next

The top-bit concern is real at one packet point.  At height `11`, take the
height-uniform clock `Q=7D_11` and the core-wall scale `c=8`.  Then

```text
v_2(14c)=v_2(Q),             r=Q/(14c) is odd.          (28)
```

At phase `r/Q=1/112`, the classes `1` and `1+Q` have the same centered
residue modulo `Q`, but their literal two-lift masks are

```text
at r:       (danger,safe),       (safe,danger).         (29)
```

So centering at one point genuinely forgets the owner bit.  At the adjacent
phase `(r+1)/Q`, however, both masks are

```text
at r+1:     (danger,safe),       (danger,safe).         (30)
```

The apparent complementary pair has disappeared.  This is exactly the
alternation in `(25)`, and explains why the two-point argument is necessary.
The control is a literal mask hostile, not a claim that this particular wall
belongs to `A_Q(H)` for every body.

## 7. Exact audit

The companion checks:

- wall-LCM divisibility and the strict threshold for every `1<=h<=256`;
- all `31,824` eleven-bodies with `11<=h<=18` for body/height divisibility,
  every possible doubled wall, and the strict threshold;
- `40,600` exact doubled-wall fractions;
- the no-wrap lemma in `2,566,368` finite cases;
- adjacent centered-owner constancy in `880,594` finite cases;
- the isolated top-bit flip and adjacent-point repair `(28)--(30)`; and
- the max-`72` body-specific versus height-uniform clocks.

Reproduce with

```text
python3 -B 04-computation/lrc14_adaptive_clock_lcm_sharpening_probe_adaptive_clock_20260903.py
python3 -B -O 04-computation/lrc14_adaptive_clock_lcm_sharpening_probe_adaptive_clock_20260903.py
PYTHONHASHSEED=947 python3 -B 04-computation/lrc14_adaptive_clock_lcm_sharpening_probe_adaptive_clock_20260903.py
```

All three streams equal
`05-knowledge/results/lrc14_adaptive_clock_lcm_sharpening_probe_adaptive_clock_20260903.out`.
The final raw-LF SHA-256 values are

```text
6eb8b09a135aa98b241562f2480d8dad9ef282a77abf536674e27ec07d3d47a9  script
b591050f7ebd533de8e12f69a0721b08ea62cea6323432a252b6eb4d6f674efe  output
```

## 8. Exact scope and next sharp question

The theorem candidate sharpens an **if and only if certificate**.  It does
not assert that `R_(Q_H)(H)` is always empty, classify unsafe bodies, improve
the inherited `12h` tail cutoff, or lower the logical LRC(14) frontier.
THM-4349 remains a completeness theorem, not an empty-seam theorem.

Within the endpoint-synchronization architecture, `Q_H` is the exact common
multiple in `(12)`.  This does not prove global numerical optimality: a
different argument could use witnesses away from endpoints or a smaller
non-divisibility clock.  The next sharp problem is whether one can exploit
the adjacent-owner rigidity structurally at practical clocks, rather than
only after the enormous LCM has synchronized every bounded-pair witness.  For
the max-`72` hostile, `Q_H` still has `378` digits and is
`1,171,280,000` times smaller than THM-4349's old explicit clock, whereas
an exhibited empty clock is only `47`; the practical-clock gap remains vast.
