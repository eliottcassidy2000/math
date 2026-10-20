# Two decompositions of ten and the involutions they retain

**PROVED:** the explicit divisor-coordinate and involution statements below.
The classification of `F=U+S`, the ten-letter map, and the obstruction to a
geometric complement action are inherited. **FINITE-EXACT:** the independent
integer-divisor and permutation controls in `arithmetic_check.py`. No
literature-priority, tournament-count, or Collatz consequence is claimed.

## 1. Inheritance, conventions, and the literal prime-shape question

The closest proved mechanism is [DB1--DB3, arithmetic braids](../arithmetic_braids_20260917_divisors.md):
divisors form an exponent box, and squarefree divisors form its Boolean
subcube. The canonical hostile is `p^2`, where `1 != 1+1`. The corrected
near miss is the identification of the squarefree and rank `3,4,3`
partitions, separated in [the tournament/clock note, section 4.3](../procgen_tcpc_20261001_tournament_clock_prime_collatz.md).
The least-used sidecar is the **full complement/prime-exchange action** in
[divisor layers, sections 2--4](../divisor_layers_pair_geometry_20261004.md).

The local board is exponent box / squarefree excess / rank / complement /
prime exchange / geometric pair carrier. The anchor is the two requested
decompositions; the niche is their common symmetry group; the wildcard is
which parts of that group a tetrahedral model can retain.

For `N>=2`, let `D(N)={d:d|N,1<d<N}`. Define `F=|D|`, `S` as the number
of squarefree elements of `D`, and `U` as the number of prime elements of
`D`. Thus the prime set is contained in the squarefree set. `F=U+S` is an
equality of counts, not the literal disjoint union of those two subsets.

For `N=prod p_i^a_i` with `r` distinct primes, the inherited formulas are

    F=prod(a_i+1)-2,
    S=2^r-1-[N squarefree],
    U=r-[N prime].

They give the complete inherited classification

    F=U+S  iff  N=p, p^3, or p^2 q r,

with three distinct primes in the last shape. For completeness, a
nonsquarefree solution needs `prod(a_i+1)=2^r+r+1`. At `r=1,2,3`, this
forces respectively `(3)`, no profile, and `(2,1,1)`; at `r>=4`,
`prod(a_i+1)>=3*2^(r-1)>2^r+r+1`. Squarefree composites fail because
`F=S` and `U>0`; primes have all three counts zero.

The expression `{p^q r,p,p^3}` must retain its literal typing. If `p,r`
are distinct primes and `q>=2` is an integer exponent, then `p^q r` has
`(F,S,U)=(2q,3,2)`, so equality would require `2q=5`. It is never a
solution. At `q=1`, the squarefree product also fails. The recovered
three-shape classification uses `p^2 q r`, with `q` a separate prime
factor. These two strings describe different families.

## 2. Several 3-4-3 partitions are already different

Fix distinct primes `p,q,r`, and let `N=p^2qr`. A divisor is uniquely
`d=p^a q^b r^c`, where `0<=a<=2` and `b,c in {0,1}`, after deleting
`(0,0,0)` and `(2,1,1)`. The ten proper divisors have the following
different partitions:

| Observable | First class | Middle class | Last class |
|---|---|---|---|
| prime / squarefree composite / nonsquarefree | `p,q,r` | `pq,pr,qr,pqr` | `p^2,p^2q,p^2r` |
| rank `a+b+c=1,2,3` | `p,q,r` | `p^2,pq,pr,qr` | `pqr,p^2q,p^2r` |
| repeated-prime exponent `a=0,1,2` | `q,r,qr` | `p,pq,pr,pqr` | `p^2,p^2q,p^2r` |
| sign of `b-c` | `r,pr,p^2r` | `p,p^2,qr,pqr` | `q,pq,p^2q` |

Every row has sizes `3,4,3`. The first two differ by exchanging `p^2`
and `pqr`. The first row does not become the second by an invariant
arithmetic identification. In the first row, the count equality is
realized by the inherited bijection

    prime ell -> p*lcm(p,ell)

from the three prime divisors to the three nonsquarefree proper divisors.
It gives a genuine repaired disjoint decomposition: the squarefree set
plus a **tagged copy** of the prime set maps bijectively to `D(N)`.
Removing the tag would identify prime divisors already present in the
squarefree set.

## 3. A precise 4-2-4 reconstruction on the same divisors

The bounded source search recovered no definition identifying the user's
`4-2-4` with a particular partition. Here is an explicit construction that
keeps the same ten divisors and the inherited symmetries.

Let `C(d)=N/d` and let `Q` exchange prime labels `q,r`, fixing `p`.
In exponent coordinates,

    C(a,b,c)=(2-a,1-b,1-c),
    Q(a,b,c)=(a,c,b),
    CQ(a,b,c)=(2-a,1-c,1-b).

These are commuting involutions. `C` reverses divisibility, `Q` preserves
it, and `CQ` reverses it. They generate a Klein-four action on the
proper divisor set. This action was already inherited; the following
two observables make the requested splits explicit:

    f=a-b-c,                  g=b-c.

The sign partition of `f` is

    f<0: {q,r,qr,pqr},
    f=0: {pq,pr},
    f>0: {p,p^2,p^2q,p^2r}.                       (1)

This is exactly `4-2-4`. Its complete height distribution, at
`f=-2,-1,0,1,2`, is `1,3,2,3,1`. The middle two elements are exactly the
fixed points of `CQ`. Its four exchanged pairs are

    p <-> pqr,    q <-> p^2q,    r <-> p^2r,    qr <-> p^2.

Thus (1) selects one member of each nontrivial `CQ` orbit on each side,
without an arbitrary transversal. The gauge is the distinguished repeated
prime direction, and the convention that it has positive sign.

The sign partition of `g` is the fourth `3-4-3` row above. Its middle
four elements are exactly the fixed points of `Q`; its three exchanged
pairs are `q<->r`, `pq<->pr`, and `p^2q<->p^2r`. Choosing which of `q,r`
is positive exchanges its two outside classes.

The transformation laws display the difference between the observables:

| Involution | `f` | `g` | Number of fixed divisors |
|---|---:|---:|---:|
| `C` | `-f` | `-g` | 0 |
| `Q` | `f` | `-g` | 4 |
| `CQ` | `-f` | `g` | 2 |

The preserved predicate is a specified sign or a specified symmetry
orbit. A sign alone forgets the actual divisor, rank, divisibility and
magnitude. For instance `p,p^2,p^2q,p^2r` have the same positive sign
but different prime support and divisibility relationships. The full
triple `(a,b,c)` plus the prime labels restores those data.

## 4. The even-exponent boundary explains the two central points

For all integers `k>=2`, replace `N` by `p^kqr`. Retain the same exponent
box and omit bottom and top. It has `4k+2` points. Define

    C(a,b,c)=(k-a,1-b,1-c),
    Q(a,b,c)=(a,c,b),
    f_k=2a-k(b+c),                 g=b-c.

The same sign transformation table holds. There are `2k` `Q`-fixed
points: those with `b=c`, after deleting the two endpoints. The two
outside sign classes of `g` each contain `k+1` points. Thus its split is

    (k+1), 2k, (k+1).                                 (2)

For `CQ` to fix a point requires `a=k/2` and `b+c=1`.
Hence its fixed set has two points when `k` is even and none when `k`
is odd. These are also exactly the points with `f_k=0`: the apparent
alternatives `(a,b,c)=(0,0,0),(k,1,1)` were removed. `CQ` exchanges the
positive and negative sets, so their sizes are forced:

    sign(f_k): 2k,2,2k                when k is even;
               2k+1,0,2k+1           when k is odd.     (3)

At `k=2`, equations (2) and (3) give `3-4-3` and `4-2-4` on one set.
The central pair in (3) is a parity boundary, not a fitted count. In
contrast, the inherited rank layers are `3,4,...,4,3`, and squarefree
balance `F=U+S` remains true only at `k=2`: its defect is `4(k-2)`.

### The two full heights together are a faithful planar chart

Although either sign partition is lossy, the pair `(f_k,g)` recovers
every proper divisor exactly. Its image consists of three rows:

    g=0:    f_k=-2k,-2k+2,...,-2, 2,4,...,2k;
    g=+1:   f_k=-k,-k+2,...,k;
    g=-1:   f_k=-k,-k+2,...,k.                          (4)

For `g=+1`, recover `(b,c)=(1,0)`; for `g=-1`, recover `(b,c)=(0,1)`.
In either case `a=(f_k+k)/2`. When `g=0`, a positive `f_k` means
`(b,c)=(0,0)` and `a=f_k/2`; a negative `f_k` means `(b,c)=(1,1)` and
`a=(f_k+2k)/2`. The origin is absent.

Equivalently, the linear projection from three exponent coordinates
has kernel generated by `(k,1,1)`. Within the full exponent box, only
bottom and top collide under that projection, both at the origin.
Removing exactly those two endpoints makes it injective. Thus the
proper-divisor convention that enters `F=U+S` also makes this planar
chart faithful. Its reflection actions are

    C:(f_k,g)->(-f_k,-g),
    Q:(f_k,g)->(f_k,-g),
    CQ:(f_k,g)->(-f_k,g).

At `k=2`, the row counts are `3,4,3`; splitting the same chart by the
sign of its other coordinate gives `4,2,4`. This is an explicit shared
geometry with exact decoding. Divisibility is recovered by comparing
the decoded exponents coordinatewise; it is not ordinary planar
coordinatewise order. The prime labels remain necessary to recover
the numeric divisor. The smallest hostile to forgetting the deleted
endpoints is the collision `1,N -> (0,0)`.

## 5. The ten-letter geometric map has an exact symmetry obstruction

The inherited [duck decoder, section 2](../duck_decoder_20260925.md)
maps these ten divisors bijectively to unordered pairs with repetition
from the four-element set `V={0,P,Q,R}`, and then to `E(K5)` by
`{u,v}->uv` for distinct letters and `{v,v}->v-infinity` for repeated
letters. Its table sends

    p,q,r -> {0,P},{0,Q},{0,R};
    pq,pr,qr -> {P,Q},{P,R},{Q,R};
    pqr -> {0,0};
    p^2,p^2q,p^2r -> {P,P},{Q,Q},{R,R}.

This retains a labelled ten-element set and its specified `7+3` split.
It does not preserve divisor order or multiplication. In this exact
table `Q` is geometrically realized by exchanging the tetrahedral
corners `Q,R`. However `C` sends `p` to `pqr`, exchanging a midpoint
with a corner, and `CQ` fixes the two adjacent edge midpoints
`{P,Q},{P,R}`. Neither is a tetrahedral affine symmetry in this table.

There is also a table-independent obstruction to making the whole
Klein-four action geometric. An involution of four corners has `t`
transpositions and `s=4-2t` fixed corners. Its induced fixed unordered
pairs with repetition number

    s(s+1)/2+t = 10,4,2.

None has zero fixed pairs, whereas `C` has none. Therefore no bijection
to the tetrahedral corners and midpoints intertwines the full divisor
action with affine symmetries. An arbitrary set bijection can match
`CQ` alone with a double transposition, because both have cycle type
`1^2 2^4`; it cannot at the same time make `Q` geometric, because their
product would make `C` geometric. This is the exact missing-coordinate
boundary behind matching decompositions of ten.

The three simple-graph types on four vertices and three edges are a
different three-element quotient: their labelled fibre sizes are
`12,4,4`, not ten. Their mere count of types does not identify them with
the three prime shapes or with either partition above. An additional
map would have to specify graph incidence, labels, gauge and predicate.

## 6. Reproduction and hostile controls

Run from the repository root:

    python -B 05-knowledge/results/trijection_reset_20261020/arithmetic_check.py
    python -B -O 05-knowledge/results/trijection_reset_20261020/arithmetic_check.py

The checker compares numeric divisors found by trial division with the
exponent model for three marked prime assignments and `2<=k<=8`; tests
all group laws, fixed sets, sign partitions and complement-order reversal
there; checks the all-`k` formulas through `k=40`; and independently
enumerates the corner involutions of the tetrahedron. Positive controls
include `60` and `p^3`; hostile controls include `p^2`, the literal
two-prime shape `p^q r`, and `p^4qr`. Normal and optimized output agree
byte for byte; the receipt is [arithmetic_check.json](arithmetic_check.json).
These controls supplement the structural proofs above.
