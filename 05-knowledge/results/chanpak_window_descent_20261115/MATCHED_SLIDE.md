# Sliding a fixed strong core makes an unbounded matched window

A proper return need not keep the core at a fixed position. The
arithmetic carry can be realized by sliding the marked strong
component along a strip of transitive vertices. The exact output
is an unbounded arithmetic window, not just congruence coverage.

## Slide theorem

Let $C$ be a fixed marked tournament with ordinary path count $h$,
two-slot count $Q$, total marked starts $s$, and total **marked**
ends $r$. For a transitive word of marks, retain its width and read
its bits in tournament order. Suppose two words $z,w$ of the same
positive width $d$ obey
$$
 (2^d-1)s+hz=Qw,\qquad
 Q\,\operatorname{rev}_d(z)=(2^d-1)r+h\,\operatorname{rev}_d(w).
$$
Then $Cz$ and $wC$ have equal full ordinary count, two-slot count,
marked starts and unmarked ends, without asserting equality of
the signed start/end channels. This is the two-ended product
factorization proved in [the five-core note](FIVE_CORE_WIDTH.md);
the same proof uses arbitrary $h,Q$, with $s,r$ taken on $C$.

Put $p$ equal to the last bit of $w$ and $q$ equal to the first
bit of $z$. The one-vertex modules with marks $p,q$ act as guards.
Then for every $k\ge0$ the $k+1$ literal marked bases
$$
 B(k,j)=p\;w^j\;C\;z^{k-j}\;q,\qquad 0\le j\le k ,
$$
have the same full switch response $(h,Q_2,E_s,O_s,E_e,O_e)$
at fixed $k$. At each slide, the two unsigned totals match and
transport unchanged through arbitrary common outer ordinal blocks.
The guards are the ultimate first and last vertices for every $j$,
so the signed start is $h$ times the last mark and the signed end
is $h$ times one minus the first mark.

Moreover, if a marked base's apex count is written $N$, then
$$
 N(B(k,j))=N(B(k,0))-j\,\Delta
$$
where $\Delta=N(pCzq)-N(pwCq)$. To check this in the presence of
arbitrarily many other tiles, compare $Cz$ with $wC$ and keep the
exterior parts fixed. All unsigned effects in the cut transfer
agree. The only remaining exterior effects are the central
signed end difference multiplied by the left signed start, and
the central signed start difference multiplied by the right signed
end. The last bit of the left exterior always equals $p$ (whether
the guard or a tile), and the first bit of the right exterior
always equals $q$. The ordinary counts of both exteriors are one.
Thus the count difference is exactly the one measured with two
guards alone. The same proof is valid if either of the two signs
vanishes. This argument is the endpoint transport in the
[cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md),
not equality of the unmarked tournaments.

Consequently $|\Delta|=2$ would give $k+1$ consecutive odd
outputs with full matched profiles. More generally $\Delta=2e\ne0$
gives a local digit interval of length $k+1$ in steps of $e$ in
the half-gap coordinates. The existence of a return with nonzero
output is an additional obligation; the all-order width theorem
implies none exists for $h=5$.

## A nonzero physical return

Take the strong five-vertex core with outgoing rows
$(6,12,24,1,11)$, in bit convention vertex $i$ is bit $2^i$.
Mark only vertex 1. Literal core enumeration gives
$$
 h=13,\quad Q=104,\quad s=r=18,\quad N(C)=29.
$$
Use $z=101000010010$, $w=010000000111$. Both guards are marked.
The two integer identities above hold and the literal
recurrence gives
$$
  N(1Cz1)=382625,\qquad N(1wC1)=382605.
$$
Thus $\Delta=20$. From the left-guarded marked $C$ we start
with $(N,S,d)=(113,122,2)$. Appending a $z$ tile updates
$S$ to $4096S+33514$ and $N$ to $N+2938S+24076$:
apply the singleton marked-cut rule to the 12 displayed bits.
The first bit is marked, so the block's contribution to $N$
is independent of the incoming $d$. The final guard is marked
and adds no $N$. Solving the affine recursion proves that the
window at each $k\ge1$ is the exact arithmetic progression
$$
 \frac{9267808\cdot4096^k+1944617}{99225}
 +\frac{9776k}{315}-20j,\qquad j=0,\ldots,k.
$$
The coefficients are rational shorthand for the integer
recursion. The ordinary bases have order $5+12k+2$, and the
apex adds one more vertex. Each fixed-height palette is a
matched source bank, not a claim that every output is strong.

The central difference $20$ has half-gap $10$ coprime to $h=13$.
At each $k\ge12$ the physical palette therefore realizes all
local residue digits modulo 13; at $k\ge25$ it realizes two
consecutive lifts of each residue **in its own step-10
coordinates**, not two adjacent lifts in the original
step-one half-gap lattice. This distinction matters in CRT:
using only these switches leaves a global step-10 defect.
The strategy is to seek an all-height return of unit half-gap,
or find further matched returns and a literal coupling that
removes the common step. Positions or heights alone do not
place the exponentially separated absolute baselines.

The [replay](replay_matched_slide.py) enumerates the core and its
single-guard profiles with literal subset paths, then applies
the exact transfer along the actual marked words to check
full-response equivalence, the difference and the closed form.

## Two coprime slide steps for composition

The following additional tiles satisfy the same two integer tests.
Replay computes the signed-corrected step with the stated guards
and confirms both complete responses are equal at each test height.

| Ordinary core count | Outgoing rows | Mark mask | $z$ | $w$ | Guarded counts at height one | $\Delta$ |
|---:|---|---:|---|---|---|---:|
| 13 | $(6,12,24,1,11)$ | 2 | 101001000010 | 010000001101 | 395585, 395575 | 10 |
| 29 | $(70,8,2,5,15,31,62)$ | 49 | 0001000011010010111100100010 | 1001010001000011000010110011 | 425498252879, 425498252873 | 6 |

The guards are $(1,1)$ in the first row and $(1,0)$ in the
second row. In a composition with the first bank before the
second, the global first vertex is marked and global last vertex
unmarked. Hence the common apex reaches every base vertex and
every base vertex reaches it, and every composite is strong.
The corresponding half-gap increments after gluing are
$5\cdot29=145$ and $3\cdot13=39$, coprime. The slide theorem
makes their capacities arbitrarily large, while each core keeps
its ordinary count. The two-denomination bounded carry, rather
than CRT residue membership alone, is what converts these
returns into contiguous windows.
