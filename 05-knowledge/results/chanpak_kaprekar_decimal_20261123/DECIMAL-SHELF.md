# Five strong decimal children and literal Kaprekar seams

## An exact decimal append on a physical family

Let $X$ be any nonempty tournament, $D=H(X)$, and $Q=Q_2(X)$ (two ordered path slots, possibly empty). Add three ordinal singleton sources to $X$, and wrap the result in the directed sandwich $X'\to b\to c\to X'$. This is a strong tournament of height $h=8Q+D$. The shell outputs below are all strong and have exact heights
$$h_d=80Q+(10+d)D=10h+dD,\qquad d\in\{1,3,5,7,9\}.$$
In particular for a nonempty transitive filler of width $m$, $D=1$ and $Q=2^m$. Thus this **literally appends any odd decimal digit** to every prefix $h=2^{m+3}+1$, $m\ge1$. The five output exteriors need not be adjacent to each other or to the parent by one arc.

Use a fixed tournament $A$ with marked set $I$, a new singleton $b$, and a lane $Y$ containing $X$ preceded by the indicated number of ordinal singleton positions. Orient $b\to Y\to A$ and $I\to b\to A\setminus I$. The row masks and exact profile coefficients are:

| digit $d$ | $A$ rows | $I$ | extra pads in $Y$ | $(G,\beta)$ |
|---:|---|---:|---:|---|
| 1 | 2,12,25,17,3 | 16 | 2 | (20,11) |
| 3 | 2,12,25,17,35,15 | 34 | 0 | (80,13) |
| 5 | 2,4,41,51,39,3 | 32 | 1 | (40,15) |
| 7 | 2,4,9,19,39,15 | 32 | 2 | (20,17) |
| 9 | 2,12,41,17,39,11 | 12 | 0 | (80,19) |

Here $G$ counts ordered two-path partitions $(R,S)$ of $A$ with $R$ nonempty and ending in $I$, and $B$ counts the subset where $S$ is empty or starts outside $I$. The intercept is $\beta=H(A)+B-G$. Every Hamilton path has at most two runs in $Y$: any run other than an initial one must follow $b$. The two-run skeleton has $G$ options; the one-run skeleton has $H(A)+G+B$. Consequently the shell has height $GQ_2(Y)+\beta H(Y)$; the extra pads multiply $Q$ by their power of two. The [independent replay](decimal-tournament-check.py) verifies every shell profile by subset path DP in both directions, verifies the strong quotient with singleton $Y$, and verifies three fully inflated examples for each digit. Replacing that singleton by any nonempty tournament preserves strength by lifting quotient paths and cycles.

**Inverse and renewal boundary.** For a prescribed ordinary count $N$, take $h=\lfloor N/10\rfloor$ and $d=N\bmod10$. This literal decimal shelf accepts exactly when $d$ is one of the five odd digits and $h-1$ is a power of two at least 16; the exponent fixes the exact filler width. No child re-enters the same prefix family: $h_d-1=10\cdot2^k+9+d$ lies strictly between $2^{k+3}$ and $2^{k+4}$ for $k\ge4$. For a general filler the correction is $dD$, not the absolute digit $d$. This is one exact step on an infinite family, not a repeatable universal decimal grammar.

There is also an absolute inverse on the larger dyadic filler language from the two-read packet. With $D=19^a217^b223^c599^j$ and $Q=2^m$ (the corresponding width budget must hold), $Q\ge2D$ places the child between $80\cdot2^m$ and $89.5\cdot2^m$; its bit length fixes $m$. Subtract $80\cdot2^m$ and test each divisor $10+d$, then factor by repeated division by those four fixed bases and read the singleton width. The digit is uniquely determined if a representation exists: the coefficients 11,13,15,17 each contribute a distinct prime support outside the four atom bases, and coefficient 19 contributes none. This is a constructive absolute selector for its accepted targets, not a claim that arbitrary targets pass.

## Minimal literal odd neighbours of 6174; 495 at two adjacent orders

The even number 6174 is not itself a tournament Hamilton height. Its two adjoining odd values do occur on one marked order-ten arc:
$$6173\leftrightarrow 6175.$$
The lower outgoing masks are $(720,733,153,129,776,31,188,304,623,236)$, and the sole flipped edge is $(2,6)$ (zero-based). The two outputs are strong and are independently checked by subset path DP.

This order is minimal even for achieving either endpoint alone. The [anchored order-eight scan](kaprekar-small-scan.cpp) verifies that 661 is the maximum across all $2^{21}$ possible complements of a fixed Hamilton frame. Every tournament can be relabelled into that frame. At order nine the number of Hamilton paths is at most the sum of the nine eight-vertex deletion heights, since each Hamilton path specifies its first vertex; that is at most $9\cdot661=5949$. Smaller orders cannot do better because adding an ordinal singleton preserves height.

At the other Kaprekar value, the exact lower seam $493\leftrightarrow495$ is realized strongly on nine vertices by $(148,181,176,391,136,89,415,32,183)$ and edge $(0,5)$. The upper seam $495\leftrightarrow497$ is realized strongly on eight vertices by $(194,180,9,19,37,77,158,60)$ and edge $(1,6)$. The lower seam does not exist at order eight: the same anchored exact scan covers all 10890 anchor presentations at 495 and reverses **every** edge, including anchor edges. Its answer is 0 for target 493 and 8910 anchored marked presentations for target 497. Smaller lower seams would pad to order eight, so the new order-nine lower seam is minimal. This is not a claim that the two seams have a common 495 center.

    python3 decimal-tournament-check.py
    g++ -O3 kaprekar-small-scan.cpp -o kaprekar-small-scan
    ./kaprekar-small-scan 495 493
    ./kaprekar-small-scan 495 497

## Uniform-shape obstruction to a universal append

A fixed exterior treating the whole input tournament \(M\) as one uniform module cannot implement \(H\mapsto aH+b\) for all nonempty fillers with a nonzero intercept. If its matching depth is positive, test arbitrarily large transitive fillers: \(H(M)=1\) but the exterior count grows at least as \(1+Q_2(M)=1+2^{|M|}\), by the [reentry lower bound](../chanpak_fork_curvature_20261123/UNIFORM-AMPLIFIER.md). If its matching depth is zero, the extension is ordinal and its count is precisely \(cH(M)\), where the multiplier \(c\) is itself an odd tournament count. Compare fillers of heights 1 and 3 to force \(b=0,a=c\). Thus even the slope 10 cannot arise as a universal exact scalar append from one fixed uniform context. The positive shelf above avoids this obstruction by attaching to the same two-read filler with distinct exteriors, not by using the parent as an unrestricted single module.
