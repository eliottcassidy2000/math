# Two physical tripling families from a four-vertex core

**PROVED.** For every \(r\ge2\) and \(w\in\{3,5\}\), let \(h=11w\,2^{r-1}+7\). There are two fixed, distinct physical arc reversals giving
\[
                    h\longrightarrow 3h\longrightarrow 3h+2 .
\]
Their complete address domain in the following construction consists of these two rays. Exactly two addresses admit any third reversal from \(3h\) to \(3h-2\). They give the forks \((117,351,349,353)\) and \((139,417,415,419)\). The coefficient 11 is the sum of twice the source's linear address coefficient and its ordinary height. It follows from the path counts of this core.

Use the [right-flank construction](DIRECT-FORK.md) with core rows \((6,4,8,3)\), apex out-mask 5, left address zero and right address \(y\ge0\). The center is \(C\). Reverse \((2,3)\) to obtain the source \(P\) or \((0,3)\) to obtain the upper neighbor \(U\). For \(o=y\bmod2\), its three tagged-core readings give
\[
    P=f(y)+5y+o+3,\qquad C=5f(y)+14y+2o+15,\qquad U=3f(y)+15y+o+13 .
\]
Here \(f(0)=0, f(2v+1)=f(v), f(2v)=f(v)+v+[v\text{ odd}]\). Subtracting the three formulas shows that \(C=3P\) and \(U=C+2\) together are equivalent to \(y\) odd and \(2f(y)=y-5\). The function \(f\) is even-valued. Hence \(y\equiv1\pmod4\), and \(y=1\) fails. Write \(y=2^r w+1\), \(r\ge2\), \(w\) odd. Iteration gives \(f(y)=f(w)+(2^{r-1}-1)w+1\), so the condition reduces exactly to \(w-f(w)=3\). For \(w=2v+1\), this is \(f(v)=2v-2\). The recurrence gives \(0\le f(v)\le v\) by induction (the odd case has strict upper bound). Consequently \(v\le2\). Only \(v=1,2\) work, giving \(w=3,5\). Substitution yields \(P=(11y+3)/2\).

## Which apex-to-flank reversal gives the second exit?

**PROVED.** On these two rays an apex-to-flank reversal takes \(C\) to \(C-2\) exactly at \(y=21,25\); in both cases it toggles binary bit 1 of \(y\), numbered from zero. No further edit confined to the four-vertex core and its apex tags supplies this exit on either ray.

An apex-to-flank edit at position \(k\) replaces \(y\) by \(y\mathbin{\mathrm{xor}}2^k\). Clearing the highest bit can leave trailing sink vertices in the tournament; their removal preserves the height. At \(k=0\) the change is \(-16\bmod5\), which excludes \(-2\). At \(k=1\), writing \(v=(y-1)/4\) gives the change \(28-5(v+[v\text{ odd}])\), equal to \(-2\) precisely for \(v=5,6\).

For \(2\le k<r\), write \(q=2^{r-k-1}w\). Applying the recurrence to the bit changed from zero to one gives the change \(33\cdot 2^{k-1}-5(q+[q\text{ odd}])\). Equality to \(-2\) forces \(k\equiv1\pmod4\) by reduction modulo 5. Now \(k\ge5\); reduction modulo 4 excludes \(4\mid q\). The remaining \(q=w,2w\) have \(q+[q\text{ odd}]\le10\), too small for the equality.

For \(k\ge r\), put \(w'=w\mathbin{\mathrm{xor}}2^{k-r}\). Since the flipped bit belongs to the flank, \(0\le w'\le7\). The height change equals
\[
 (33\,2^{r-1}-5)(w'-w)+5\bigl(f(w')-f(w)+[w'\text{ odd}]-1\bigr).
\]
The first term has magnitude at least 61 and the second at most 40. Finally, each of the ten edits internal to the core or its apex tags gives an affine function of \(y\) for twice the change, using \(2f(y)=y-5\). The [replay](eleven_rays.py) constructs those functions and checks their slope and value at 13 to exclude \(-4\) at every \(y\ge13\).

## All ordinary arc reversals on both infinite families

**PROVED by finite sign certificates.** On these prescribed hosts an individual reversal from \(C\) to \(C-2\) exists if and only if \(y=21,25\). The preceding apex-to-flank arc is unique at each of the two addresses. This allows the candidate arc and its position to vary with \(r\).

The right flank has one low digit, a uniform transitive run of size \(s=r-1\ge1\), and \(\ell\) high digits from the fixed binary expansion of \(w\), where \(\ell=2,3\). Let \(\Delta\) be the change in height of the center. If both endpoints of an ordinary arc lie in the run, reversing it gives \(\Delta\ge0\). Indeed every path using its original arc \(i\to j\) maps to a distinct path by swapping the adjacent vertices. Every other vertex that beats \(i\) beats \(j\), and every vertex beaten by \(j\) is beaten by \(i\).

For an ordinary edge outside the run, keep the run as a single variable module. For an edge with one endpoint in the run, distinguish that vertex and write \(s=a+1+b\). Keep the two residual runs as separate modules. In either orientation one Hamilton path can visit each module at most three times. Before the flip, ordinary vertices outside such a module lie on its left or on its right, with all corresponding arcs directed forward. An excursion from the module back to it must traverse the apex or the single reversed ordinary arc, each available at most once. This proof also applies when the reversed arc stays on one side of a module.

For a module of size \(t\ge1\), the ordered path-cover counts for one, two and three visits are respectively \(1,2^t-2,3^t-3\cdot2^t+3\). Deleting the interiors of the visits shows that the entire height is a linear combination of these counts, separately in each variable module. With two modules the slots for both can be filled independently after the outside path pattern is fixed. Treat a module of size zero by deletion. Consequently on each positive variable \(t=m+1\) the function \(\Delta+2\) lies in the span of the three nonnegative functions
\[
        1,\qquad 2^m-1,\qquad 3^m-2\cdot2^m+1.
\]
The last counts maps to three colors that use each of two specified colors at least once. If the values at \(m=0,1,2\) are \(A,B,C'\), twice the coefficient vector in this basis is \((2A,2(B-A),C'-3B+2A)\).

There are \(d=5+\ell\) ordinary vertices outside the run. For the \(\binom d2\) outside edges use one interpolation. For each of the \(d\) outside endpoints of an edge meeting the run, use the four faces \(a=0\) or \(a\ge1\), and \(b=0\) or \(b\ge1\). Interpolate separately in every positive coordinate. The [exact replay](eleven_rays.py) computes the samples by cutting paths at the candidate arc on hosts with \(1\le s\le7\). In all \(\binom d2+4d=49,60\) charts, respectively, the constant coefficient is nonzero and every other coefficient has the same weak sign as the constant. Thus \(\Delta+2\) never vanishes for an ordinary reversal. Together with the algebraic apex calculation this exhausts the entire arc set for every width.

The [unrestricted center census](PORTS.md) gives the ordinary ten-vertex isomorphism classes. From this directory run `python3 -B eleven_rays.py` to reconstruct the readings, the two literal ten-vertex witnesses, and all 109 sign certificates.
