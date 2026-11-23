# Two physical tripling families from a four-vertex core

**PROVED.** For every \(r\ge2\) and \(w\in\{3,5\}\), let \(h=11w\,2^{r-1}+7\). There are two fixed, distinct physical arc reversals giving
\[
                    h\longrightarrow 3h\longrightarrow 3h+2 .
\]
Their complete address domain in the following construction consists of these two rays. Two particular addresses also admit a third reversal from \(3h\) to \(3h-2\). They give the forks \((117,351,349,353)\) and \((139,417,415,419)\). The coefficient 11 is the sum of twice the source's linear address coefficient and its ordinary height. It follows from the path counts of this core.

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

The [unrestricted center census](PORTS.md) gives the ordinary ten-vertex isomorphism classes. From this directory run `python3 -B eleven_rays.py` to reconstruct the readings and the two literal ten-vertex witnesses.
