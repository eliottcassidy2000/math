# A physical two-path carry with the smallest possible puncture count

The previous fixed-profile unit switch used a puncture with h=15.
Here a movable four-vertex core gives h=5. This is **minimal in h**
among all fixed-six-coordinate switches, even if order is unbounded
and the two punctures need not be isomorphic. It is not asserted to
minimize the vertex order.

Take the four-vertex core A with outgoing masks (8,1,3,6). We use
two ways to build a nine-vertex puncture around A, then add an apex.
Place each transitive wing before or after A in the natural ordinal
order, with individual apex contact bits as specified. In the table
x is the ordinary binary left address, y the reverse-complement right
address, and U the apex's outgoing contact mask inside A.

| x | y | U | puncture order | augmented K |
|---:|---:|---:|---:|---:|
| 1 | 12 | 4 | 9 | 579 |
| 11 | 1 | 2 | 9 | 581 |

The whole ten-vertex outgoing-mask rows, with the apex last, are
\[
\begin{split}
 &(510,1008,994,486,1004,448,384,768,512,105),\\
 &(510,1020,504,496,896,272,816,864,512,45).
\end{split}
\]
Both alternatives have exactly the same exposed profile
\[
          (h,f,P,Q,S,T)=(5,1024,716,810,0,0),
\]
where S,T are the signed endpoint responses (both zero here).
Thus every context glued along the same apex sees a true gain of
two when this is the only nontrivial puncture block (more generally
two times the product of other puncture counts). The witnesses and
the six coordinates are independently verified by literal subset DP.

To see the moving support, first compare A followed by a three-bit
right wing with address 4, using U=4, against a three-bit left wing
with address 3 followed by A, using U=2. In these mini-punctures the
unsigned (P,Q) pairs already agree. Namely the intrinsic core
profiles (f,P,Q) for those masks are (32,9,21) and (32,6,18),
and the two unsigned equalities are
\[
\begin{split}
 8\cdot9+5(7-\operatorname{rev}_3(4))&=32\cdot3+6=102,\\
 32\cdot4+21&=8\cdot18+5(7-\operatorname{rev}_3(3))=149.
\end{split}
\]
The wing which moved positions is not a formal number: its
declared three-bit address includes its leading zero. Put a singleton before all seven old vertices which the apex beats,
and a singleton after all seven old vertices which beats the apex. The two new boundary contacts fix
S=T=0 in both punctures and make both augmented graphs strong.
The [ordinal start formula](SUPPORT_DECODER.md) gives exposed
responses 716 and 810. The local count difference can either be
computed from the [two-wing formula](../chanpak_movable_reader_returns_20261114/bilateral.md)
or checked directly.

More generally, prepend and append common nonempty transitive
wings of external addresses X,Y instead of the two singletons.
The two final charts are (X,8Y+4;U=4) and (8X+3,Y;U=2).
The bit addresses are physical for every X,Y>0. The exposed
puncture interfaces still agree (signed coordinates are fixed by
the common boundary vertices). The second augmented count minus
the first is exactly
\[
                       2(p(X)+p(Y)-1),
             \qquad p(u)=u\bmod 2 .
\]
This follows by substituting the binary suffix recurrence in the
two-wing identity, once for each fixed three-bit suffix. The X and Y
coefficients cancel by the two unsigned equalities; the four boundary
parity values are −2,0,0,2 in the even/even, even/odd, odd/even,
odd/odd order. Thus there are both exact full-profile neutral changes
and unit digit changes on arbitrarily large literal graphs.

## Minimality and what it buys

A puncture with h<5 must have h=1 or 3. The component product and
odd-cycle packing expansion force all nonsingleton components in
these two cases to be directed triangles: a strong component with
H=3 has precisely one directed triangle, and every vertex belongs
to one. The [signed support decoder](SUPPORT_DECODER.md) proves that
for every such puncture, (h,f,P,Q) determines K. In particular there
is no fixed-profile unit return, or any other switched digit, with
smaller h at any order. The four-vertex core above is the first place
where a directed contact can move while retaining the unsigned
support signatures but changing the actual count.

The new h=5 digit can be paired with the older h=17 double digit
(their f values are coprime to the product 5·17). It reduces the
integer conductor of the normalized two-generator carry to 144.
This is a way to make more physical carries at fixed usable
baselines, not a guarantee that every prescribed N enters such a
baseline window.
