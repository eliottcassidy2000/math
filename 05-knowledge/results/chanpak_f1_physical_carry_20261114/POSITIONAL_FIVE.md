# Physical base-five digits with a single final absolute check

A background follow-up found a complete five-symbol digit on a
five-path puncture. We independently checked its literal tournament
profiles, extended it to **every** positive left address, and proved
a minimum block width and a quantitative obstruction to holding the
left source fixed. This is an exact positional decoder, not a proof
that its final quotient always takes the required value.

## Finite, labelled digits

Take the four-vertex core of count five with out-masks (8,1,3,6).
As in the [smallest switch](MINIMAL_SWITCH.md), surround it by two
transitive pieces of lengths a,b and add an apex. The physical
contacts are the binary word of value x0 at declared width a, an
out-mask U inside the core, and the binary word of value z0 at
declared width b. Set D=a+b. These choices have h=5 and f=32·2^D.
For the original alphabet the five rows are listed below; Q and K
are outputs, not assigned weights:

| a | b | x0 | z0 | U | Q | K |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 6 | 1 | 32 | 12 | 4000 | 3169 |
| 7 | 0 | 101 | 0 | 0 | 3676 | 2917 |
| 6 | 1 | 50 | 0 | 12 | 3292 | 2613 |
| 2 | 5 | 3 | 0 | 1 | 4048 | 3205 |
| 4 | 3 | 12 | 0 | 7 | 3724 | 2951 |

They all share (h,f,P,S,T)=(5,4096,3232,0,0). Both Q and K
separately take all residues modulo five. The other four explicit
alphabets are in the replay table. They let the same decoder work
at any positive left-wing address X: choose its residue r modulo
five, and the corresponding alphabet. In each row of this table
P,S,T are shared, S·T=0, and the five Q residues are distinct.

| r | f | P | S | T |
|---:|---:|---:|---:|---:|
| 0 | 4096 | 3232 | 0 | 0 |
| 1 | 4096 | 3491 | 5 | 0 |
| 2 | 8192 | 1438 | 0 | 5 |
| 3 | 4096 | 600 | 0 | 5 |
| 4 | 4096 | 1438 | 0 | 5 |

In the appropriate bank the five residues K+rQ are also distinct,
and f r+P is nonzero modulo five. All 25 choices are actual masks
checked by a literal subset-path count. The field-with-one-element
lens here has a precise role: a place has five physically named
choices, rather than a formal polynomial coefficient promising
five choices. The ordered geometry retains which choice was used.

## The absolute decoder with independently movable wings

Write g for the usual transitive-spine one-apex count and p(t)=t mod 2.
Use a left wing of positive binary address X, then k≥1 alphabet
blocks in ordinal order sharing one apex, and finally a right wing
of positive reverse-complement address Y. The first vertex is
beaten by the apex and the last beats it, so the full graph is strong.

Set H_0=g(X) and P_0=X. For digit d_j at position j, 1≤j≤k,
the common-apex formulas give
\[
\begin{split}
 P_j&=fP_{j-1}+P5^{j-1},\\
 H_j&=5H_{j-1}+5^{j-1}(K_{d_j}-5)+P_{j-1}Q_{d_j}
                       +{\bf1}_{j=1}p(X)T .
\end{split}
\]
The adjacent-boundary correction after the first digit vanishes
because the two signed endpoint counts have product S·T=0.
The total after the right wing is
\[
 N=H_k+5^k(g(Y)-1)+P_kY+5^{k-1}S\,p(Y).
\]
In particular the right correction is independent of the digit
word, even though it depends on the **actual** reader Y.

Subtract that correction from the prescribed N, then peel from
j=k down to 1. For j≥2, P_{j-1} is nonzero modulo five, because
P_1≡f r+P and later factors multiply by the unit f. Thus Q across
the five choices gives the unique residue-matching digit. At j=1
the choices K+rQ give the five residues; the first-boundary
correction is shared. Subtract the exact increment and divide
by five. The target is attained for (X,Y,k) **if and only if** the
remaining integer equals g(X). The procedure uses O(k) carry
steps (on integers of growing bit length), and returning a positive
witness gives its actual contacts. In particular there is exactly
one word per residue class modulo 5^k, for every fixed X,Y,k.

The original X=10,Y=1 chart gives, at k=1, 87,401; 83,909;
79,765; 87,917; 84,423. Literal 17-vertex counts verify the list.
Its blocks have eleven vertices, so its order is 11k+6.

## Optimal size in this model, and a fixed-source obstruction

A count-five puncture has just one nonsingleton component, the
[unique strong four-vertex count-five core](../chanpak_movable_reader_returns_20261114/SMALL_ENDPOINT_BANK.md).
Every possible contact pattern is therefore captured by the
(a,b,x0,z0,U) enumeration at a fixed f=32·2^D. For each D group
the masks by their common (P,S,T), and record the actual Q residues.
The complete finite test gives:

| D | exposed groups | groups containing all five Q residues |
|---:|---:|---:|
| 0 | 16 | 0 |
| 1 | 62 | 0 |
| 2 | 176 | 0 |
| 3 | 416 | 0 |
| 4 | 896 | 0 |
| 5 | 1856 | 0 |
| 6 | 3784 | 0 |
| 7 | 7624 | 40 |

All graph orientations in this finite test are covered by the
unique-core fact and the ordinal-position split. Thus eleven
vertices is minimum for a count-five block giving a full Q
alphabet at fixed (h,f,P,S,T). This is a limitation on this
particular stationary interface, not on other positional schemes.

For fixed X,Y the outputs at width k lie between positive constants
times f^k (lower bound P_kY≥X Y f^k; upper bound by the geometric
sum in the recurrence). Hence the number of distinct fixed-wing
outputs up to B grows as Θ(B^(log 5 / log f)), not linearly.
Even allowing **every** positive Y cannot cover all odd values
if X is held fixed. At width k, there are at most B/(X f^k)
possible positive Y and at most 5^k words per Y. Summing gives
at most
\[
                         \frac{5B}{X(f-5)}
\]
distinct values up to B. This is a rigorous absolute obstruction
even though the residues at every width are complete. A spacing
proof based solely on increasing the right reader cannot close
the remaining quotient check. Changing the left address as well
is allowed by the decoder, but no proof guarantees a suitable
pair of wings at every prescribed target.

Run **python3 replay_positional.py** for the 25 physical digit profiles,
literal controls, the inverse checks, and the complete minimum-width test.
