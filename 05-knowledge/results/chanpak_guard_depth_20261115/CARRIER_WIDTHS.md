# Which fixed-bit carrier windows can have longer companions?

The [published small-core carrier table](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_waves_20261113/scripts/carrier_forms_orders3_5.txt)
has 168 physical forms, each
[
 F(X,Y)=qXY+tX+sY+a[p(X)+p(Y)]+f[X\text{ odd}]+e[Y\text{ odd}]+n
]
on positive integer addresses, with
(p(2x+1)=p(x),\ p(2x)=p(x)+x+[x\text{ odd}]).
The simultaneous absolute placement conjecture for the whole bank is
[open](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/hypotheses/HYP-9722-finite-carrier-conjecture-every-odd-n-from-333-is-a-value-of-one-of-168-digit-quadratic-carrier-forms.md).
This note studies a narrower matter: whole dyadic cylinders with
matched nonconstant terms. The copied table has Git blob
`f58b1e2679edd0740a12247acfef423af1170b43`.

## Exact positive window and controlled limits

The independent background derivation found three consecutive physical
odds with both addresses free and only four fixed bits per lane. If
(F_{c,m}) denotes order-five class (c), mask (m), put
[
 D=1344xy+427x+245y+9[p(x)+p(y)+(x\bmod2)+(y\bmod2)].
]
For (x,y\ge1),
[
 F_{4,4}(4y+1,4x)=D+89,\quad
 F_{4,10}(8x+1,2y)=D+91,\quad F_{4,17}(2x,8y+2)=D+93.
]
The first triple is 2123, 2125, 2127. Also (F_{5,3}) and
(F_{5,5}) differ by exactly two with **no** fixed bits.
These are changing-core windows, not a claim that a desired target is
inside one of them.

The direct coefficient scan checks every cylinder with at most nine
fixed input bits in total (allowing transposed addresses). The maximum
run is three, first at four bits; no four-point run is found. This is
exact for the stated menu: for any (0\le u<2^d) expansion gives
[
 p(2^d x+u)=p(x)+(2^d-1-\operatorname{rev}_d u)x
    +[d>0,\ u<2^{d-1}](x\bmod2)+p(u).
]
This follows by appending one low bit to an existing suffix. The
coefficient signature contains the two linear, two parity, mixed,
and common surcharge coefficients; only the constant is omitted.
These terms determine equality up to a constant: take mixed
differences to read the mixed coefficient, then vary one coordinate.
The surcharge is not linear plus parity and constant, since it is zero
at all (2^j-1) but nonzero at (2^j) for (j\ge1).

## Transfer bits between the two inputs, rather than enumerating all words

Two different base forms with the same signature must have the same
(a), and the same (q) and total suffix width in this table: at fixed
(a) no pair of distinct (q)'s differ by powers of two. Consider the
shift from widths ((d,D)) to ((d+j,D-j)), where the first subscript
labels the shorter-X form. Equality of the X and Y slopes, reduced
modulo the shorter respective powers of 2, forces
[
 u_2=2^j u_1+U,\qquad v_1=2^jv_2+V.
]
Here (u_i,v_i) are the X and Y suffixes of form (i) and (U,V)
have width (j). The coefficient (a) is odd, so after bit reversal
there is no ambiguity in these reductions. Removing the common suffix
parts, the necessary and sufficient equations for the two **linear**
slopes to agree are
[
 \begin{split}
 qV+a\operatorname{rev}_j U &=2^j(t_2+a)-(t_1+a),\\
 qU+a\operatorname{rev}_j V &=2^j(s_1+a)-(s_2+a).
 \end{split}
]
There is no dependency here on the bits inherited by both forms.
Their remaining disagreement is a constant plus two parity terms;
each can be frozen by spending one more bit on the corresponding
address. If the original shared suffix is nonempty, those bits have
already been frozen. Same-width slope matches force the same two
suffixes; changing widths is what enabled the displayed triple.

The equations have a finite mirrored carry decoder at **every** width.
Read the second equation from low to high and the first from high to
low. At step (i), choose bits (u=U_i,v=V_{j-1-i}).
Writing (K_1,K_2) for the right sides, the carries update by
[
 c'_2=(c_2+q u+a v-(K_2)_i)/2,\qquad
 c'_1=(K_1)_{j-1-i}+2c_1-q v-a u.
]
Start ((c_2,c_1)=(0,\lfloor K_1/2^j\rfloor)) and finish at
((\lfloor K_2/2^j\rfloor,0)); all carries are integers in
(0,\ldots,q+a-1). Reject a negative right side outright.
These rules are just the ordinary bit equations solved in opposite
directions, so they both construct and certify the suffixes.

For each base pair the *accepted widths for slope matches alone* are
ultimately periodic. For large (j), the low end of each (K_i) is
a fixed-width two's-complement tail of a constant and the remaining low
(j) bits are all ones. The high carry is constant. Thus a fixed entry
transfer and a fixed exit transfer enclose a power of a single finite
Boolean transition matrix. The powers are eventually periodic by
pigeonhole in the finite matrix semigroup. This assertion does not
include the exact offset; the checksum can still grow across the
pumped center.

## Completeness of the bounded span test

If several lanes agree in their nonconstant signatures, choose
extremes of the X width and call their span (j). Strip off the same
outer suffixes, leaving the two extremal encodings
((x,2^j y+V)) and ((2^j x+U,y)).
At intermediate width (i) the only candidate suffixes are the
high (i) bits of (U) and the high (j-i) bits of (V), by the same
reduction of slopes. Test every base form at these inputs, and each
parity choice of the common (x,y). This includes all window offsets
of the original matched cylinders, even when their total width was
much larger than (j+2). Conversely every returned offset group becomes
physical after fixing the two parity bits. The span test through 16
finds no four-point run either. This is **not** an all-width exclusion,
not a restriction on unrelated carrier identities, and not a result
about absolute placement.

The next exact bottleneck is to attach the offset evolution and the
absolute baseline to the mirrored carry state. A successful pump must
sell new baselines, not just new ways to match their slopes. This is
the same distinction the [guard clock](GUARD_DEPTH.md) makes in a
different physical class.

Reproduce with `python scan_carrier_windows.py 9` and
`python scan_carrier_transports.py 16`; compare the saved receipts.
