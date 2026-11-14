# The rank-13 glue has an actual transverse repair of size two

The Leech comparison in [MRB_LEECH_SOURCES.md](MRB_LEECH_SOURCES.md)
asks for the coordinate that an exact reset discards, not just for
a shared determinant or group order. Here are two **physical**
prefixes with the same retained coordinate and different transverse
coordinates:
\[
 \begin{array}{c|c}
 \text{terrace gap word}&\text{two-state column}\\ \hline
 (1,0,0,3,0,0,1)&z_1=(291,758)^t\\
 (2,0,0,0,0,0,0,2)&z_2=(179,852)^t .
 \end{array}
\]
Gap widths are the same normalized positive/zero gap lengths used
in the terrace proof. Start at \((1,2^{r_0})^t\), then multiply
by \(D_{2^r}C_q\) at each successive positive gap. The two evaluations
are immediate from \(C_2=(3,3;3,4)\) and
\(C_6=(31,37;37,44)\). Their difference is \(2(56,-47)^t\),
and both have \(47a+56b=56125\). In particular they give *identical
ordinary counts after the length-13 corridor*, for every choice of
next width and every subsequent continuation. Their normalized
gap words are different, not just reversals of one another.

## A neighboring seam converts a hidden difference into exactly two

For a terminal width \(r\ge1\), the total after inserting \(q\)
zeros and then appending that positive width is
\[
                    H_i(q,r)=(1,2^r)C_qz_i .
\]
Let \(R_q=(56,-47)C_q\). Each of its two integer entries obeys
the tribonacci recurrence. In particular
\[
 \begin{array}{c|rrrrrrrr}
 q&8&9&10&11&12&13&14&15\\\hline
 R_q^{(1)}&5&-2&-1&2&-1&0&1&0\\
 R_q^{(2)}&-3&2&0&-1&1&0&0&1 .
 \end{array}
\]
Thus, without changing any terminal widths,
\[
 \begin{array}{rcl}
 H_1(10,r)&=&423939+505125\cdot2^r,\\
 H_2(10,r)&=&423941+505125\cdot2^r,\\
 H_1(14,r)&=&4851811+5780875\cdot2^r,\\
 H_2(14,r)&=&4851809+5780875\cdot2^r.
 \end{array}
\]
These give two repeatable width rails, with opposite signs of the
correction depending on the seam, and **no width search or approximate
cancellation**. The second entry of \(R_q\) is zero only for
\(q=10,13,14\): check \(0\le q\le17\), then use positivity and the
recurrence. So these are precisely the constant-in-\(r\) seams
for this physical pair, with length 13 the only one giving zero.
An independent interval-packing replay checks literal gap words
in addition to the two-state arithmetic.

This is the count analogue of carrying a transverse coordinate
through an identified seam: the genuine equality at 13 can be
turned on either side into a controlled signed difference at
other seams. It gives an absolute repair on two unbounded but
sparse rails, **not** a way to bring an arbitrary prescribed
target onto one of those rails. In particular the rails have
fixed affine baselines; arbitrary even gaps and arbitrary
primitive inputs are two different requirements. The next useful
question is how to vary the physical left pair and the outer
coefficient jointly while preserving the exact transverse difference.
