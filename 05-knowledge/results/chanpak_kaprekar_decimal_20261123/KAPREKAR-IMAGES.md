# Exact Kaprekar gap sidecar and physical binary image

## All bases: digits, fibres, and image size

Fix a base \(B\ge2\) and width \(n\ge1\), including leading zeros.
Write \(K\) for descending-digit value minus ascending-digit
value. For descending digits \(d_0\ge\cdots\ge d_{n-1}\), let
\(c_i=d_i-d_{n-1-i}\) for \(0\le i<d=\lfloor n/2\rfloor\).
The gaps are weakly decreasing nonnegative integers at most
\(B-1\). The next output depends only on this gap partition.

If there are \(k>0\) positive gaps, the exact output digits are
\[
(c_0,\ldots,c_{k-2},c_{k-1}-1,
 (B-1)^{[n-2k]},B-1-c_{k-1},\ldots,B-1-c_1,B-c_0),
\]
where the superscript in brackets means repetitions and the
first and mirrored lists are empty if \(k=1\). Subtraction
starts a borrow in the last column; it crosses the remaining
negative differences and center zeros and is repaid at the
innermost positive difference. The zero partition gives zero.

Every valid gap partition occurs: take its entries as the top
paired digits and put zeros below them. No two partitions give
the same output: at the highest differing digit place \(B^j\),
the difference has magnitude at least \(B^j\), and the entire
lower tail has magnitude at most \(B^j-1\), because every
difference coefficient is between \(-(B-1)\) and \(B-1\).
Thus the fixed-width image includes exactly
\[
 \binom{\lfloor n/2\rfloor+B-1}{B-1}
\]
distinct numbers, including zero. Every output is divisible by
\(B-1\). Odd bases give no odd outputs; for even \(B\), a
positive output is odd exactly when \(c_0\) is odd (its last
digit is \(B-c_0\)).

There is an exact sorted-multiset fibre count. Put
\(y_i=d_{n-1-i}\), so paired top digits are \(y_i+c_i\).
The choices are \(0\le y_0\le B-1-c_0\) and, independently,
\(0\le y_{i+1}-y_i\le c_i-c_{i+1}\). These ensure both sorts
and all digit bounds. If \(n\) is odd, the center may be any
digit between the last pair. Thus the number of sorted input
multisets giving this output is exactly
\[
 (B-c_0)\prod_{i=0}^{d-2}(c_i-c_{i+1}+1)
    (c_{d-1}+1)^{n\bmod2} .
\]
Width one has \(B\) sorted inputs and only zero. Weight each
multiset by its multinomial coefficient for the unsorted
word fibre.

## Binary: idempotent projection with exact inverse

Let \(a\le b\) be the unordered counts of zeros and ones; put
\(M_j=2^j-1\). The positive binary output is
\[
 K=M_aM_b,
 \qquad\text{word }1^{a-1}0\,1^{b-a}0^{a-1}1.
\]
This follows either by subtracting the packed ones or from
the borrow formula. It has \(b\) ones and \(a\) zeros, so it is
already fixed after the first step. The full input fibre at its
width is the union of complementary Hamming spheres, counted
once if their weights are equal.

The pair is unique globally, even allowing leading zero
padding. If \(a=1\), the output is Mersenne and determines
\(b\) immediately. Otherwise \(K\ge(9/16)2^{a+b}>2^{a+b-1}\),
so \(n=a+b\) is the bit length, and
\(2^n+1-K=2^a+2^b\). Decode the latter from its two set bits
or from its single bit one place above equal exponents.
Check \(a\ge2\) and \(a+b=n\); such a product is not Mersenne.
Below \(2^L\) there are exactly \(\lfloor L^2/4\rfloor+1\)
positive images: all balances at widths \(2,\ldots,L\), plus
the Mersenne output at width \(L+1\).

### Every positive binary image except 7 and 21 is physical

For odd \(j\ge5\) a transitive spine under an apex with contact
word \(OI^a O I OI^a\), \(j=2a+3\), has height \(M_j\).
This is the symmetric Mersenne-offset rail with middle run
length one. Here O is an apex out-contact, I an in-contact;
the odd-cycle count gives \(H=2^{2a+3}-1\).
Also \(M_1=1\), \(M_2=3\), and \(M_6=63\) is a literal
strong seed in the replay. For even exponents other than
six, factor \(M_j=M_{j/2}(2^{j/2}+1)\) and induct; the plus
factor is the height of a nearly transitive tournament
(a forward order with just its last-to-first arc reversed).
The only missing Mersenne is the known hole \(M_3=7\).

Use ordinal joins to multiply two attained Mersennes. If one
exponent is three, the remaining exponent \(b\) gives:
\(b=1,2\) the established holes 7 and 21; \(b=3\) the
seed 49; \(b=4\) the product \(3\cdot35\); \(b=5\) the
seed 217; and every \(b\ge6\) a signed shell. In the shell
take rows \((8,1,3,166,15,23,63,119)\), tag mask 80, a
singleton \(z\) and transitive filler \(X\) of size \(s=b-5\).
Use \(z\to X\to A\), and \(I\to z\to A\setminus I\).
The exact shell functional is \(GQ+(H(A)-L)D\), where
\(G\) counts ordered splits \((R,S)\) of \(A\) into directed
paths with \(R\) nonempty ending in \(I\), and \(L\) counts
the same splits with \(S\) nonempty starting in \(I\).
Every two-filler-run path has skeleton \(X_1,R,z,X_2,S\),
and the one-run coefficient is \(H(A)+2G-L\).
The replay certifies \((G,H(A)-L)=(224,-7)\), giving
\(224\cdot2^s-7=7M_b\). This establishes exact coverage of the
whole binary image other than the two holes, not of its
preimages or of arbitrary tournament heights.

For the perfect target \(E_p=2^{p-1}M_p\), use width
\(2p-1\) and \(p\ge3\). Direct Hamming weights give
\[
 K(E_p-1)=M_{2p-2},\quad K(E_p)=M_{p-1}M_p,
 \quad K(E_p+1)=M_{p-2}M_{p+1}.
\]
All three images at \(p=61\) are physically attainable,
but center minus upper image is \(3\cdot2^{p-2}\);
the operation does not transport a unit seam. At \(p=3\)
the center itself maps to the hole 21.

## Width two in an even base: exact cycle lengths

An input gap \(a>0\) yields the two digits \((a-1,B-a)\),
so the next gap is \(f(a)=|2a-(B+1)|\), always nonzero odd
for even \(B\). Odd representatives \(1,3,\ldots,B-1\)
represent each nonzero sign class modulo \(m=B+1\) exactly
once. On them \(f\) is conjugate to multiplication by two
on the sign classes. It is consequently a permutation.
The exact period on the orbit reached from initial gap \(a\)
is
\[
 \min\{r>0: 2^r\equiv\pm1\pmod{m/\gcd(a,m)}\}.
\]
This follows by cancelling the common divisor in
\(2^ra\equiv\pm a\pmod m\). The gap reaches a cycle after
at most one step; the numerical word reaches a cycle after
at most two. In decimal, \(m=11\) yields the single five-cycle
of positive odd image values \(9,81,63,27,45\), which are
actual tournament heights by 63 and ordinal products of 3
and 5. The formula applies to each of the five lucky-prime
bases \(B=2,4,10,16,40\), without an implication for arbitrary
widths in those bases.

Replay: `python3 kaprekar-sidecar-check.py` (standard library).

## How many trailing digits are actually enough?

The signed radix formula is \(K=\sum_{j=0}^{d-1}(B^{n-1-j}-B^j)c_j\). As \(n-1-j\ge d\), reduction modulo \(B^d\) gives exactly \(-\sum_j B^j c_j\). Read the ordinary base-\(B\) digits of the negated residue to recover the entire gap vector and hence the numerical image. Therefore on the image **last \(\lfloor n/2\rfloor\) digits are a sufficient and sharp suffix length** for injective coding, for \(n\ge2\). To see sharpness, compare gaps all ones to the same with last coordinate zero. They collide modulo \(B^{d-1}\) but give distinct outputs. The suffix code applies only after making one Kaprekar step; for arbitrary starting numbers, equal ending digits need not give the same next result.

In decimal, the units digit records \(c_0\), not the old units digit. It gives the full nonzero width-three state \(K=99c_0\) (including 495), but width four needs two digits: \(K=999c_0+90c_1\). For example 6084 and 6174 are both in the width-four image and both end in 4, but their next states are 8172 and 6174. The repo's [padded census](../kaprekar_cuboid_20260928.out) confirms the nonrepdigit width-three and width-four basins; width five already allows other cycles. This is consistent with [Prichett, Ludington and Lapenta's fixed-point analysis](https://www.fq.math.ca/Scanned/19-1/prichett.pdf). The exact quotient is not a universal convergence claim for longer lengths.

For the tournament connection, [the exact image replay](kaprekar-image.py) builds every odd numerical image at decimal widths \(2,3,4,5\), respectively **5, 5, 30, 30** distinct values. It enumerates gap vectors, synthesizes a small apex or an ordinal product, and independently verifies every count by subset endpoint DP. The even images cannot occur as ordinary Hamilton counts. In particular \(495\) and \(3087=6174/2\) occur, whereas 6174 itself is not an odd count. The [actual single-arc seam](DECIMAL-SHELF.md) at 6173 and 6175 is separate from Kaprekar iteration.

## An operational return to all five Heegner primes

With width two in even base \(B\), the nonzero next-gap map \(a\mapsto |2a-(B+1)|\) sends \(a\) to the odd representative of twice its class up to sign. At \(B=p-1\) for the odd lucky primes, this is precisely the *same doubling clock* on \(\mathbb F_p^\times/\{\pm1\}\) that appeared in the [nine-vertex Heegner registers](../chanpak_heegner_absolute_targets_20261123/NINE-REGISTERS.md). Thus the five bases \(2,4,10,16,40\) supply the cycle types \(1;\ 2;\ 5;\ 4+4;\ 10+10\). In base ten the orbit of odd heights is \(09\to81\to63\to27\to45\to09\). The replay verifies every state in all five bases as an actual Hamilton count, with the correct integer \(K(a)=(B-1)a\), not a residue proxy.

These exact height witnesses do **not** assert a cycle of single arcs on one labelled host. The labelled single-arc metagraph is a cube: every edit flips the product of the edge orientation bits; every closed walk has even length. In the unlabelled space an odd closed path of single-arc edits would have to return with an odd vertex relabelling. Indeed orientation-bit sign changes by the sign of the relabelling, and an automorphism of a tournament is even because its automorphism group has odd order (an involution would reverse a paired edge). This is a precise sidecar that a putative five-cycle lift would have to keep.
