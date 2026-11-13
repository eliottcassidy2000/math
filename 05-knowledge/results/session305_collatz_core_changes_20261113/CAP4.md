# Four-valuation compiler by an exact cylinder partition

**Proved by a finite certificate:** Every positive odd $n$ outside
`{1,5,7,21}` admits a strictly decreasing block of at most 26 mixed-sign
odd Collatz edges, each with exact valuation in `{1,2,3,4}`. Every
odd output and every count in the reverse dyadic ladders avoids 7
and 21. This sharpens the cap 5 in [the elementary argument](ARITHMETIC.md)
at a cost in block length. Whether cap 3 suffices for uniform
bounded-length descent remains open here.

The [certificate](cap4_certificate.txt) has 1068 pairwise disjoint
binary cylinders, moduli at most $2^{43}$. Each line contains a
hexadecimal residue $r$, a decimal precision $k$, and a word.
Digits 0–3 encode a minus sign and exact valuation 1–4;
digits 4–7 encode plus with valuation 1–4. A cylinder is
$n\equiv r\pmod{2^k}$. The replay neither trusts nor invokes the
search that produced the words.

## Why this finite proof applies to unbounded integers

A word acts on its cylinder by
$$
 x_j=(p_j n+b_j)/2^{A_j},\quad
 (p_0,b_0,A_0)=(1,0,0),\quad
 (p_{j+1},b_{j+1},A_{j+1})
  =(3p_j,3b_j+s_j2^{A_j},A_j+a_j).
$$
At each prefix, the verifier checks $A_j<k$ and
$p_jr+b_j\equiv 2^{A_j}\pmod{2^{A_j+1}}$. Thus the value is
odd for *every* member of the cylinder, not only the representative;
the stated valuation is exact by comparing consecutive cumulative
exponents. It also checks $p_j>0$ by construction, $p_\ell<2^{A_\ell}$,
and the following maxima over all relevant rows or prefixes:
$$
 \max \frac{b_\ell}{2^{A_\ell}-p_\ell}
  = \frac{144967391}{5077565} <29,\qquad
 \max \frac{12\cdot 2^{A_j}-b_j}{p_j}
  = \frac{102660143}{1594323} <65.
$$
Consequently, for $n\ge65$, the output is less than $n$ and every
odd endpoint in the word is larger than 12. By the bad-ladder
classification in [the elementary argument](ARITHMETIC.md), the only
possible reverse-ladder obstructions involve an endpoint at most 11.
The other forbidden endpoint, 21, cannot be an odd Collatz output
since $3n+s\not\equiv0\pmod3$.

Pairwise disjointness follows by checking that no shallower cylinder
is a prefix of any deeper cylinder. The sum of cylinder measures is
exactly $1/2$. Since their moduli all divide $2^{43}$ and every
residue is odd, they cover exactly all odd residues modulo $2^{43}$.
Finally, the verifier replays explicit short words on every admissible
odd $n<65$; none is omitted, and the longest such word uses 10 steps.
This is symbolic certification of the infinite region plus an exact
finite boundary, not testing to a cutoff.

## Comparison and next use

Both this certificate and the elementary compiler use only the
arithmetic seeds 1 and 5 and the physical operations one would
like to lift: $D_+,D_-$ and a certified peel at the end of each
ladder. The four-valuation proof allows a smaller inventory of
dyadic ladder depths, but longer odd-edge blocks. It neither
provides tournament lifts nor addresses convergence for either
single sign.

Replay: `python verify_cap4.py` or `python -O verify_cap4.py`.
