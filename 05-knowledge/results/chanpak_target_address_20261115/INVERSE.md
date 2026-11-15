# A target-first inverse for an exact two-address repair

**Proved:** a complete constructive inverse in
(O(\sqrt N\log^2 N)) arithmetic operations for the literal six-vertex
center family described here. The family has a three-count repair and
allows its left and right baselines to be chosen independently. This
is not an inversion of the unrestricted spectrum.

## Physical parameters

Use an ordinal join L⇒C⇒R and one apex, which beats each marked
vertex. Orientations in C use a bit for each lexicographic (i,j), i<j,
where 1 means i→j. Use three marked centers (code, mask):
(8850,54), (8881,54), (8881,53), with 1-bits of the mask marked.
Each has the same tuple

    (h,q,s,e,t,f) = (45,432,258,30,129,15),

and its own apex count is n_C=135,137,139, respectively.
Here q counts optional two-slot covers; s,e are the total and signed
odd-minus-even partial path counts starting in the marks, weighted
by Hamiltonian counts of their complements. t,f are the dual totals
ending in the complement of the marks. These are exact literal
packets: verify_inverse.py rebuilds the orientations and recomputes them.
A common tuple makes the three-choice repair valid in *any* fixed
ordinal surroundings (with spacing multiplied by H(L)H(R)).

For transitive L and R, encode L's marks as the ordinary binary digits
of x≥1, and R's marks as the reversed complement of the ordinary
binary digits of y≥1. The inherited
[exact cut transfer](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_selective_target_20261114/CUT_STATE.md)
then gives

    N = 432xy + 129x + 258y + 15(x mod 2) + 30(y mod 2)
                 + 45(p(x)+p(y)) + k,      k in {135,137,139}.

The mixed term comes from a crossing apex cycle with optional
traversal of C; both signs on either end of such a crossing cancel
over the nonempty center. The neighboring seams supply t·x,f·parity(x)
and s·y,e·parity(y). Each exterior alone has apex count 1+p;
its increment scales by h. The unsigned binary checksum is

    p(0)=0, p(2m)=p(m)+m+(m mod 2), p(2m+1)=p(m).

This follows from the weighted pair formula for one apex (ending in
0 creates all new 1→0 inversions; ending in 1 creates none).
Induction gives 0≤p(m)≤m, with strict upper bound m−1 for odd m.
All coefficients except k are divisible by 3, so each prescribed
N selects exactly one of the three center choices before any search.
No symmetry condition between x and y is necessary.
The resulting tournament is strong: the apex reaches the first
site of L, and the last site of R beats the apex. Every site of
the ordinal chain can be reached from that first site and can reach
the last.

## Peel low bits, not candidate graphs

The residual equation on fixing either address and k is
`F(a,d;y)=a y+45p(y)+d(y mod 2)=R`. Its digit recurrences are

    F(a,d;2m)   = F(2a+45,45;m),
    F(a,d;2m+1) = a+d+F(2a,0;m).

At any node with coefficient u and parity coefficient v, any answer
m obeys

    ceil((R−v)/(u+45)) ≤ m ≤ floor(R/u),   m≥0.

A negative residual is impossible, and residuals never increase.
Independently of the chosen low bits, u grows by at least 2^depth.
Once u≥2 sqrt(45N)+2·45+2v+4, the real interval has width less
than one, by
`R/u−(R−v)/(u+45)≤45N/u²+v/u`.
The search then tests at most one m and recovers its literal bit
address by reinstating the peeled bits. One need never use any
monotonicity of p. Without pruning, the branch count is already
(O(1+\sqrt N/a)); computing/checking a leaf costs O(log N)
operations in integers of O(log N) bits.

As every coefficient is nonnegative, xy≤N/432. Hence at least one
address is ≤sqrt(N/432). Enumerate that side and solve its rows by
the digit inverse, then do the same with the sides interchanged.
The slopes are 432x+258 or 432y+129. Summing their inverse slopes
gives a harmonic sum: the total cost is (O(\sqrt N\log^2 N))
arithmetic operations, including duplicate elimination and literal
reconstruction. This returns every address/center witness in the
bank and certifies non-membership in the bank when there are none.

Examples: x=y=1 gives 999,1001,1003. Target 35,099 selects k=137
and the exact inverse returns x=5,y=14. The resulting literal
tournament has 14 vertices and 35,099 Hamiltonian paths. The same
decoder rejects 36,833, which the previous marked-host bank can
reach. No universality follows from either a growing candidate
count or a complete inverse.

## Replay

    python verify_inverse.py
    python inverse.py 999 1001 1003 35099 36833

The first command independently counts 15 full literal tournaments,
checks 36,000 scalar inverses and compares the complete decoder
against direct pair enumeration for every odd target through 12,345.
