# Growing the middle without losing the signed switch

An alternate way to keep an exact +2 is to **change the center itself**,
while keeping two ordinary flank addresses independent. This gives a
four-state all-depth grammar with an exact spacing obstruction. It is
not a finite list of precomputed centers.

Use the signed profile (h,P,L,R,l,r,K) of the
[free seam](../../chanpak_wheel_square_selector_20261117/FREE_SEAM.md);
l and r are the signed companions of L and R. For ordinal composition
of marked blocks X→Y, direct partition of the selected path gives

    h = h_X h_Y,              P = P_X P_Y,
    L = L_X P_Y + h_X L_Y,    R = R_X h_Y + P_X R_Y,
    l = h_X l_Y,              r = r_X h_Y,
    K = K_X h_Y + h_X K_Y - h_X h_Y + L_X R_Y + l_X r_Y.

For example, the contribution when a selected path starts in X
picks an arbitrary ordered optional cover of Y. Its signed sum is
zero for nonempty Y by the following involution: compare the first
vertices of the two paths and move the one that beats the other to
the front of the other path (if one path is empty, move to it). On
the next application that same vertex beats the new first vertex
of its original path, or the original path is now empty, so it moves
back. The sign changes with the first slot's parity. A selected
path confined to Y instead gives h_X l_Y.
The last line is the adjacent marked-interval law.

The old deletion core already splits as X→Y with
X=TT3 marked on its first vertex, and Y=(8,1,3,6) marked mask5,
of profiles (1,8,4,6,0,0,5) and (5,32,14,14,-2,-2,15).
Insert any nonempty `TT_n` in the middle of X→Y, with marked
contact word U in n bits, and write V for its reverse complement.
At transitive size n, h=1, P=2^n, L=U, R=V, K=F(U), l is minus
the final bit, and r is minus the complement of the first bit.
The ordinary readings L=R of the new whole center are balanced iff

    4U - 5V = -2*2^n + 2.

The signed values stay (-2,0) for every such word. The whole
center has h=5, P=256*2^n and common L=128*2^n+32U+14. Its K is

    5 F(U) + 30U + 64*2^n + 22 + 2[U odd].

This follows from the composition formulas and substitution for
5V. Both outside flanks of the free seam may be varied exactly
as before. The full host remains strong via its positive flanks.

A structural constraint on this approach is worth making explicit.
Between fixed nonempty terminal blocks, the signed difference of
any ordinal center is multiplied by h of every interior block.
Thus each interior block must have h=1 for an exact unit switch.
A nontransitive tournament contains a directed triangle; a triangle
has three Hamiltonian paths, and insertion of the remaining
vertices preserves at least three paths by the usual path-extension
injection (for each old path choose a legal insertion, then delete
that vertex to recover the path). So h=1 means transitive.
A large nontransitive correction would have to be put at a terminal
block or use an architecture beyond this simple ordinal center.

## All depths by paired radix-8 carries

At length n=3k the balance equation becomes

    4U + 5 rev(U) = 3(8^k - 1)

with bit reversal in 3k bits. Write `rho` for 3-bit digit
reversal. If a,b are the next low and high octal digits, and
l,h the incoming low and high carries, respectively, the new
carries are ((4a+5rho(b)+l-21)/8, 21+8h-4b-5rho(a)).
The equivalent redundant RHS has digit 21 at each position.
Beginning with zero carries at both ends, all truthful carries
lie between -3 and 6 (local sum between 0 and 63 and induction).
The reachable paired states, with labels (high,low), are exactly:

    A=(0,0)    : (4,0) -> B, (4,4) -> A ; center digit 4
    B=(-2,5)   : (7,3) -> C, (7,7) -> D ; center digit 7
    C=(3,3)    : (3,3) -> C, (3,7) -> D ; center digit 3
    D=(5,-2)   : (0,0) -> B, (0,4) -> A ; center digit 0

Start A; if no middle digit remains accept A,C only. With one
digit left there is exactly the stated choice at every state.
At even length each state has one successor that can terminate
and one that cannot. Thus the exact number of distinct balanced
words at length 3k, for k>=1, is 2^floor((k-1)/2).
For example U=100100100 and 100111000 give two true centers at
k=3, translating both exterior counts together by
730+640(u+w), for the same arbitrary positive odd u and even w.
The octal digits are (4,4,4) and (4,7,0); the displayed U's
here are binary.

## Why this particular growing bank cannot make a tiny adjustment

In each state the two outgoing edges have the same high digit;
their next high digits (or their middle digits on the final odd
step) differ by 3 in absolute value. If two words first take
different edges after j common outer pairs, the first high
difference is at octal place k-2-j. The total subsequent digits
can cancel strictly less than one unit of that place. Hence the
positive difference DeltaU is at least 2*8^(k-2-j), and in
particular at least 2*8^floor(k/2) (an even-length final pair
cannot give two legal endings).

The two addresses share a suffix of 3j bits. Applying the exact
suffix law, the nonconstant contribution of this suffix to the
F difference is nonnegative except for a possible -1 in its
parity term. The high input of the smaller word is strictly
less than 5*8^(k-j-1), since its initial octal digit is 4.
Using F(x)<=x+1, the F difference is therefore at least
-20*DeltaU-1. When u,w are fixed opposite-parity positive
addresses, the full lower-count difference is at least

    (32(u+w)-70)*DeltaU - 7 >= 26*DeltaU-7.

So the count bank for each fixed height is strictly ordered,
but growing the center in this very structured way does not fill
the remaining short window.

## Fixed-scale check and exact two-center displacement

At h=5 and P=256 all centers are ordinal products of exactly one
strong 4-vertex H5 block and three singletons. Indeed strong
components multiply their h and P, so there is just one
nontrivial component; [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
bounds its order by 4; direct 4-vertex enumeration gives the
unique H5 block with P=32. Scanning its four SCC positions and
all marked masks leaves only two free-unit profiles up to reversing
the marked center:

    (5,256,142,142,-2,0,91)
    (5,256,109,109,-3,-5,55).

Their lower members at the same unordered opposite-parity
positive addresses differ by 33(u+w+1). This is another physical
second move, with no way to tune a smaller displacement at fixed
u,w in the h=5,P=256 balanced-center class. It does not rule out
changing the h,P scale or admit any unsupported all-core no-go.
