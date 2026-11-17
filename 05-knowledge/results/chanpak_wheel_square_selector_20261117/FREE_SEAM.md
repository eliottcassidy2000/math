# A free reflected seam and a sub-square-root absolute selector

The [singleton reflection law](REFLECTION.md) obtains an arbitrary exact
gap by setting w-u equal to the gap. That ties up one address. To get
a fixed +2 gap *without* tying up an address, balance the unsigned
ports and leave a signed imbalance. This does occur physically.

## The marked tournament

Take the seven-vertex deletion tournament with zero-based outgoing
masks (54,32,2,119,6,20,55), apex beats-mask 26. Its full signed profile
is (H,P2,L,R,Lminus,Rminus,kappa)=(5,256,142,142,-2,0,91).
The reflected runs therefore give

    H(u,w) = 91 + 5(F(u)+F(w)-2) + 142(u+w) + 256uw + 2[w odd].

For every odd positive u and even positive w, transposing the address
roles raises the count by exactly 2. All outputs are strong: the apex
reaches the first vertex of the left flank and the last vertex of the
right flank reaches the apex, and the ordinal cross-edges go forward.
For any vertex inside the fixed middle block there is a route through
these flanks from and back to the apex, even when the middle block
itself is not strong. The two tournaments have equal order. The unsigned
boundary stays balanced, and the parity sidecar alone gives the gap.

Since gcd(5,256,142)=1, [the protected-fringe theorem](CONGRUENCE.md)
gives every odd residue as a **lower** count of this exact +2 switch
even with independent frozen leading and trailing address words
(of the stated opposite parities). Two fully free flanks replace the
previous strip u≈w from the singleton argument. This does not alone
prove a density result. For example, addresses (81,2) give the actual
twin prime counts 53549,53551. No infinitude of twin primes is implied.

## Choose both addresses from the absolute target

There is an unconditional approximation guarantee for the *same* +2
switch: for every sufficiently large prescribed N, choose an exact
+2 pair whose lower member lies within O(N^(2/5)) of N. This improves
the elementary O(sqrt(N)) placement bound for this free-seam family.
It is a theorem about absolute placement, not merely matching a residue.

Take a power of two Y=2^t comparable to N^(1/10), t>=2. Write
u=Ya+1, w=Yb, where a,b>0 will be chosen as integers. The exact
suffix expansion reduces the lower read to

    H = A a b + B a + C b + 223 + 5(F(a)+F(b)+[a odd]+[b odd]),
    A = 256Y^2,
    B = 142Y + 5(Y/2 - 1),
    C = 398Y + 5(Y - 1).

The last term is positive and O(a+b+1), not O(u+w).
The main part equals A(a+p)(b+q)+E, where
p=C/A, q=B/A, E=223-BC/A. Thus p,q=O(1/Y), E=O(1),
with absolute constants. Set m to the smallest integer at least
2 sqrt((N-E)/A)-p-q, and take a+b=m. Let s=m+p+q. The
real solutions for a are

    a = s/2 - p ± sqrt(s^2/4 - (N-E)/A).

By minimality of m the square root is O(sqrt(m)), and for large N
both a and b lie in the interior of the positive interval. Round
a to its nearest integer and use b=m-a. The product error is
O(A sqrt(m)) because the derivative along that line at the real
root is at most 2A times that square root; the second-order term
is bounded by A/4. Evaluating the retained small term costs only
O(m). Now m=Theta(N^(2/5)), A=Theta(N^(1/5)), so both errors are
O(N^(2/5)). This is constructive with integer comparisons on the
quadratic and the exact recurrence for F; exact square roots need
not be computed numerically. Choosing nearby test inputs with a
larger constant in the error term also brackets a given target
from both sides.

Here the two flanks do not play the same role as two independent
remainder knobs: fixed low suffixes make the absolute equation
nearly a shifted product, and a near-square factor split chooses
its center. The fixed odd/even signs retain the exact +2 while
making that approximation. The remaining Chan–Pak step would
replace the residual of order N^(2/5) by an exact correction,
without consuming the signed switch. No such universal correction
is shown here.
