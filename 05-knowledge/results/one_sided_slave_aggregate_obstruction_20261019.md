# Actual obstruction to a four-upper-slave aggregate extension

Date: 2026-10-19. Status: explicit actual-poset family with exact counts. This disproves a proposed aggregate inequality, not the balanced-pair conjecture.

For r>=1, take the poset

    x<a1<a2<...<ar,
    z<ar,
    y isolated,

with no other relations except transitivity. The ordered incomparable pair (x,y) is critical. It has no lower slave, exactly r upper slaves a1,...,ar, and exactly one common incomparable neighbor z. There is no exterior.

Delete y. The remaining word is the fixed (r+1)-chain x,a1,...,ar with z inserted before ar. There are exactly r+1 such words. Reinsert y into r+3 gaps. Consequently

    E=(r+1)(r+3),
    N(x<z)=r(r+3),
    N(z<y)=(r+1)(r+4)/2.

The last count follows by summing the possible rank of z, which is uniformly 1,...,r+1 before y is inserted. Therefore

    p(x<z)+p(z<y)
      = r/(r+1)+(r+4)/(2(r+3)),

    4E-3[N(x<z)+N(z<y)]
      = -(r+4)(r-3)/2.

At r=4 this is an actual seven-point counterexample to the proposed 4c/3 aggregate extension:

    E=35, N(x<z)=28, N(z<y)=20,
    p(x<z)+p(z<y)=48/35>4/3.

It has root degrees (2,6) and c=1. The pair (z,y) remains balanced, with probability 4/7; this is not a counterexample to the 1/3–2/3 conjecture or to a suitably different witness-selection argument.

The aggregate equals 4/3 at r=3 and exceeds it for every r>=4. Thus r=3 is the largest possible universal one-sided-slave threshold for this particular aggregate witness.

The stronger conjecture that common-to-upper attachments only decrease the disjoint baseline also fails earlier. The baseline is 3/2-1/(r+2) for c=1. This family's excess over it is

    (r^2+r-4)/[2(r+1)(r+2)(r+3)],

which is positive for every r>=2. In particular r=2 gives E=15, aggregate numerator 19, and 19/15>5/4. A monotonicity or FKG argument for that baseline would therefore be invalid.
