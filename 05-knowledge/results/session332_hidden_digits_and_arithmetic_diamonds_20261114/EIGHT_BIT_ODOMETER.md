# A three-instruction copy of the naturals with a hidden 3+5-bit section

The two odd-index conventions are (a_i=2i+1) and (b_i=2i-1);
(b_{i+1}=a_i). In (a)-coordinates the known omissions 7 and 21
are indices 3 and 10; in (b)-coordinates they are 4 and 11.
The transported multiplication laws are respectively
(istar j=i+j+2ij) and (istar'j=2ij-i-j+1), conjugate by
the translation (imapsto i+1). They are not ordinary addition
or multiplication of the index integers.

## A literal FRACTRAN carry

Represent k by a hidden triple

    k = 256h + 8y + x,     0 ≤ x < 8,    0 ≤ y < 32,
    N(k) = 2^h 3^x 5^y.

The first register h is unbounded; x and y are the low 3 and high
5 bits of the hidden byte. Use these fractions in the given order:

    2/(3^7 5^31),    5/3^7,    3.

If x=7 and y=31, the first fraction clears both and increments
h. If only x=7, the second clears x and increments y. Otherwise
the third advances x. These are the three exact carries of the
represented natural index, so one microstep sends N(k) to N(k+1).
For example the hidden indexes 8, 9, 10, 11 have raw numbers
5, 15, 45, 135 in order. An observer who starts at displayed
number 45 can therefore be starting at index 10. No extra
state-marker primes are needed. The three ratios are distinct
and all required by the specified one-step map, so fewer than
three instructions cannot implement it on this encoding.

This particular copy even self-normalises on arbitrary
nonnegative exponents of 2, 3 and 5: the weighted index formula
increments by one whether or not the counters initially satisfy
the bounds. Whenever x>7, the applied instructions strictly
decrease x by seven until x≤7, and subsequent steps preserve
this upper bound. If y still exceeds 31, the top carry occurs
within at most eight steps and decreases y by 31; repeat until
y≤31. Counters then remain in their canonical ranges forever.
This differs from the d=0 trap of the original PRIMEGAME.
If the first numerator 2 is replaced by 1, the same low
registers make a literal 256-cycle on the fixed-h section,
while the encoded index advances modulo 256.

This cyclic translation must not be confused with the Galois
action of odd units modulo 512. The additive index translation
i→i+1 has order 256, whereas an odd multiplier has order at
most 128 modulo 512: square once to get 1 mod 8, then each
subsequent squaring gains at least one 2-adic digit. The bound
128 is attained by 5 (inductively v₂(5^(2^s)−1)=s+2).
Thus a single Galois multiplier cannot be the whole odometer
in different coordinates. The transported multiplication laws
preserve a distinct structure from this carry.

## Physical hidden-state lift, and why it remains sparse

Translate 2 to a singleton SCC, 3 to a cyclic triangle and 5 to
the unique strong order-four tournament, in an ordinal stack.
Their (H,Q(2)) pairs are (1,2),(3,12),(5,32), respectively. The
lift of N(k) therefore has

    H = 3^x 5^y,
    Q(2) = 2^h 12^x 32^y.

This pair determines the hidden index exactly: recover x,y by
prime valuations of H, then h=v₂(Q)−2x−5y. Each of the three
fractions becomes a literal block carry at the same threshold.
Neither the Hamiltonian count alone nor the raw integer order
is the copied natural number.

If only strong outputs are wanted, wrapping the ordinal graph
in C3[T,1,1] gives a strong graph with H' = H + Q(2). For each
fixed byte (x,y) this gives exactly the shape

    3^x 5^y + 3^x · 2^(h+2x+5y).

(The empty ordinal graph at k=0 can instead be diagnosed
directly or started later.) The wrapped readouts form only 256
exponential families, hence O(log X) distinct values below X.
This is a physical natural index and regeneration mechanism,
not a positive-density spectrum construction. The needed
missing ingredient is a genuinely additive readout of the
weighted hidden index rather than a monomial readout of its
prime-counter encoding.
