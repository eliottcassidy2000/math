# The actual eleven-block orbit saturates simultaneous local conditions

**Proved.** For every odd integer N prime to 3 and all nonnegative
a,b,c,d, a literal tournament in the Wave-20 eleven-block orbit has
count congruent to N modulo 2^a 3^b 5^c 7^d.
It suffices to use one of seven fixed short seed words, followed only
by repeated empty blocks J. The proof gives a digit-by-digit inverse.
All sufficiently large primes and seven times such primes meet the
hypothesis on N, including arbitrarily prescribed 7-adic phase.
No asserted equality to N or scalar cofiniteness follows.

## The physical alphabet

Use the rooted triple (H,W,d), with H the ordinary count, W the sum
of H(T−P) over paths P starting at the root, and d the same sum with
sign (−1)^(|P|+1). From the physical source (1,1,1), a block Jw
has the [THM-5148 law](../../../01-canon/theorems/THM-5148-one-apex-chain-semigroup-matrix-law-and-its-incompleteness-in-range.md)

    M_w=[g(P),Q+1,-[Q even]; P+1,2^|w|,0; [P even],0,-[w empty]],

where P is the binary contact word (O=1) and Q is its reversed
complement with the declared word width. J denotes the blank w.
Wave 20 uses w empty,O,OO,OZ,OOO,OOZ,OZO,OZZ,OOOZ,OZOZ,OZZZ
after a first run of length at most four. All seeds below have
empty first run. The [inherited finite certificate](../chanpak_waves_20261113/README.md#the-two-root-change-conjecture-hyp-9677-to-224-and-a-run-cap-trichotomy)
and the 13-port family are distinct.

## Seven seeds and finite base phases

Construction-order lists of w and their resulting triples:

    (OO,OZ)       (27,35,1)
    (OO,OOOZ)     (81,143,1)
    (O,OZOZ)      (53,75,1)
    (OOOZ,OZZ)    (301,333,17)
    (O,OZZ)       (33,37,1)
    (OZZ)         (11,13,1)
    (empty,OZZZ)  (39,41,1)

Let K=(J^8−I)/4, and for odd p=3,5,7 use K_p=(J^ell−I)/p
with ell=13,31,48, respectively. These are integer matrices.
Their first rows, respectively over the displayed small modulus:

    K exactly: (14,17,-5), with first row of K² even
    K_3 modulo 3: (1,0,0)
    K_5 modulo 5: (4,1,4)
    K_7 modulo 7: (6,2,5).

The following small finite check is part of the proof and is replayed
by local_11.py. Iterate each of the seven seeds through J^j, 0<=j<624,
modulo 16*9*7. Keep states s with first-row K*s equal to 2 modulo 4,
H prime to 3, and first-row K_7*s a unit modulo 7. They supply all
8*6*7 possible (H mod16,H mod9,H mod7). For the chosen seed, the
31 J phases modulo 5 supply all H mod5 with first-row K_5*s a
unit. Pick the phase for the desired H mod5, and combine its exponent
mod31 with j mod624 by the Chinese remainder theorem. Moving along j
mod624 preserves the other phases and their unit tests. This gives
1,680 proved base cases. The base phase at 3 is mod9, **not mod3**:
the 7-period 48 contributes a factor 3 to the common exponent.

## Infinite lifting and one common witness

Let L=lcm(48,39,31)=19344 and B=J^L. Let y=J^j s be the chosen
base state after combining the phase constraints. At p=2, square
J^8=I+4K to get I+8(K+2K²). Since L/16 is odd and the first
row of K² is even, B=I+8C with first-row C*y = 2 modulo 4.
For every i>=0, B^(2^i)=I+2^(i+3) C_i with C_i=C mod4 by
repeated squaring. Moreover B^t y=y mod4. It follows that
adding 2^i to the exponent t changes H modulo 2^(i+5)
by exactly 2^(i+4). The binary choices therefore lift uniquely
from the selected phase modulo 16.

For an odd p on this list, B=I+p^e C with (p,e)=(3,2),(5,1),(7,1),
where the first row C*y is a unit. At 3 use
J^39=I+9K_3 modulo 27 and L/39=496 prime to 3; at 5 and 7
the respective powers L/31 and L/48 are prime to the modulus.
The binomial expansion proves B^(p^i)=I+p^(e+i) C_i with
C_i=C modp. Addition of k*p^i to t changes H by
k*p^(e+i)*(first row C*y) modulo p^(e+i+1). Thus each new
base-p digit of the target is selected uniquely by one digit of t.

The resulting local conditions on t have pairwise coprime
moduli (powers of 2,3,5,7). Combine them by CRT. The single
physical word is the chosen seed followed by j+Lt empty J blocks;
it is not a list of four unrelated local witnesses. If an exponent
requested in the theorem is smaller than its base phase, lift to
the base phase anyway and then reduce.

## Why this is still far from absolute equality

The locally admissible 67 is already a certified hole of the exact
eleven-block language; nevertheless the decoder finds a common
congruence match at every precision. For instance with b=2,c=d=1,
the least canonical exponent in the chosen phase column is
17,021 at a=8, and 214,503,293 at a=20. These are counts of literal
J steps, not logarithmically sized witnesses. Since J² has
first row (1,2,0), the count grows by at least a factor 3 every
two J steps. More generally, for a fixed base phase and any
forever unattained target, the compatible binary exponent has
infinitely many nonzero digits or else its least exponent would
stabilize at a literal equality. At those depths the positive
word length is a fixed positive fraction of the binary modulus.
This is an absolute-height obstruction to using the modular
inverse as a selector for the original target, not an orbit
obstruction at primes other than the four tested here.

Replay: python3 local_11.py verifies all base phases, every binomial
lifting guard, simultaneous deep sample targets and the 67 examples.
