# Complementary finite suffixes preserve an exact sum

Let q(x)=p(x)/2 as in the single-address note. For any odd m there is
a fixed length L such that, from **any** two high addresses x,y,
appending complementary L-bit suffixes selects X mod m and
q(X)+q(Y) mod m independently, while always forcing the exact sum

    X+Y = 2^L(x+y+1)-1.

In particular for fixed P≡−1 mod 2^L, every full aligned block
of X between 0 and P supports all prescribed modular pairs under
the exact constraint Y=P−X. The high initial sum need not have
any residue condition.

**Proof.** After one complementary digit, the running sum x+y is
odd. Keep w=x+y+1, and at a position where x is odd, set

    z=q(x)+q(y)+(x+1−w)/2 mod m.

An address gets an increment to q only when appending 0, of
(x+(x mod 2))/2 to that address. Under complementary append,
w doubles. A=1 and B=01 both end at odd x, and act by

    A(x,w,z)=(2x+1,2w,z),
    B(x,w,z)=(4x+1,4w,z+x+1−w/2).

For B the first digit 0 adds (x+1)/2 to q(x), and the second digit 1
(on x, hence 0 on y) adds w−x to q(y); the state formula follows by
substituting the final x and w. All displayed divisions are in the
ring mod odd m. As in the single-address compiler, A has order
o=ord_m(2) and D=B∘A^-2 is a physical block:
    D(x,w,z)=(x−2,w,z+(x+1)/4−w/8).
Its m-fold repeat acts identically, as can be checked by summing all
x arguments and observing that m is odd.

Set E=A∘D∘A^-1. Both E and D² change x by −4 and preserve w. Their
increment difference provides a pure vertical return:

    D^-2∘E : z ↦ z+(1−3x)/8+3w/16.

At x=−1+w/2 this is the unit 1/2 mod m. All inverse operations
can be read as positive repetitions modulo m. In particular the
physical word for D has length divisible by o: its A^-2 part has
length −2 mod o and its B part length two. Hence the whole vertical
return word also has length divisible by o.

Choose an initial fixed length l with 2^l≥2m. The end w of this
normalizing suffix is known before choosing the bits (it depends
only on l). Its binary digits freely choose x≡−1+w/2 mod m, with
actual last bit 1, by choosing the odd representative mod 2m.
This also ensures the actual running sum is odd. A final fixed
l-bit suffix freely sets the desired x residue, with previously
computable q-sum drift: modular q increments depend on x,y mod 2m
and on the bits. If one uses the final suffix with the same residues
but a different parity of x at a stage, the drift can change; the
normalization and every subsequent programmed block end at odd x,
so the fixed x and w mod m and the known opposite endpoint parities
uniquely determine both x,y mod 2m before the final suffix.

Use up to m−1 vertical returns to adjust z for the desired final q
sum. Pad to the full budget of (m−1) returns with copies of A^o,
which preserve x,w,z mod m as well as the known endpoint parities.
Append the final suffix. Independent q-sum and x control as claimed,
with a fixed O(m³+log m) suffix; m=1 is trivial.

## A boundary of the construction

At even modulus this odd-unit argument does not directly lift.
If x≡−1,y≡0 mod 2M and a complementary reset suffix has defect
address R_final≡0 mod 2M (write 1 for a defect, i.e. a zero on x),
its paired surcharge change is

    Δ[q(x)+q(y)] ≡ R_final/2 − Σ_{defects} (R+(R mod 2)) mod M,

where R is the preceding defect address. Indeed the per-step change
is (1−2e) ceil(R/2) for defect bit e, and summation uses
Σ R=R_final−Σ e, and Σ(R mod 2)=Σ e as long as the
last bit is 0 (which follows from the reset). Thus all those
reset shifts are even modulo an even M. This is a limitation of
that normalization and fixed-sum reset, not a no-go theorem for
all sum-preserving methods.

## The parity coupling pays its own cost: all moduli for placement

**Theorem.** For every M=2^a m (m odd), there is a fixed-length
complementary suffix compiler from any two high prefixes
independently prescribing (X mod m, q(X)+q(Y) mod M), with the same
exact sum law. It does not prescribe X mod 4 independently.

For a>0 begin with two fixed complementary bits, e.g. 11 on x.
Thereafter w=x+y+1 is divisible by 4. If e is appended to y,
the increment of the combined q is
    (1-2e) ceil(y/2)+e w/2.
On the other hand floor(y_new/2)=y, hence modulo 2 its increment
is ceil(y/2). Thus q(x)+q(y)-floor(y/2) mod 2 is an invariant after
these two bits. The final y mod 4 can always be chosen to set the
desired parity of the combined q, independently of final x mod m.

The odd-state normalization above can now set x=-1+w/2 mod m and
x=-1 mod 2^(a+1) simultaneously, using a sufficiently long fixed
length of freely chosen complementary bits, with w known from the
length. The corresponding y is 0 mod 2^(a+1). Do the odd vertical
programming as before, padded to fixed length. Then pad with a fixed
word of all ones of length divisible by ord_m(2) and >=a+1. This
preserves the odd state (both final x and w and z), and restores x=-1
and y=0 modulo 2^(a+1).

To program the remaining a−1 dyadic digits of combined q, use the
isolated two-zero-defect reset words from the single-address note.
Use d=j+2 for each j=0,...,a-2. Append enough trailing ones to
make each word length divisible by ord_m(2); the defect address
remains a multiple of 2^(a+1)m. Thus x,w and endpoint parities all
reset modulo 2M, so repeating each word m times kills its odd
q change. By the paired defect formula its dyadic change is
−2m·2^j mod 2^a. Those steps realize every even residue. At the
normalized state, padding by all ones in multiples of the odd
period preserves all modular outputs and allows fixed total
length. (For m=1 omit the odd stage and use period 1.)

The final fixed-length free suffix chooses x mod m and y mod 4
jointly (since the exact final sum is known and m is odd).
Its surcharge on q is deterministic because the state x,y mod
2M before it is fixed by the preceding construction. Precompensate
the odd and dyadic q adjustments for that drift. The targeted
pre-final parity is automatically correct by the parity invariant:
every valid execution has that one pre-final parity at y=0 mod 4,
and the chosen suffix changes it to the target parity. This proves
the independent placement coordinates for arbitrary M.

Consequently the common-low-bits two-tape tangent bracket can
also be made congruent modulo 2M for arbitrary fixed M. In that
scheme choose the low-bit radix q=2^r with r>=a+1; all ordinary
bilinear and linear terms then lose their dependence on the
programmed high address mod 2^a. For the odd part m, choose
u=0 mod m and the high sum P=-1 mod m. Choose the common low
address v by B+(T+S)v-Sq≡N mod g, g=gcd(45,m), using
gcd(T+S,45)=1. After division by two of the full congruence
mod 2M, the combined high q has coefficient 45 and is freely
programmable mod M. The fixed alignment cost in P and the
absolute N^(2/5) bound are unchanged for each fixed M.
