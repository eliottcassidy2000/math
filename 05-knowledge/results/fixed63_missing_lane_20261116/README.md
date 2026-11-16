# Recovering the missing mod-three lane on the literal 63 core

Keep the eight-vertex tournament (2,16,3,23,165,79,31,111)
induced, with root 1 and normalized (H,g,E,O)=(63,11,26,37).
The singleton-only near-sink rail from [the 7p packet](../forty_two_7p_plane_20261116/README.md)
never has quotient 2 mod 3. Two enlarged alphabets do.

## A cheaper exact recovery: one sink followed by four near-sinks

From the inherited binary controller on this root,
H=52(z+1)+11 f(z), with z the integer read from the reverse
complement of the append word ending in a near-sink.
Here f(0)=1, f(2u+1)=f(u), and f(2u)=f(u)+u+(u mod 2)
for u>0. Use the word 1 followed by k zeros (one ordinary
sink then k near-sinks, k>=1). Now z=2^(k+1)−2 and
f(z)=1+2^k, hence

    H = 115·2^k−41.

This is divisible by 7 exactly when k=1 mod 3; its
quotient is 2 mod 3 exactly when additionally k is even.
At k=4 the literal count is 1799=7·257. This is the
least quotient in the missing class even among **all**
fixed-root binary words: a trailing ordinary sink does
not alter H, and H<=1799 forces z<=33 from f(z)>=1.
Exact checking gives just z=30, i.e. word 10000.
Direct subset DP verifies the full profile after these five
vertices. All old arcs remain unchanged.

On k=4+6t the quotient is (1840·64^t−41)/7.
Since v7(64−1)=1 and 7 does not divide 1840, the valuation
of the quotient difference equals that of the address
difference. Thus this cheap branch already hits every
7-adic phase with quotient 2 mod 3. Also the stronger
[all-moduli selector](../forty_two_7p_plane_20261116/ALL_MODULI.md)
applies directly to this fixed core: its proof only needs
even A and odd g with gcd(A,g)=1, and here A=52,g=11.
Prescribe H=7q modulo 7M, with a parity-compatible q,
to select the quotient modulo M. No exact size bound
follows from this argument.

## Uniform-block recovery and an all-order minimum

For a single uniform near-sink block B the count is
52v+11u, for u=H(B), v=P_B(2) (sum of products of Hamilton
counts over ordered two-colour partitions, empty count 1).
It is divisible by 7 iff v=u mod 7. The identity
v=(-1)^|B| u mod 3 means the quotient is 2 modulo 3
iff B has odd order and u=2 mod 3. The first quotient in
that lane is 6245 (block profile (41,832)); the first prime
is 6353 (block profile (53,844), producing H=44471).
Literal rows are in replay.py.

These minima are all-order **within the uniform-block
architecture**, not minima across arbitrary extensions.
To exclude any candidate through 6353 we only need v<=855.
Under any insertion v at least doubles: extend each
old two-colouring to either colour and use induced
monotonicity on the enlarged colour class. Hence every
block below the cap can be enumerated by canonical
one-vertex augmentation with pruning. Its parent is in
the previous layer. Canonicalizing within degree classes
loses no isomorphism types. The exhaustive standalone
replay checks every such layer and the prime test.

The recorded bank of 15 physical block profiles offers
all quotient residues modulo 15. For each selected block,
append 12t ordinary near-sinks. This multiplies its
two-slot coordinate by 4096^t, leaving u unchanged.
Since (4096−1)/7=585 is divisible by 15 but not by 7,
and v is a 7-adic unit, that branch retains its mod-15
residue and is a bijective isometry modulo 7^a in t.
A fixed branch has only logarithmically many possible
ordinary addresses below an absolute bound. Hensel
inversion to a modulus wider than that address interval
and then literal testing is exact; it is not global
absolute coverage.

## A general uniform amplification bound

If an active normalized source of count h and even
response E can be multiplied by q>1 through a single
uniform block (u,v), then (q−u)h=2(v−u)E.
Let d=gcd(h,E), h0=h/d and e0=E/d. The coprime
integer coordinates satisfy v−u=k h0 and q=u+2k e0
for some positive k (h0 is odd). The inherited
h<=4E−1 gives h0<=4e0−1<=2q−3. Thus fixed q
sees only finitely many primitive ratios under this
architecture. In particular any active prime count p
produced by a near-sink step has E<p and reduced
numerator p, so p>=2q−2 excludes a uniform q lift.
If q<=p, k>=2 violates q>=1+4E>=p+2, leaving
u=q−2E and v=u+p. Except for the singleton block
(which has v−u=1), the path-colour polynomial gives
v=0 mod 4, an extra necessary condition u=−p mod 4.
This applies also to any sequence of uniform near-sink
blocks by ordinal concatenation. It does not constrain
the nonuniform two-vertex sevenfold attachments.

Run python replay.py to verify both minima, literal
witnesses and all bank profiles.
