# One physical comparison word for three simultaneous congruences

This adapts a distinction in [Sun, §§6–7](https://arxiv.org/html/2611.14605v1#S6).
The common comparison point is chosen against the joint congruence
conditions, not chosen independently for each detector. The real and
height bounds are a separate ingredient; they are what make exact
arithmetic conclusions possible. In this tournament chart we can prove
the common-witness part completely, and prove a quantitative obstruction
to simply adding a small uniform height bound.

## What can be prescribed together

Use the one-apex family whose puncture is a transitive spine.
The contact word w has n bits, with 1 meaning the apex beats that
spine vertex. Require first bit 1 and last bit 0, so every spine
vertex lies in a directed cycle through the apex and the whole graph
is strong. Write P for the unsigned binary address and rev_n(P)
for its n-bit reversal. Then the three observables are

    H=g(P),     R=P+1,     R_c=2^n−rev_n(P),

where R and R_c are the positive root responses of the apex in
the tournament and its converse. For the first response the
empty old part contributes 1 and a first contacted selected
vertex at exponent i contributes 2^i. Converse gives the
reverse-complement word. The ordinary count obeys

    g(0)=1,  g(2x+1)=g(x),  g(2x)=g(x)+x+(x mod 2).

**Joint residue theorem.** Fix any modulus M≥1 and a desired
spine-width parity. At a single width of that parity, one physical
strong tournament exists for every triple of residues (H,R,R_c)
subject exactly to:

* If 2 divides M, all three residues are odd.
* If 3 divides M, R_c≡R (mod 3) for even width, and R_c≡−R for
  odd width.

There are no further local constraints, including for moduli with
powers of 3. The required width is O((m+k)(k+ord_m(2)+1)) for
M=2^k m with m odd, treating ord_1(2)=1. In particular it is
O(k²) for pure powers of two. One common graph works also for a
finite list of moduli by using their least common multiple.

Necessity of parity follows from Rédei and the endpoint bits.
Modulo 3, 2 is its own inverse, so rev_n(P)≡2^(n−1)P.
Substitute above to obtain the stated relation.

## A modular switch with fixed exposed addresses

Use exponents counted from the right. In a fixed-width block compare
the word with 1s at b,b+d and the word with 1s at b,b+d+q.
Choose q≥max(k,2) a multiple of ord_m(2), b≥max(k,1), and
enough leading zeros to put all the moving bits above position k−1
in both readings. Use enough room for every d≤max(k,1).
The two blocks have equal P and equal reverse address modulo M.
In the direct reading the difference is 2^(b+d)(2^q−1);
in the reverse reading the same odd factor is multiplied by
a power of two at least 2^k. Their first and last bits agree.

The difference in ordinary g depends on their internal pair, not
on the surrounding word. One way to see this is the weighted
inversion form of g−1. A 1 before a 0 at distance j contributes
2^(j−1), with one extra if j=1. Cross-block contributions factor
as the binary address of the left block times the reverse-complement
address of the right, with a power of two for any gap between
them. The extra adjacency term depends only on unchanged boundary
bits. Thus replacing a block changes the whole H modulo M by
just its own change Δ, even in an arbitrary surrounding word.

For a block with two 1s the sum of terms from one 1 at i to
every lower position, before removing the other 1, is 2^i.
The internal pair would be subtracted with weights 2 at distance
1 or 2^(j−1) at distance j≥2. Therefore:

* When d=1, Δ=2^(b+1)(2^q−1)−2^q+2 is even and equals 1 mod m.
* For d=s+1≥2, Δ=2^(b+d)(2^q−1)+2^(d−1)−2^(d+q−1)
  is 0 mod m and equals 2^s mod 2^k when 1≤s<k.

After other coordinates are chosen, at most m−1 adjacent switches
set the odd part of H. Then the k−1 weighted switches correct every
remaining bit above parity, without changing the odd part. These
are literal words, not hypothetical additive corrections.

## Set both responses first, in that same word

Reserve max(k,1) fixed bits at each end. Choose the last k bits
of P and the last k bits of its reversal according to the two
desired response residues. This fixes first bit 1 and last bit 0.
For k=0 just reserve those two endpoint bits. Reserve switches in
the middle, initially in their first options, and at least qm
additional interior positions for address control. Add at most
one inert zero position to achieve the prescribed width parity.

In the odd part, report the second address after dividing the
reversed value by the unit 2^(n−1). Use positions with exponents
0 or 1 modulo q. A 1 placed at either kind contributes respectively
(1,1) or (2,1/2) to the two addresses. There are m available
positions of each kind. After subtracting the frozen bits, choose
numbers x,y of the two kinds using

    x+2y = A,      x+y/2 = B        (mod m).

The only consistency condition is that 3y=2(A−B) be solvable,
and that is exactly the mod-3 relation if 3 divides m.
Use least residues to stay within the available slots. All these
positions are beyond k bits in both readings, so the even parts
are unchanged. Finally choose the switches for H as above. They
preserve both response residues. All choices take place at the
one fixed width. This proves the theorem constructively.

## Why this cannot substitute for absolute equality

For strong words the first-to-last inversion gives H≥1+2^(n−2)
(with n=2 immediate); hence R<4(H−1). At each exact R the
entire normalized word is forced. Thus at height H≤X there are
at most 4X possible absolute triples in this chart. The joint
residue theorem supplies (M/gcd(2,M))³ compatible residue
triples, divided by 3 if 3 divides M. To supply every such triple,
some output must have H≥M³/96. In particular a uniform O(M)
bound on ordinary count is impossible for this joint construction,
regardless of how cleverly the words are encoded.

For example (H,R,R_c)=(63,65,65) passes this entire local
classification at even width at every modulus. But at exact R=65
the word has P=64 and g(64)=65, not 63. This contrast uses an
attainable ordinary target, not just the known scalar gaps 7 and 21.
We can also impose all these local conditions on H=7 without
getting a scalar witness. This shows where the analogy stops:
a physical solution to a simultaneous residue problem is still
missing the mechanism that controls absolute error.

Run **node replay_joint_word.js** for 992 exhaustive compatible
triples through modulus seven, larger composite-modulus checks,
and literal subset-path replays on small blocks. For example
**node replay_joint_word.js 6 63 65 65 0** constructs one common
word and prints the actual counts and responses.
