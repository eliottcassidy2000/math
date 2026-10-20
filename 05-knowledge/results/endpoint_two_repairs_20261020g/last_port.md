# Last-port repair and the two-step orbit bound

Research continuation, October 20, 2026. **PROVED** for the marked
presentations and orders stated below. This is not yet an all-order
black-core connectivity theorem.

## Objects and inheritance

Let P=(a,M,z) be a directed Hamiltonian word in a tournament T of order n.
Write R_P(T)=T xor E(P) for reversal of its path arcs and
K_P(T)=T^op xor E(P). A path is blue exactly when its positional reversal
r_P is an anti-automorphism. Blackness is invariant under a move and
under global converse. Three R words Q,S,V therefore give three K words
Q,reverse(S),V with the corresponding complement target.

The inherited mechanism is the two-block repair in
[strong_blue_games_20261020f/detours.md, Section 2](../strong_blue_games_20261020f/detours.md).
The hostile case is a reversed suffix whose prescribed first vertex can
only reach its prescribed last vertex immediately. The corrected near
miss is to mistake an internal Hamiltonian-path relocation for a
metagraph edge. The retained sidecars are the actual extra endpoint
ports and the permutation obtained by composing two anti-automorphisms.

The concept board is: endpoint ports; exact path-edge XOR; local
anti-automorphism witnesses; orbit score capacity; and module entry/exit
parity. The anchor is the (n-3,2) score layer, the niche is port-aware
module lifting, and the wildcard is the orbit capacity obstruction.

For a closed blue P with d_T(a)=n-3 and d_T(z)=2, there is a unique
u in M with u->a and a unique v=r_P(u) in M with z->v. Every other
middle vertex satisfies a->w->z. The closing arc is z->a. This gives a
coordinate description, not an identification of arbitrary tournaments
by their degree sequence.

## L1. A cyclic-cut word cannot be blue in this layer

**PROVED.** Suppose n>=7 and P is closed and blue with d_T(a)=n-3.
Every word with the same positional reversal as the one-step cyclic cut
of P is black, whenever it is a directed Hamiltonian word in T.

If such a word Q were blue, r_Q r_P would be an automorphism shifting
positions by two modulo n. For odd n this is transitive, so T would be
regular. The equation n-3=(n-1)/2 forces n=5.

For n=2m, the automorphism is transitive on the even positions and on
the odd positions. Each induced m-vertex tournament is regular, giving
(m-1)/2 outneighbors within its orbit. The closing arc z->a and its
translates give every even-position vertex at least one incoming arc
from the odd-position orbit. Hence

    d_T(a) <= (m-1)/2 + m-1 = 3(m-1)/2.

But d_T(a)=2m-3, so m<=3, contradicting n>=7. This argument also
excludes even m automatically, since an even-order tournament cannot
be regular. The shift is TWO, not one; ignoring its two orbits would
invalidate the even-order argument.

## L2. Repair when the last middle vertex is the extra port

**PROVED.** Let n>=7 and let

    P=(a,b,x,U,c,z)

be closed and blue, with endpoint scores (n-3,2), u=c, and x->c.
Then the marked blue edge has a three-black detour, reaching its
literal labelled target. No hypothesis that the middle is open or
closed, or that the partner is strong, is required.

The R words are

    Q=(a,x,U,c,z,b),
    S=(a,b,x,c,reverse(U),z),
    V=(b,z,U,c,x,a).

Equivalently, use K words Q,reverse(S),V. The source-to-target map
retains every labelled arc: its path-edge XOR is exactly E(P).

To prove legality, apply the inherited two-block construction with
pivot b, left block (a), and right block (x,U,c,z). In the right block
after reversing its displayed path, a same-endpoint repair is

    H=(x,c,reverse(U),z).

Its first arc is the assumed x->c; its interior arcs are reversed
displayed-path arcs; its last arc is unchanged because z beats only a
and b. These also directly verify the three successive R words.

Q has the same positional reversal as the cyclic cut
(b,x,U,c,z,a), so L1 proves it black. The endpoint score sum of S at
its source is n-2, so S is black. Finally evaluate reverse(V) at the
other end of that move, namely R_P(T); it has the same positional
reversal as V. The arc z->a is unchanged there, whereas
b->x is false because it was a P arc. The reversal of V pairs z with
x and a with b. An anti-automorphism would therefore send z->a to
b->x, a contradiction. Thus V is black, including the cases in which
its endpoint score test vanishes.

## Scope and information losses

The theorem preserves the exact marked complement edge and supplies
three genuine black moves, rather than auxiliary path relocations.
The degree pair alone forgets which vertex u is the extra port and the
chord x->c; those data are indispensable to this displayed construction.
Failure of that chord does not certify a bridge: it calls for the
block transport repairs developed in the companion notes. Orders below
seven require separate controls. Global connectivity and arbitrary
strong-to-strong blue edges remain separate obligations.
