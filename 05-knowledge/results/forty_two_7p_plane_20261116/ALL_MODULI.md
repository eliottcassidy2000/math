# Mixed dyadic-odd modular steering of the exact sevenfold lifts

The same physical induced lifts in INTERTWINER.md can realize every
allowed residue of the parent H modulo **any** M, including even
moduli (allowed means odd if M is even). This strengthens the odd
modulus statement. It is still not exact absolute placement.

For either normalized parent, p=A+13 with A=28 or 30 and g=13.
The inherited binary controller, with the word chosen from the
reverse complement of the binary expansion of any z>=1, reads

    H = p + A z + 13(f(z)−1)
      = A(z+1) + 13 f(z),

where f(0)=1, f(2u+1)=f(u) and f(2u)=f(u)+u+(u mod2) for u>=1.
These are physical words ending in a near-sink; the induced
larger count is identically 7H.

Put M=2^a m where m is odd (a may be zero). Work with an even
nonnegative u and v=u+1, and for k>=1 take
T_k(u)=2(2^k(u+1)−1). The f-law gives

    v' = 2^(k+1) v−1,   f(u') = f(u)+2^k v.

In order k=1,...,a−1, apply or skip T_k so the kth bit of
f(u) takes any chosen value. The increment toggles that bit and
fixes all lower bits since v is odd. Later use only k>=max(1,a),
preserving f modulo 2^a. It therefore has any chosen odd residue
there.

For m>1, let d be the order of 2 mod m. A k>=max(1,a) with
k+1=0 mod d makes v'=v−1 mod m. Another with k+1=1 mod d
fixes v=1 and adds 1 to f modulo m at v=1. Translate to v=1,
adjust f by this second move, then translate to any chosen v,
precompensating the deterministic f-change in that last transport.
This produces independently any f and v residues modulo m while
preserving the previous dyadic programming. If m=1 skip this.

Finally let L>=max(1,a), z=2^L v−1. Then f(z)=f(u), since
the last L binary digits are all 1. Thus H=2^L A v+13 f(u).
Modulo m, gcd(A,13)=1 and 2^L is a unit, so the independently
selected two registers reach every H residue. Modulo 2^a the
result reads 13 f(u), which reaches all the required odd classes.
CRT completes the physical residue selector. Nothing in this
argument bounds the size of the integer H near a target; e.g.
the controlled v is an unbounded integer while only its phase
is specified.
