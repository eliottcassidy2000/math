#!/usr/bin/env python3
"""
SUPER-UNIFORMITY of Paley tournaments.

For Paley T_p (p = 3 mod 4 prime):
- f(d,j) = constant for d in QR (quadratic residue)
- f(d,j) = 0 for d in NQR (quadratic non-residue)

For Paley T_p (p = 1 mod 4 prime):
- f(d,j) = constant for d in QR
- f(d,j) = constant for d in NQR
(since -1 is QR, so both QR and NQR edges exist from vertex 0)

Wait, T_p only exists at p=3 mod 4 (otherwise the tournament is not well-defined
since QR would include half the elements and -1 being QR would make T symmetric).

Actually T_p at p=1 mod 4 IS a tournament: QR has (p-1)/2 elements. The issue is
that -1 is a QR, so d in QR implies -d in QR (since -1 * d = -d and QR is a group).
This means T[i,j]=1 iff T[n-i, n-j]=1, which means the tournament is self-converse.

For p=3 mod 4: -1 is NQR. So d in QR implies -d in NQR.
T[i,j]=1 implies T[j,i]=0 (antisymmetric) - always true for tournaments.

Actually: T_p is a tournament for ALL odd primes p. QR mod p always has (p-1)/2 elements.
For i != j: exactly one of j-i and i-j is in QR (since one is the negative of the other,
and -1 relates QR and NQR differently based on p mod 4).

Wait: for p=1 mod 4, -1 is QR. So if d is QR, then -d = (-1)*d is also QR.
So j-i in QR implies i-j = -(j-i) in QR. But then T[i,j] = T[j,i] = 1, contradiction
for a tournament!

So T_p is a tournament ONLY for p = 3 mod 4.

For p = 1 mod 4: the "Paley graph" is an undirected graph (edges when j-i in QR).

So our Paley tournaments are at p = 3, 7, 11, 13 (wait, 13 = 1 mod 4!).

Let me check: 13 mod 4 = 1. So T_13 should NOT be a tournament!
But we computed H(T_13) = 1,579,968 earlier...

Actually, QR mod 13 = {1,3,4,9,10,12}. Check: is -1 = 12 in QR?
12^2 = 144 = 11*13 + 1, so 12^2 mod 13 = 1. Is 12 a QR? Yes: 12 = -1, and
(-1)^2 = 1. So 12 is QR if and only if there exists x with x^2 = 12 mod 13.
5^2 = 25 = 12 mod 13. Yes! So -1 IS a QR mod 13.

But then d and -d are both QR or both NQR. So for d in QR, both j-i and i-j are QR.
This means T[i,j] = T[j,i] = 1, which is NOT a tournament!

So our T_13 computation earlier must have been using a different generator set.
Let me check...

In our code, we used S = {1, 3, 4, 9, 10, 12} which has 6 = (13-1)/2 elements.
And QR mod 13 = {1, 3, 4, 9, 10, 12}.

For d=1: 1 in S, so T[0,1]=1. For d=12=-1: 12 in S, so T[0,12]=1 and T[12,0]=1.
But also T[1,0]: 0-1 = -1 = 12, which IS in S, so T[1,0]=1. So T[0,1]=T[1,0]=1.
This is NOT a tournament!

So the T_13 computation was wrong — it wasn't computing a tournament.
Let me verify this and understand what went wrong.

kind-pasteur-2026-03-06-S25d
"""

def qr_set(p):
    """Quadratic residues mod p."""
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    return qr

# Check which Paley primes give tournaments
for p in [3, 5, 7, 11, 13, 17, 19, 23]:
    qr = qr_set(p)
    print(f"p={p}: QR = {sorted(qr)}, |QR| = {len(qr)}, p mod 4 = {p % 4}")
    # Check if -1 is QR
    neg1_qr = (p - 1) in qr
    print(f"  -1 = {p-1} in QR? {neg1_qr}")
    # Check if it's a tournament
    is_tournament = True
    for d in range(1, p):
        d_in = d in qr
        neg_d_in = (p - d) in qr
        if d_in and neg_d_in:
            is_tournament = False
            break
    print(f"  Is tournament? {is_tournament}")
    print()
