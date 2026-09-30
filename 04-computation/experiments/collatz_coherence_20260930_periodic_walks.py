"""(d) corrected: the orbit of the Banach periodic point x_w (a rational with odd denominator) reduced mod 2^k
is a closed walk of the parity graph G_k: consecutive residues r_j, r_{j+1} satisfy r_{j+1} = T(r_j) mod 2^(k-1)
(r_{j+1} is one of the two lifts), the parities follow w, and the walk closes after p steps."""
from fractions import Fraction
def T(x): return x/2 if x.numerator % 2 == 0 else (3*x+1)/2
def periodic_point(w):
    p=len(w); a=sum(w); c=0
    for j,b in enumerate(w):
        if b: c = 3*c + 2**j
    return Fraction(c, 2**p - 3**a)
def res(x, m): return (x.numerator * pow(x.denominator, -1, m)) % m
words = {"0 ({0})":[0], "10 ({1,2})":[1,0], "1 ({-1})":[1], "110 ({-5,-7,-10})":[1,1,0],
         "11110111000 ({-17})":[1,1,1,1,0,1,1,1,0,0,0], "100 (x=1/5)":[1,0,0], "1010 (x=?)":[1,0,1,0], "11100 (x=?)":[1,1,1,0,0]}
for name,w in words.items():
    x = periodic_point(w); p=len(w); orbit=[x]
    for j in range(p-1): orbit.append(T(orbit[-1]))
    assert T(orbit[-1]) == x, "not periodic?"
    assert all((o.numerator%2==1) == bool(w[j]) for j,o in enumerate(orbit)), "parity word mismatch"
    checks=[]
    for k in range(1,9):
        m=2**k; r=[res(o,m) for o in orbit]
        ok = all(((3*r[j]+1)//2 if r[j]%2 else r[j]//2) % (m//2) == r[(j+1)%p] % (m//2) for j in range(p))  # lift adjacency
        ok = ok and all((r[j]%2)==w[j] for j in range(p))
        checks.append(ok)
    print(f"  word {name:24s} x_w = {x!s:>8}: orbit mod 2^k is a closed walk of G_k with parities w, k=1..8: {checks}")
