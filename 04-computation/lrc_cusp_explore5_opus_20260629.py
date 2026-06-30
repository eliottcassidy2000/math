"""
The Galois/cyclotomic structure of the binding. (Z/14)* = units, ~ Z/6 (cyclic). The binding PAIRS at
the 6 witnesses = the 3 cosets of {+-1}. Residues mod 14 split by 14=2*7: units(bind), even(slack,2-part),
mult-7(slack,7-part), mult-14(KILLER). Confirm.
"""
from math import gcd
n=14
units=[k for k in range(1,n) if gcd(k,n)==1]
# (Z/14)* cyclic? find a generator
def order(g):
    x=g%n;o=1
    while x!=1: x=(x*g)%n;o+=1
    return o
gen=[g for g in units if order(g)==len(units)]
print(f"(Z/{n})* = {units}, |.|={len(units)}=phi(14); cyclic generator(s): {gen} (so ~ Z/{len(units)})")
g=gen[0]; powers=[(pow(g,i,n)) for i in range(len(units))]
print(f"   as powers of {g}: {powers}  => Galois group Z/6 (cyclotomic Q(zeta_14))")
# binding pairs {+-k^{-1}} at witness k
print("\nbinding pairs (the 3 cosets of {+-1} in (Z/14)*):")
pairs=set()
for k in units:
    kinv=pow(k,-1,n)
    pair=tuple(sorted({kinv, n-kinv}))
    pairs.add(pair)
for p in sorted(pairs): print(f"   {{{p[0]},{p[1]}}} = {{+-{min(p, key=lambda x: min(x,n-x))}}}")
print(f"   => {len(pairs)} binding pairs = 3 conjugation cosets {{+-1}},{{+-3}},{{+-5}}")
# residue structure by 14=2*7
print("\nresidue structure mod 14 (the 2*7 split):")
for r in range(n):
    if r==0: tag="KILLER (mult-14, origin, kills the edge)"
    elif gcd(r,n)==1: tag="BIND (unit, holds the razor's edge)"
    elif r%2==0 and r%7!=0: tag="slack (2-part, even)"
    elif r%7==0: tag="slack (7-part, mult-7)"
    else: tag="?"
    print(f"   {r:>2}: {tag}")
print("\n=> 14=2*7: units BIND, mult-14 KILLS, the 2-part & 7-part are SLACK. The covering forces the killer.")
