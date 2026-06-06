def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
P=[n for n in range(2,60) if isprime(n)]
def goldbach_pairs(E): return sorted({(p,E-p) for p in P if p<=E-p and isprime(E-p)})
def lemoine_pairs(O): return sorted({(p,(O-p)//2) for p in P if (O-p)>0 and (O-p)%2==0 and isprime((O-p)//2)})
print("=== even = p+q (unordered/sigma-sym) ; odd = p+2q (ordered/sigma-broken) ===")
for E in [6,8,10,12,14]:
    print(f"  even {E}: Goldbach {{p,q}} = {goldbach_pairs(E)}")
for O in [7,9,11,13,15,21]:
    print(f"  odd  {O}: Lemoine (p,q) p+2q = {lemoine_pairs(O)}")
print("\n=== diagonal (p,p): 2p (double, even) and 3p (triple, odd) — the 2 vs 3 ===")
for p in P[:8]:
    note = " (p=2 special: 3p=6 even)" if p==2 else " (odd, p odd)"
    print(f"  p={p}: double 2p={2*p} (even), triple 3p={3*p}{note}")
print("\n=== the meeting of doubling and tripling: 2p = 3q for primes p,q ===")
meet=[(p,q,2*p) for p in P for q in P if 2*p==3*q]
print(f"  solutions (p,q,n) of 2p=3q: {meet}")
print("  => 6 = 2*3 = 3*2 is the UNIQUE number that is both a doubled prime and a tripled prime.")
print("     6 = the 2-3 commutation point = lcm(2,3) = the hexagonal / pi-3 / dZ=1/6 (S623-S628).")
print("\n=== shared-pair link: even E=p+q linked to odds E+p, E+q (Lemoine reps of the SAME pair) ===")
for E in [6,8,10,12]:
    for (p,q) in goldbach_pairs(E):
        odds=sorted({p+2*q, 2*p+q})
        print(f"  even {E}=({p}+{q}) -> odds {{p+2q, 2p+q}} = {odds}  (O-E in {{{p},{q}}})")
print("\n=== sigma(swap): p+q=q+p ALWAYS; p+2q=q+2p IFF p=q (diagonal=apex). Doubling breaks swap-symmetry => ordered. ===")
