"""
The 7-adic structure of the strong H-ATOM spectrum (kind-pasteur S31f).

Owner lead (developed by mac-mini HYP-2879: 'weight'=atom-count, 49&75 are SINGLE
strong atoms).  My apex-prime value-add: where does the apex prime 7 first enter the
strong-ATOM spectrum?  H=7 and H=21 are FORBIDDEN (not atoms, THM-200/079).  But
49=7^2 IS reported as a single atom.  So the apex prime enters the atom spectrum first
as its SQUARE 49=7^2 -- and at which m?  This ties HYP-2879 (atoms) to HYP-2878
(apex-7 / E_7 odd holes, n=7 = the apex prime).

Method: fixed-base-path enumeration (2^C(m-1,2), covers all iso classes since every
strong tournament has a Hamiltonian path by Redei), filter STRONG (single strong
component = irreducible atom), collect the atom H-spectrum per m.
"""
from collections import defaultdict

def hcount(n, out):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    full = (1<<n)-1
    for mask in range(1<<n):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            w = out[v] & ~mask
            while w:
                u = (w & -w).bit_length()-1
                dp[mask|(1<<u)][u] += c
                w &= w-1
    return sum(dp[full][v] for v in range(n))

def is_strong(n, out):
    def reach(s):
        seen = 1<<s; stack=[s]
        while stack:
            v = stack.pop(); w = out[v]
            while w:
                u=(w&-w).bit_length()-1; w&=w-1
                if not (seen>>u)&1: seen|=1<<u; stack.append(u)
        return seen
    # strong iff from 0 you reach all AND all reach 0 (use reverse)
    rin = [0]*n
    for v in range(n):
        w=out[v]
        while w:
            u=(w&-w).bit_length()-1; w&=w-1; rin[u]|=1<<v
    full=(1<<n)-1
    if reach(0)!=full: return False
    # reverse reachability
    def reachR(s):
        seen=1<<s; stack=[s]
        while stack:
            v=stack.pop(); w=rin[v]
            while w:
                u=(w&-w).bit_length()-1; w&=w-1
                if not (seen>>u)&1: seen|=1<<u; stack.append(u)
        return seen
    return reachR(0)==full

def atom_spectrum(m):
    base=[0]*m
    for i in range(1,m): base[i]|=(1<<(i-1))
    tiles=[(x,y) for y in range(m) for x in range(y+2,m)]
    T=len(tiles)
    atoms=set()
    for cfg in range(1<<T):
        out=base[:]
        for t in range(T):
            x,y=tiles[t]
            if (cfg>>t)&1: out[y]|=(1<<x)
            else: out[x]|=(1<<y)
        if is_strong(m,out):
            atoms.add(hcount(m,out))
    return sorted(atoms)

if __name__=="__main__":
    print("Strong H-ATOM spectrum per m (single strong component = irreducible):")
    all7=[]
    for m in range(3,8):
        sp=atom_spectrum(m)
        sevens=[h for h in sp if h%7==0]
        print(f"  m={m}: atoms={sp}")
        if sevens: print(f"         7-MULTIPLE atoms at m={m}: {sevens}  (49=7^2? {'YES' if 49 in sevens else 'no'})")
        all7+=[(m,h) for h in sevens]
    print("\n7-adic atom structure:")
    print("  H=7 a strong atom anywhere m<=7?", any(h==7 for _,h in all7))
    print("  H=21 a strong atom anywhere m<=7?", any(h==21 for _,h in all7))
    print("  first 7-multiple strong atom(s):", sorted(set((m,h) for m,h in all7))[:6])
    print("  is 49=7^2 a strong atom, and at which m?", [m for m,h in all7 if h==49])
    print("  is 75 a strong atom, and at which m?", [m for m in range(3,8) if 75 in atom_spectrum(m)] if False else "(see spectra above)")
