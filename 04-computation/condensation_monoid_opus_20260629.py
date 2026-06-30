"""
THE deep extension: tournament iso classes = FREE MONOID on strongly-connected "primes"
under transitive composition (X (+) Y = X dominates Y, the condensation/module product).
  - complement tau = REVERSE the word AND complement each letter.
  - SC (self-converse) iso classes = tau-palindromic words.
  - the strongly-connected tournaments are the IRREDUCIBLES (primes); the regular/resonant
    ones (max residual) are large primes => incompressible.
Recursive compression: iso class -> condensation word (sequence of SC primes) + each prime.
Verify: #iso(n) = free-monoid count over SC(m); tau=word reversal; the prime ladder.
"""
from itertools import permutations, combinations

def pairs(n): return list(combinations(range(n),2))
def adj_from_bits(n,bits):
    pr=pairs(n); A=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(pr):
        if bits[k]: A[i][j]=True
        else: A[j][i]=True
    return A
def canon(n,bits):
    pr=pairs(n); idx={p:k for k,p in enumerate(pr)}; best=None
    for perm in permutations(range(n)):
        out=[]
        for (a,b) in pr:
            oi,oj=perm[a],perm[b]
            out.append(bits[idx[(oi,oj)]] if oi<oj else 1-bits[idx[(oj,oi)]])
        t=tuple(out)
        if best is None or t<best: best=t
    return best
def complement(bits): return tuple(1-b for b in bits)

def sccs(n,bits):
    A=adj_from_bits(n,bits)
    reach=[[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n): reach[i][i]=True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                rk=reach[k]
                ri=reach[i]
                for j in range(n):
                    if rk[j]: ri[j]=True
    comp=[-1]*n; c=0
    for i in range(n):
        if comp[i]<0:
            comp[i]=c
            for j in range(i+1,n):
                if comp[j]<0 and reach[i][j] and reach[j][i]: comp[j]=c
            c+=1
    # order components by domination (transitive condensation)
    groups={}
    for i in range(n): groups.setdefault(comp[i],[]).append(i)
    # comp X dominates Y if some/all edges X->Y; tournaments => condensation transitive
    order=sorted(groups, key=lambda g: -sum(A[i][j] for i in groups[g] for j in range(n) if comp[j]!=g))
    return [groups[g] for g in order]

def condensation_word(n,bits):
    comps=sccs(n,bits)
    word=[]
    pr=pairs(n); idx={p:k for k,p in enumerate(pr)}
    for grp in comps:
        m=len(grp);
        if m==1: word.append((1,(0,))); continue
        sub=[grp[t] for t in range(m)]; relabel={v:t for t,v in enumerate(sub)}
        subbits=[]
        for (a,b) in pairs(m):
            va,vb=sub[a],sub[b]
            if va<vb: subbits.append(bits[idx[(va,vb)]])
            else: subbits.append(1-bits[idx[(vb,va)]])
        word.append((m,canon(m,tuple(subbits))))
    return tuple(word)

def all_iso(n):
    pr=pairs(n); m=len(pr); seen={}
    for x in range(1<<m):
        bits=tuple((x>>k)&1 for k in range(m))
        c=canon(n,bits)
        if c not in seen: seen[c]=bits
    return seen

print("n | #iso | #SC-primes | #decomposable | free-monoid check")
SC={}  # n -> set of SC canons
isoall={}
for n in range(1,7):
    iso=all_iso(n); isoall[n]=iso
    scset=set(); dec=0
    for c,bits in iso.items():
        comps=sccs(n,bits)
        if len(comps)==1: scset.add(c)
        else: dec+=1
    SC[n]=len(scset)
    # free monoid: #iso(n) = sum over compositions of n of prod SC(parts)
    from functools import lru_cache
    print(f"{n} | {len(iso):>3} | SC(n)={len(scset):>2} | decomposable={dec:>2}")
print(f"\nSC-prime ladder SC(n) for n=1..6: {[SC[n] for n in range(1,7)]}  (the IRREDUCIBLES)")
# verify free-monoid identity #iso(n)=sum_compositions prod SC
def free_count(n):
    from functools import lru_cache
    @lru_cache(None)
    def f(m):
        if m==0: return 1
        return sum(SC[k]*f(m-k) for k in range(1,m+1))
    return f(n)
print("free-monoid predicted #iso:", [free_count(n) for n in range(1,7)], " vs actual:", [len(isoall[n]) for n in range(1,7)])

# complement = reverse condensation word + complement each letter
print("\ncomplement = reverse-word + complement-letters (check on n=4,5):")
for n in [4,5]:
    ok=True
    for c,bits in isoall[n].items():
        w=condensation_word(n,bits)
        cc=complement(bits); wc=condensation_word(n,cc)
        # reverse w and complement each SC letter
        wrev=tuple((m, canon(m, complement(letter))) for (m,letter) in reversed(w))
        if wrev!=wc: ok=False
    print(f"  n={n}: complement==reverse+complement-letters : {ok}")

# n=4 explicit condensation words (the Klein four as words)
print("\nn=4 classes as condensation words (free-monoid over primes a=triv, C3, S):")
NM={(0,1,2,3):'T',(0,2,2,2):'+',(1,1,1,3):'-',(1,1,2,2):'S'}
def scores(n,bits):
    s=[0]*n
    for k,(i,j) in enumerate(pairs(n)):
        if bits[k]: s[i]+=1
        else: s[j]+=1
    return tuple(sorted(s))
for c,bits in isoall[4].items():
    w=condensation_word(4,bits)
    sizes=[m for m,_ in w]
    print(f"  {NM[scores(4,bits)]}: condensation sizes {sizes}  ({len(w)} prime-syllables)")
