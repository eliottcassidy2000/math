"""Closure-to-denominator representation, prime-minimum certificates, type A."""
from itertools import combinations, permutations
from math import lcm, factorial
import json

def require(c,m):
    if not c: raise ValueError(m)

def closure(S,C,g):
    x=(1<<g)-1
    for T in C:
        if S&T==S:x&=T
    return x

def meetirreducibles(C,g):
    full=(1<<g)-1
    out=[]
    for M in C:
        if M==full:continue
        supers=[K for K in C if K!=M and K&M==M]
        covers=[K for K in supers if not any(L!=K and L&K==L for L in supers)]
        if len(covers)==1:out.append(M)
    return out

def chain_antichain_certificate(M):
    adj={x:[y for y in M if x!=y and x&y==x] for x in M}
    match={}
    def aug(x,seen):
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                if y not in match or aug(match[y],seen):match[y]=x;return True
        return False
    for x in M:aug(x,set())
    nxt={x:y for y,x in match.items()}
    chains=[]
    for x in M:
        if x not in match:
            C=[x]
            while C[-1] in nxt:C.append(nxt[C[-1]])
            chains.append(C)
    left={x for x in M if x not in nxt};right=set();todo=list(left)
    while todo:
        x=todo.pop()
        for y in adj[x]:
            if nxt.get(x)==y:continue
            if y not in right:
                right.add(y)
                if y in match and match[y] not in left:
                    left.add(match[y]);todo.append(match[y])
    anti=[x for x in M if x in left and x not in right]
    require(len(anti)==len(chains),'equal certificate size')
    require(sorted(x for c in chains for x in c)==sorted(M),'chain partition')
    require(all(x&y==x for c in chains for x,y in zip(c,c[1:])),'chain order')
    require(all(x&y not in (x,y) for x,y in combinations(anti,2)),'antichain')
    return chains,anti

def primes(k):
    out=[];x=2
    while len(out)<k:
        if all(x%p for p in out if p*p<=x):out.append(x)
        x+=1
    return out

def code(chains,g):
    ds=[1]*g
    for p,C in zip(primes(len(chains)),chains):
        for i in range(g):
            e=next((t for t,M in enumerate(C) if M>>i&1),len(C))
            ds[i]*=p**e
    return ds

def test(C,g):
    require((1<<g)-1 in C,'top closed')
    require(all(x&y in C for x in C for y in C),'intersection closure')
    M=meetirreducibles(C,g)
    chains,anti=chain_antichain_certificate(M)
    ds=code(chains,g);sq=code([[x] for x in M],g)
    for S in range(1<<g):
        goal=closure(S,C,g)
        for dd in (ds,sq):
            L=lcm(*(dd[i] for i in range(g) if S>>i&1))
            actual=sum((L%d==0)<<i for i,d in enumerate(dd))
            require(actual==goal,'represented closure')
    return M,chains,anti,ds

# All closure systems on a labeled three-element ground set, including nontrivial c(empty).
g=3;full=(1<<g)-1;cnt=0
for sub in range(1<<(1<<g)):
    C={x for x in range(1<<g) if sub>>x&1}
    if full not in C or any(x&y not in C for x in C for y in C):continue
    test(C,g);cnt+=1
print('All',cnt,'closure systems on 3 labeled elements: exact optimal and squarefree codes passed')

# Type A: closure under (i,j),(j,k) => (i,k).
records=[]
for n in range(2,6):
    pairs=list(combinations(range(n),2));g=len(pairs);ix={e:i for i,e in enumerate(pairs)}
    C={m for m in range(1<<g) if all(not(m>>ix[i,j]&1 and m>>ix[j,k]&1) or m>>ix[i,k]&1 for i,j,k in combinations(range(n),3))}
    M,chains,anti,ds=test(C,g)
    top=(1<<g)-1;B={S for S in C if top^S in C}
    inversion={}
    for perm in permutations(range(n)):
        pos={x:i for i,x in enumerate(perm)}
        mask=sum((pos[i]>pos[j])<<k for k,(i,j) in enumerate(pairs))
        inversion[perm]=mask
    require(B==set(inversion.values()) and len(B)==factorial(n),'exact inversion sets')
    for perm,mask in inversion.items():
        expected=set()
        for i in range(n-1):
            pp=list(perm);pp[i],pp[i+1]=pp[i+1],pp[i]
            expected.add(inversion[tuple(pp)])
        actual={T for T in B if (mask^T).bit_count()==1}
        require(actual==expected,'permutohedron adjacency')
        require(top^mask==inversion[perm[::-1]],'complement is total-order reversal')
    row={'r':n,'guards':g,'closed_sets':len(C),'meet_irreducibles':len(M),'minimum_primes':len(chains),'biclosed':len(B),'denominators':ds,'chains':chains,'antichain':anti}
    records.append(row)
    print(json.dumps(row))
with open('universal_certificates.json','w') as f:json.dump(records,f,indent=2);f.write('\n')
print('ALL UNIVERSAL CONTROLS PASS')
