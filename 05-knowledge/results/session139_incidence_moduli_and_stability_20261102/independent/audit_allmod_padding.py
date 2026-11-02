#!/usr/bin/env python3
"""Independent all-modulus/padding audit; imports no research implementation.
Profiles are counted via permutations and mandatory/optional path cuts,
not the producer's Hamiltonian/subset-partition dynamic program.
"""
from itertools import combinations, permutations
from math import comb, factorial, lcm
import hashlib, json


def need(ok, why):
    if not ok: raise RuntimeError(why)


def graph_profiles(n, bits):
    edges=list(combinations(range(n),2));a=[[0]*n for _ in range(n)]
    for k,(u,v) in enumerate(edges):
        if bits>>k&1:a[u][v]=1
        else:a[v][u]=1
    counts=[0]*(n+1)
    for order in permutations(range(n)):
        mandatory=sum(not a[u][v] for u,v in zip(order,order[1:]))
        optional=n-1-mandatory
        for extra in range(optional+1):
            counts[mandatory+extra+1]+=comb(optional,extra)
    for c in range(1,n+1):
        need(counts[c]%factorial(c)==0,('unordered divisibility',n,bits,c))
        counts[c]//=factorial(c)
    return tuple(counts)


def monomials(pc):
    fall=[1];poly=[0]*len(pc)
    for c in range(1,len(pc)):
        nxt=[0]*(len(fall)+1)
        for j,z in enumerate(fall):nxt[j]-=(c-1)*z;nxt[j+1]+=z
        fall=nxt
        for j,z in enumerate(fall):poly[j]+=pc[c]*z
    return tuple(poly)


def psi_values(pc,M):
    result=[]
    for x in range(M):
        fall=1;total=0
        for c in range(1,len(pc)):
            fall*=x-c+1;total+=factorial(c-1)*pc[c]*fall
        result.append(total%M)
    return tuple(result)


def weighted(pc,M):
    return tuple(factorial(c)*factorial(c-1)*pc[c]%M for c in range(1,len(pc)))


def clock(M):
    if M==1:return 1,1
    n=M;data=[];p=2
    while n>1:
        e=0
        while n%p==0:e+=1;n//=p
        if e:data.append((p,e))
        p+=1
    return max(a for p,a in data),lcm(*( ((p-1)*p**(a-1)) if p!=2 else (1 if a<3 else 2**max(1,a-3)) for p,a in data))


def cutoff(M):
    c=1
    while factorial(c)*factorial(c-1)%M:c+=1
    return c


def padding_states(pc,M,steps):
    k=cutoff(M);v=list(pc[:k])+[0]*max(0,k-len(pc));out=[]
    for r in range(steps+1):
        out.append(weighted(v,M))
        v=[0]+[(c*v[c]+v[c-1])%M for c in range(1,k)]
    return out


def main():
    profiles=set();total=0;checks=0
    classes={M:{} for M in range(1,65)}
    for n in range(1,6):
        for bits in range(1<<(n*(n-1)//2)):
            pc=graph_profiles(n,bits);profiles.add(pc);total+=1
            need(pc[1]%2==1,('Redei parity',n,bits))
            poly=monomials(pc)
            need(all(z==0 for k,z in enumerate(poly) if k%2!=n%2),('polynomial parity',n,bits))
            for M,bank in classes.items():
                f=psi_values(pc,M);s=weighted(pc+(0,)*(6-len(pc)),M)
                if f in bank:need(bank[f]==s,('Psi merge',n,bits,M))
                else:bank[f]=s
                # Recover each forward difference from the cyclic value table.
                for c in range(1,len(pc)):
                    delta=sum((-1)**(c-j)*comb(c,j)*f[j%M] for j in range(c+1))%M
                    need(delta==s[c-1],('Newton recovery',n,bits,M,c))
                checks+=1
    for M,bank in classes.items():need(len(set(bank.values()))==len(bank),('Psi split',M))
    mods=list(range(1,129))+[625,1024,1331,2048];pchecks=0;firsts=[]
    for M in mods:
        A,L=clock(M)
        for pc in sorted(profiles):
            delta=next(k for k,z in enumerate(monomials(pc)) if z)
            bound=max(0,A-delta)
            states=padding_states(pc,M,A+2*L)
            seen={};first=None
            for r,s in enumerate(states):
                if s in seen:first=(seen[s],r-seen[s]);break
                seen[s]=r
            need(first is not None,('missing repeat',M,pc))
            onset,period=first
            need(onset<=bound and L%period==0,('padding onset/period',M,pc,first,bound,L))
            for r in range(bound,len(states)-L):
                need(states[r]==states[r+L],('padding identity',M,pc,r));pchecks+=1
            if pc==(0,1):need(first==(A-1,L),('sharp K1',M,first));firsts.append((M,*first))
    # Literal transient detector row; obtain L21 by repeated singleton padding.
    pc=[0,1]
    for n in range(2,22):pc=[0]+[c*(pc[c] if c<len(pc) else 0)+pc[c-1] for c in range(1,len(pc)+1)]
    E4=(117,1368,5448,8064,2880)
    response=sum(a*b for a,b in zip(E4,pc[1:]))%275
    need(response==172 and (response-117)%275==55,('actual E4 response',response))
    result={'status':'PASS','profile_algorithm':'permutations with mandatory and optional cuts, divided by c!','labelled_tournaments':total,'distinct_integer_profiles':len(profiles),'Psi_tournament_modulus_checks':checks,'Psi_moduli':64,'padding_moduli':len(mods),'padding_period_comparisons':pchecks,'sharp_K1_clocks_sha256':hashlib.sha256(json.dumps(firsts).encode()).hexdigest(),'E4_L1_L21_mod275':[117,response]}
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
