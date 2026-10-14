#!/usr/bin/env python3
"""Exact, dependency-free checks for the Wieferich/tournament lift note.
No large tournament Hamiltonian enumeration is claimed. Output is JSON.
"""
import json
from collections import Counter
from itertools import permutations, product
from pathlib import Path


def factors(n):
    ans=[];q=2
    while q*q<=n:
        if n%q==0:
            ans.append(q)
            while n%q==0:n//=q
        q+=1
    if n>1:ans.append(n)
    return ans


def order_prime_power(a,p,k):
    n=p**k;d=(p-1)*p**(k-1)
    for q in factors(d):
        while d%q==0 and pow(a,d//q,n)==1:d//=q
    assert pow(a,d,n)==1
    assert all(pow(a,d//q,n)!=1 for q in factors(d))
    return d


def hp(A):
    n=len(A);incoming=[sum(1<<i for i in range(n) if A[i][j]) for j in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for j in range(n):dp[1<<j][j]=1
    for mask in range(1,1<<n):
        for j in range(n):
            if mask>>j&1:
                rest=mask^(1<<j)
                if rest:
                    ins=rest&incoming[j];s=0
                    while ins:
                        bit=ins&-ins;ins-=bit;s+=dp[rest][bit.bit_length()-1]
                    dp[mask][j]=s
    return sum(dp[-1])


def cayley(n,S):
    return [[int(i!=j and (j-i)%n in S) for j in range(n)] for i in range(n)]


def orbit_pairs(n,a):
    left=set(range(1,n));pairs=[]
    while left:
        x=min(left);o=set();y=x
        while y not in o:o.add(y);y=y*a%n
        neg={(-i)%n for i in o}
        assert o.isdisjoint(neg)
        pairs.append((o,neg));left-=o|neg
    return pairs


def odd_cycles(A):
    n=len(A);out=[]
    def visit(path,used):
        x=path[-1]
        if len(path)>=3 and len(path)%2 and A[x][path[0]]:out.append(tuple(path))
        for y in range(path[0]+1,n):
            if not(used>>y&1) and A[x][y]:visit(path+[y],used|1<<y)
    for x in range(n):visit([x],1<<x)
    return out


def ocf(A,tau,p):
    n=len(A)
    assert all(A[i][j]==A[tau[i]][tau[j]] for i in range(n) for j in range(n))
    for i in range(n):
        j=i
        for _ in range(p):j=tau[j]
        assert j==i
    cs=odd_cycles(A);ix={c:i for i,c in enumerate(cs)}
    bits=[sum(1<<x for x in c) for c in cs]
    def canon(c):
        i=c.index(min(c));return c[i:]+c[:i]
    perm=[ix[canon(tuple(tau[x] for x in c))] for c in cs]
    cols=[]
    def rec(col,used,start):
        cols.append(tuple(col))
        for i in range(start,len(cs)):
            if not used&bits[i]:rec(col+[i],used|bits[i],i+1)
    rec([],0,0)
    seen=set();F=D=O=0;profile=Counter();free_count=0
    for col in cols:
        im=tuple(sorted(perm[i] for i in col))
        if im==col:
            u=sum(perm[i]==i for i in col)
            assert (len(col)-u)%p==0
            t=(len(col)-u)//p
            F+=2**(u+t);D+=t*2**(u+t);profile[(u,t)]+=1
        elif col not in seen:
            cur=col
            for _ in range(p):seen.add(cur);cur=tuple(sorted(perm[i] for i in cur))
            assert cur==col
            O+=2**len(col);free_count+=1
    H=sum(2**len(col) for col in cols)
    assert H==hp(A)
    q=((pow(2,p-1,p*p)-1)//p)%p
    assert (H-F-p*(q*D+O))%(p*p)==0
    return dict(odd_cycles=len(cs),collections=len(cols),H=H,F=F,D=D,O=O,
                free_collection_orbits=free_count,
                fixed_profile={str(k):v for k,v in sorted(profile.items())},
                q=q,mod_p2=H%(p*p),F_mod_p2=F%(p*p),
                pqD_mod_p2=p*q*D%(p*p),pO_mod_p2=p*O%(p*p))


def matrix_mask(A):
    return sum(A[i][j]<<k for k,(i,j) in enumerate((i,j) for i in range(len(A)) for j in range(i+1,len(A))))


def from_mask(n,mask):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+1,n)):
        A[i][j]=(mask>>k)&1;A[j][i]=1-A[i][j]
    return A


def small_spectrum(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    ix={pair:i for i,pair in enumerate(pairs)};total=(1<<len(pairs))-1
    counts=[0]*(1<<len(pairs))
    # Each directed ordered path fixes n-1 orientations. Enumerate all completions.
    for path in permutations(range(n)):
        positive=known=0
        for i,j in zip(path,path[1:]):
            b=1<<ix[tuple(sorted((i,j)))];known|=b
            if i<j:positive|=b
        free=total^known;s=free
        while True:
            counts[positive|s]+=1
            if not s:break
            s=(s-1)&free
    strong_values=set()
    for mask,H in enumerate(counts):
        A=from_mask(n,mask)
        def reaches(transpose):
            seen={0};todo=[0]
            while todo:
                i=todo.pop()
                for j in range(n):
                    if j not in seen and (A[j][i] if transpose else A[i][j]):seen.add(j);todo.append(j)
            return len(seen)==n
        if reaches(False) and reaches(True):strong_values.add(H)
    assert 7 not in counts and 21 not in counts
    return dict(n=n,tournaments=len(counts),all_H=sorted(set(counts)),strong_H=sorted(strong_values))


def main():
    arithmetic=[]
    for p in [3,5,7,13,31,1093,3511]:
        d=order_prime_power(2,p,1);e=d&-d;a=2**e;delta=d//e;m=(p-1)//(2*delta)
        orders2=[order_prime_power(2,p,k) for k in [1,2,3]]
        ordersa=[order_prime_power(a,p,k) for k in [1,2,3]]
        q=((pow(2,p-1,p*p)-1)//p)%p
        bits=[sum((p**j-p**(j-1))//(2*order_prime_power(a,p,j)) for j in range(1,k+1)) for k in [1,2,3]]
        row=dict(p=p,ord_2=orders2,q2=q,normalized_base=a,odd_orders=ordersa,m=m,orientation_bits=bits,lex_bits_p2=2*m)
        assert bits[1]==(m*(p+1) if q==0 else 2*m)
        if q==0:
            row['second_Fermat_digit']=((pow(2,p-1,p**3)-1)//p**2)%p
            row['second_order_digit']=((pow(2,d,p**3)-1)//p**2)%p
            H={pow(a,j,p*p) for j in range(delta)}
            assert len(H)==delta and {h%p for h in H}=={pow(a,j,p) for j in range(delta)}
            assert H.isdisjoint({(-h)%(p*p) for h in H})
            assert 1 in H and (1-p)%(p*p) not in H and 1+p not in H
            row['single_orbit_flip_witness']={'vertices':[0,p,1], 'arcs':['1 -> 0',f'{p} -> 1'],'flipped_unit_orbit_size':len(H)}
        arithmetic.append(row)
    small=[small_spectrum(n) for n in range(3,7)]
    near=[]
    for n in range(3,11):
        A=[[int(i<j) for j in range(n)] for i in range(n)];A[0][-1]=0;A[-1][0]=1
        H=hp(A);assert H==2**(n-2)+1;near.append({'n':n,'H':H})
    cyclic9=[]
    pairs=orbit_pairs(9,1)
    for t in product(range(2),repeat=len(pairs)):
        S=set().union(*(pairs[i][v] for i,v in enumerate(t)));H=hp(cayley(9,S));u=H;v=0
        while u%3==0:v+=1;u//=3
        cyclic9.append(dict(S=sorted(S),H=H,v3=v,invariant4={4*x%9 for x in S}==S))
    controls=[]
    tau=[(x+3)%9 for x in range(9)]
    A=from_mask(9,22024714438)
    z=ocf(A,tau,3);z.update(n=9,p=3,tau=tau,mask=matrix_mask(A),adjacency=A);controls.append(z)
    assert (z['H'],z['F'],z['D'],z['O'])==(2715,39,4,888)
    A=[[0,0,1,0,0,0,1,0,1],[1,0,1,1,0,0,1,1,1],[0,0,0,0,0,1,1,1,0],[1,0,1,0,0,1,0,0,0],[1,1,1,1,0,1,1,0,0],[1,1,0,0,0,0,0,0,1],[0,0,0,1,0,1,0,0,1],[1,0,0,1,1,1,1,0,1],[0,0,1,1,1,0,0,0,0]]
    tau=[(x+3)%9 for x in range(9)]
    assert all(A[i][i]==0 for i in range(9)) and all(A[i][j]+A[j][i]==1 for i in range(9) for j in range(i+1,9))
    z=ocf(A,tau,3);z.update(n=9,p=3,tau=tau,mask=matrix_mask(A),adjacency=A);controls.append(z)
    assert (z['H'],z['F'],z['D'],z['O'])==(1557,39,4,502)
    result={'arithmetic':arithmetic,'small_spectra':small,'near_transitive':near,'all_16_cayley_order9':cyclic9,'ocf_controls':controls,
            'scope':'All exact integers; no HP enumeration at orders1093^2 or3511^2; orientation counts are proved orbit formulas.'}
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
