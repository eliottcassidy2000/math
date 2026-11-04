#!/usr/bin/env python3
"""Exact finite tests of the mod-12 sigma tournament universality construction."""
import itertools,json,math
from pathlib import Path

def require(condition):
    if not condition:
        raise AssertionError("Exact verification failed")

def primes(n):
    sieve=bytearray(b'\1')*(n+1);sieve[:2]=b'\0\0'
    for p in range(2,math.isqrt(n)+1):
        if sieve[p]:sieve[p*p::p]=b'\0'*((n-p*p)//p+1)
    return [p for p in range(2,n+1) if sieve[p]]
P=primes(200000); pool=[p for p in P if p%12==11]
def factors(n):
    out=[]
    for p in P:
        if p*p>n:break
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
    if n>1:out.append(n)
    return out

def order(a,p):
    d=p-1
    for f in factors(d):
        while d%f==0 and pow(a,d//f,p)==1:d//=f
    return d

def sigma_mod(p,a,r):
    if p%r==1:return (a+1)%r
    return ((pow(p,a+1,r)-1)*pow((p-1)%r,-1,r))%r

def labels(n,mask):
    edges={pair:bool(mask>>k&1) for k,pair in enumerate(itertools.combinations(range(n),2))}
    ps=[]
    for j in range(n):
        for q in pool:
            if ps and q<=ps[-1]:continue
            if all((pow(p,(q-1)//2,q)==1)==edges[i,j] for i,p in enumerate(ps)):
                ps.append(q);break
        else:raise RuntimeError('Finite candidate pool exhausted')
    exponents=[]
    for i,p in enumerate(ps):
        outgoing=[j for j in range(n) if j!=i and (edges[i,j] if i<j else not edges[j,i])]
        L=math.lcm(*((ps[j] if p%ps[j]==1 else order(p,ps[j])) for j in outgoing)) if outgoing else 1
        require(L%2==1)
        m=L if L>1 else 3
        exponents.append(m-1)
    common=math.lcm(*(a+1 for a in exponents))-1
    for (i,j),forward in edges.items():
        require((sigma_mod(ps[i],common,ps[j])==0)==forward)
        require((sigma_mod(ps[j],common,ps[i])==0)==(not forward))
        require((sigma_mod(ps[i],exponents[i],ps[j])==0)==forward)
        require((sigma_mod(ps[j],exponents[j],ps[i])==0)==(not forward))
    return ps,exponents

def main():
    require(order(23,11)==1)
    require(sigma_mod(23,2,11)==3 and sigma_mod(23,10,11)==0)
    counts={};maxp=0;maxdigits=0
    for n in range(1,6):
        count=1<<(n*(n-1)//2)
        for mask in range(count):
            ps,exps=labels(n,mask);maxp=max(maxp,max(ps));maxdigits=max(maxdigits,max(len(str(e)) for e in exps))
        counts[n]=count
    local=0
    for p in P:
        if p>500:break
        if p<=3:continue
        for e in range(1,61):
            # independent direct finite sum against residue formula
            require(sum(pow(p,k,12) for k in range(2*e+1))%12==(1+e*(p+1))%12)
            local+=1
    cycle=[(11,40,83),(83,88,179),(179,4,11)]
    for p,a,r in cycle:
        require(sigma_mod(p,a,r)==0)
        value=(p**(a+1)-1)//(p-1)
        require(value%r==0 and value%(r*r)!=0)
        require((value//r)%12==11)
        require(sigma_mod(r,2,p)!=0 or pow(r,(p-1)//2,p)==p-1)
        require(pow(r,(p-1)//2,p)==p-1)
    result={'status':'PASS','labelled_tournaments_by_n':counts,'total_tournaments':sum(counts.values()),'max_prime_label':maxp,'max_exponent_decimal_digits':maxdigits,'mod12_local_tests':local,'cycle':[{'source':p,'exponent':a,'target':r,'order':order(p,r)} for p,a,r in cycle], 'scope':'Finite exact verification, not proof of OPN nonexistence or exhaustive exponent search.'}
    Path(__file__).with_name('certificate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
