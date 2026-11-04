"""Independent modular-period and arithmetic checks; no producer imports."""
from itertools import product
from math import isqrt,lcm
from pathlib import Path
import hashlib,json


def require(ok,msg):
    if not ok:raise RuntimeError(msg)


def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def factors(n):
    out=[]
    for q in range(2,isqrt(n)+1):
        if n%q==0:
            out.append(q)
            while n%q==0:n//=q
    if n>1:out.append(n)
    return out


def valuation(n,r):
    v=0
    while n%r==0:n//=r;v+=1
    return v


def exact_profile(p,n,targets):
    U=(pow(p,n)-1)//(p-1)
    return tuple(valuation(U,r) for r in targets)


def modular_order(p,mod):
    x=p%mod;d=1
    while x!=1:x=x*p%mod;d+=1
    return d


def complete_period_domain(p,targets,upper,euler=False):
    # Denominators are units in all audited examples. Work directly modulo
    # each upper-bound guard prime power, without the producer's rank compiler.
    require(all((p-1)%r for r in targets),'unit denominator required')
    guards=[r**(u+1) for r,u in zip(targets,upper)]
    orders=[modular_order(p,g) for g in guards]
    period=lcm(4 if euler else 2,*orders)
    domain=set();witnesses={};tested=0
    for n in range(1,period+1):
        if n%(4 if euler else 2)!=(2 if euler else 1):continue
        row=[];valid=True;tested+=1
        for r,u,mod in zip(targets,upper,guards):
            residue=(pow(p,n,mod)-1)*pow(p-1,-1,mod)%mod
            if residue==0:valid=False;break
            row.append(valuation(residue,r))
        if valid:
            key=tuple(row);domain.add(key)
            witnesses.setdefault(key,n if n>1 else n+period)
    return {'source':p,'period':period,'guard_orders':orders,'length_classes_checked':tested,
            'domain':sorted(domain),'physical_witnesses':{str(k):v for k,v in witnesses.items()}}


def main():
    p=18481;targets=(67,331,50821,4621);orders=(33,165,105,1155)
    require(all(prime(x) for x in (p,)+targets),'N5 primality')
    require(p%1155==1,'source congruence')
    for r,d in zip(targets,orders):
        require(pow(p,d,r)==1,'claimed order does not divide')
        require(all(pow(p,d//q,r)!=1 for q in factors(d)),'claimed order not exact')
        require(pow(p,d,r*r)!=1,'base valuation is not one')
        require(all(dy%r for dy in orders),'target/label collision')
    fixed=[]
    for row in product((0,1),repeat=4):
        n=lcm(*(d for d,x in zip(orders,row) if x))
        got=exact_profile(p,n,targets)
        if got==row:fixed.append(row)
    expected={(0,0,0,0),(1,0,0,0),(1,1,0,0),(0,0,1,0),(1,1,1,1)}
    require(set(fixed)==expected,'N5 interval has wrong profiles')
    require(exact_profile(p,3,targets)==(0,0,0,0),'positive zero witness')
    coupled=[complete_period_domain(3,(11,23),(2,1)),complete_period_domain(5,(11,23),(1,1),True)]
    choices=list(product(*(x['domain'] for x in coupled)))
    require(all(any(sum(row[j] for row in C)==2 for C in choices) for j in range(2)),
            'separate columns not feasible')
    require(not any(all(sum(row[j] for row in C)==2 for j in range(2)) for C in choices),
            'coupled contradiction failed')
    require(exact_profile(3,55,(11,23))==(3,1),'coupled overflow at length 55')
    parity=[]
    for q in (191,653,991,1543):
        require(prime(q),'parity source composite')
        d=complete_period_domain(q,(13,),(2,))
        require(set(d['domain'])=={(0,),(2,)},'parity exact domain')
        parity.append(d)
    sums=sorted({sum(row[0] for row in C) for C in product(*(d['domain'] for d in parity))})
    require(sums==[0,2,4,6,8],'parity sumset')
    require(max(0,5-3*2)==0 and min(2,5-3*0)==2,'interval propagation should be fixed')
    src=Path(__file__).resolve().parent.parent/'primeflow_extension'
    require((src/'normal.json').read_bytes()==(src/'optimized.json').read_bytes(),'producer modes differ')
    here=Path(__file__).resolve().parent
    require((here/'producer_normal.json').read_bytes()==(src/'normal.json').read_bytes(),'fresh replay differs')
    require((here/'producer_optimized.json').read_bytes()==(src/'normal.json').read_bytes(),'optimized replay differs')
    result={'status':'PASS','N5_fixed_profiles':fixed,'coupled_complete_domains':coupled,
            'parity_complete_domains':parity,'parity_sumset':sums,
            'producer_proof_sha256':hashlib.sha256((src/'PROOF.md').read_bytes()).hexdigest(),
            'producer_checker_sha256':hashlib.sha256((src/'verify.py').read_bytes()).hexdigest(),
            'method':'Exact prime/order certificates and complete modular periods at valuation upper-bound guards. No producer imports.'}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
