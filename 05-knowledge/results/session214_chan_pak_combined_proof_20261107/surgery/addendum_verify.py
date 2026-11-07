"""Separate certificate for the two-component structure and 27371 gap."""
import hashlib,json
from counting import adj,count,dp,scan

def need(test,message):
    if not test:raise RuntimeError(message)

def components(out):
    n=len(out);reach=[out[i]|(1<<i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    unseen=set(range(n));comps=[]
    while unseen:
        i=min(unseen);c={j for j in unseen if reach[i]>>j&1 and reach[j]>>i&1}
        unseen-=c;comps.append(c)
    return comps

def induced(out,vs):
    vs=sorted(vs)
    return [sum(1<<j for j,v in enumerate(vs) if out[u]>>v&1) for u in vs]

def run():
    target=27371
    need(target==101*271,'factorization')
    for p in (101,271):need(all(p%d for d in range(2,int(p**0.5)+1)),('prime factor',p))
    residues={str(k):target%((1<<k)+1) for k in range(1,15)}
    need(all(residues.values()),'nearly transitive divisor hostile')
    need(9*5**5==28125 and (1<<15)+1>target,'cutoff arithmetic')
    component_checks=two_component_cases=0
    for L in range(1,7):
        for r in range(-1,L):
            a=[int(i==r) for i in range(L)]
            for ab in (0,1):
                for q in range(1<<L):
                    b=[q>>i&1 for i in range(L)];out=adj(a,b,ab)
                    cs=[c for c in components(out) if len(c)>1]
                    need(len(cs)<=2,('component count',L,r,ab,q))
                    if len(cs)==2:
                        ac=next(c for c in cs if L in c);bc=next(c for c in cs if L+1 in c)
                        need(ac!=bc and all(out[u]>>v&1 for u in bc for v in ac),'a component later')
                        achain=sorted(ac-{L});first=achain[0]
                        need([i for i in achain if out[L]>>i&1]==[first],'one-exit component shape')
                        av=dp(induced(out,ac));bv=dp(induced(out,bc));whole=dp(out)
                        need(av==1+(1<<(len(ac)-2)) and whole==av*bv,'component product')
                        two_component_cases+=1
                    component_checks+=1
    bound_checks=0
    for L in range(3,13):
        for mid in range(1<<(L-2)):
            b=[1]+[mid>>i&1 for i in range(L-2)]+[0]
            h=count([0]*L,b,0) # a is a global sink, so this is the one-apex count.
            need(h>=1+(1<<(L-2)),('one-apex floor',L,mid,h));bound_checks+=1
    witnesses={1:(1,0,1,0)};cases=0;levels=[]
    for L in range(1,18):
        for r in range(L):
            for ab in (0,1):
                for h,q in scan(L,r,ab,target):witnesses.setdefault(h,(L,r,ab,q))
                cases+=1<<L
        levels.append({'L':L,'covered_through_27371':len(witnesses)})
    missing=[h for h in range(1,target+1,2) if h not in witnesses]
    need(missing==[7,21,7477,27371],('complete omitted list',missing))
    for h,(L,r,ab,q) in witnesses.items():
        need(count([int(i==r) for i in range(L)],[q>>i&1 for i in range(L)],ab)==h,('witness replay',h))
    out=adj([int(i in (0,3)) for i in range(11)],[804>>i&1 for i in range(11)],0)
    need(out==[6142,8188,4088,6128,8160,4032,8064,7936,3584,3072,6144,9,2852],'control rows')
    need(dp(out)==27371 and count([int(i in (0,3)) for i in range(11)],[804>>i&1 for i in range(11)],0)==27371,'positive control')
    raw=json.dumps([[h,*witnesses[h]] for h in sorted(witnesses)],separators=(',',':')).encode()
    return {'status':'PASS','target':target,'busch_f20':28125,'factorization':[101,271],
      'power_plus_one_remainders':residues,'omitted_odd_values_through_27371':missing,
      'parameter_cases':cases,'positive_values':len(witnesses),'positive_witness_sha256':hashlib.sha256(raw).hexdigest(),
      'component_checks':component_checks,'two_component_cases':two_component_cases,'one_apex_floor_checks':bound_checks,
      'two_exit_control':{'outgoing_rows':out,'H':27371},'levels':levels}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
