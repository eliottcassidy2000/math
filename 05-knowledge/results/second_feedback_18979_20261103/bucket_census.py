"""Combine inherited packed projective responses with degree-bucket scheduling."""
from collections import defaultdict
from math import factorial,gcd
import importlib.util,time,sys,hashlib
from support import HERE,ROOT,need,save
from integer_census import census

def shift(x,e):
    if e>=0:return x<<e
    need(x&((1<<(-e))-1)==0,'exact power-of-two shift');return x>>(-e)

def count(n,placements=False):
    facts=[factorial(i) for i in range(n+1)];twos=[(f&-f).bit_length()-1 for f in facts];odds=[f>>v for f,v in zip(facts,twos)]
    width=max(1,n.bit_length());digit=(1<<width)-1;masks=[(1<<(width*c))-1 for c in range(n//2+2)]
    states={0:{0:1}};rows=[];neg=0;transitions=0
    for p in range(n if n%2 else n-1,2,-2):
        idx=(p-3)//2;row=sum((j-gcd(p,j))<<(width*i) for i,j in enumerate(range(3,p,2)))
        new={};tr=0
        for d,bucket in states.items():
            r=n-d
            for k in range(r//p+1):
                nd=d+p*k;rr=n-nd;c=max(0,(min(p-2,rr)-3)//2+1);mask=masks[c];offset=k*row
                target=new.setdefault(nd,defaultdict(int));tr+=len(bucket)
                if k==0:
                    for response,w in bucket.items():target[response&mask]+=w
                    continue
                if placements:
                    factor,rem=divmod(facts[r],facts[rr]*p**k*facts[k]);base=k*(p-1)//2+p*k*(k-1)//2+rr*k
                else:
                    factor,rem=divmod(odds[nd],odds[d]*p**k*odds[k]);base=k*(p-1)//2+p*k*(k-1)//2-twos[k]+rr*k
                need(rem==0,'integer cycle placement transfer')
                for response,w in bucket.items():
                    e=base-k*((response>>(width*idx))&digit);neg+=e<0
                    z=shift(w*factor,e) if placements or e>=0 else shift(w,e)*factor
                    target[(response+offset)&mask]+=z
        states=new;transitions+=tr;rows.append(sum(map(len,states.values())))
    total=0
    for d,bucket in states.items():
        r=n-d;need(set(bucket)=={0},'all future coordinates consumed');w=bucket[0]
        if placements:total+=w<<(r*(r-1)//2)
        else:
            factor,rem=divmod(odds[n],odds[d]*odds[r]);need(rem==0,'integer unit-tail transfer')
            total+=shift(w*factor,r*(r-1)//2-twos[r])
    ans,rem=divmod(total,facts[n] if placements else odds[n]);need(rem==0,'integral Burnside output')
    return ans,{'peak_states':max([1]+rows),'stage_states':rows,'transitions':transitions,'negative_shift_checks':neg}

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    incoming=module('incoming_projective',ROOT/'05-knowledge/results/session175_tournament_spectrum_20261104/a000568/projective.py')
    prior=module('prior_bucket',ROOT/'05-knowledge/results/session87_four_vertices_enumeration_20261028/fast_counts.py')
    for n in range(31):
        a,_=count(n);b,_=count(n,True);c=census(n)[0];need(a==b==c,'independent small integral gauges')
    exact=[];timings=[]
    for n in (100,150,200):
        results=[]
        for label,f in [('incoming_packed_scalar',lambda:incoming.packed_count(n)),('bucket_odd',lambda:count(n)),('bucket_placements',lambda:count(n,True)),('prior_bucket_full_prefix',lambda:prior.degree_bucket_dp(n))]:
            t=time.perf_counter();v,stats=f();elapsed=time.perf_counter()-t
            if isinstance(v,list):v=v[-1]
            results.append((v,stats));timings.append({'n':n,'engine':label,'seconds':elapsed})
            print(n,label,round(elapsed,3),flush=True)
        need(len({a for a,s in results})==1,'all four exact engines agree')
        a,ss=results[1];b,tt=results[2]
        need(ss['stage_states']==tt['stage_states'] and ss['transitions']==results[0][1]['transitions'],'same inherited projective carrier')
        exact.append({'n':n,'value_sha256_hex':hashlib.sha256(hex(a).encode()).hexdigest(),'bucket_odd':ss,'bucket_placements':tt,
            'incoming_peak':results[0][1]['peak_states'],'incoming_transitions':results[0][1]['transitions']})
    save('bucket_census.json',{'status':'FINITE_EXACT_FOUR_ENGINE_COMPARISON','small_orders_checked':list(range(31)),'large_orders':exact,
        'scope':'new implementation combination of inherited projective quotient and inherited degree-bucket scheduling; two integral gauges, no new state quotient; prefix baseline does more output work'})
    save('bucket_timings.json',{'status':'OBSERVED_WALL_CLOCK_NOT_A_THEOREM','runs':timings})
    print('PASS four-engine bucket census',flush=True)
if __name__=='__main__':main()
