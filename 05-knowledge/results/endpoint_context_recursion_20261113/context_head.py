"""Complete neutral-context census through8 from the inherited ordinary bank."""
from pathlib import Path
import sys,json,gzip
from fractions import Fraction
from collections import Counter
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'surplus_selector_20261113'))
from context_sequences import canonical
from repair23 import profile,hp2
from inflation_probe import relabel

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def save(name,obj):
    suffix='' if __debug__ else '.optimized'
    (HERE/(name+suffix+'.json.gz')).write_bytes(gzip.compress(json.dumps(obj,separators=(',',':')).encode(),mtime=0))

def ordinary_bank(limit=8):
    raw=json.loads(gzip.decompress((HERE.parent/'multicore_metagraph_20261113/rooted_growth.json.gz').read_bytes()))
    by={}
    for r in raw['records']:
        if r['n']<=limit:by.setdefault((r['n'],r['ordinary_code']),tuple(r['rows']))
    return by

def response(A,r):
    h,s,q,F=profile(A,r);full=(1<<len(A))-1;rev=tuple(full^(1<<v)^A[v] for v in range(len(A)))
    _,e,_,G=profile(rev,r);P=hp2(relabel(A,[v for v in range(len(A)) if v!=r]))[1]
    delta=G-h-q+e;alpha=G+e-2*q;beta=F-P-s
    return dict(profile=(h,s,q,F),ends=e,G=G,P=P,response=(alpha,delta,P,beta),determinant=alpha*beta+delta*P)

def run():
    ordinary=ordinary_bank();levels=[];all_examples={};up=[]
    for n in range(2,9):
        graphs=[T for (size,code),T in ordinary.items() if size==n]
        need(len(graphs)==[0,1,1,2,4,12,56,456,6880][n],'ordinary bank cardinality')
        rooted=set()
        for A in graphs:
            for r in range(n):rooted.add(canonical(A,r))
        records=[];transitions=Counter()
        for A in sorted(rooted):
            r=n-1;data=response(A,r)
            if data['determinant']:continue
            h,s,q,F=data['profile'];alpha,delta,P,beta=data['response'];rho=Fraction(-beta,P);sigma=Fraction(s,q)
            key=str(rho),str(sigma);transitions[key]+=1
            rec=dict(rows=A,root=r,input_slope=str(rho),output_slope=str(sigma),**data)
            records.append(rec);all_examples.setdefault(key,rec)
            if rho<sigma:up.append(rec)
        level=dict(n=n,ordinary=len(graphs),rooted=len(rooted),neutral=len(records),slopes=len(set(x['input_slope'] for x in records)),transitions=[dict(input=a,output=b,classes=c) for (a,b),c in sorted(transitions.items())],records=records)
        levels.append(level)
        print('LEVEL',n,'ROOTED',len(rooted),'NEUTRAL',len(records),'SLOPES',level['slopes'],'CUMULATIVE_RULES',len(all_examples),'UP',len(up),flush=True)
        for target in ('1/3','4/9','6/11'):
            hits=[x for x in records if x['input_slope']==target]
            if hits:print('DEAD_END_REPAIRED',target,'outputs',sorted(set(x['output_slope'] for x in hits)),'first',hits[0],flush=True)
    save('context_head',dict(status='FINITE-EXACT',scope='All ordinary tournament types through8 from inherited complete end-rooted bank; all roots canonicalized independently; neutral response equations evaluated exactly',levels=levels,first_rules=list(all_examples.values()),increasing_slope_rules=up))
    print('ALL_RULES',sorted(all_examples));print('INCREASING',up[:3])

if __name__=='__main__':run()
