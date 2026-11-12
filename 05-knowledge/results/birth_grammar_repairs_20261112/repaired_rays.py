"""Actual adjacency constructors closing five formerly even higher-mode families."""
from common import *
from construct import tail
from hamilton_clocks import skeleton,trim
import argparse

def setup():
    even={(r['birth_id'],r['port']):r for r in read('even_valuations.json.gz')['records']}
    old=read('odd_axis_bank.json.gz');direct=read('direct_head_profiles.json.gz')['hits']
    cases=[]
    for bid,p,shift in [(122,2,0),(334,1,0),(741,3,1)]:
        match=next(r for r in old['matches'] if (r['birth_id'],r['port'],r['original_parameter_shift'])==(bid,p,shift))
        axis=old['records'][match['odd_axis_index']];source=even[bid,p]
        need(all(c*(i+1)**shift==2*axis['exponential_coefficients'][i] for i,c in enumerate(trim(source['exponential_coefficients'].copy()))),'whole normalized odd-axis identity')
        cases.append({'family':bid,'source':source,'method':'odd_seed','shift':shift,'repair':axis,'start':1+shift})
    for bid,p,shift in [(14,1,1),(42,5,0)]:
        match=next(r for r in direct if (r['birth_id'],r['port'],r['shift'])==(bid,p,shift))
        need(match['division_exponent']==0 and match['factor']==1,'whole lowest-head identity')
        source=even[bid,p];cases.append({'family':bid,'source':source,'method':'Hamilton_head','shift':shift,'repair':match['profile'],'start':max(1+shift,source['valuation_start'])})
    band={r['c']:r for r in read('coefficient_band.json.gz',True)['records']};small={}
    for r in read('hamilton_clocks_9.json.gz')['records']:
        h=value(r['exponential_coefficients'],1)
        if h not in small:small[h]=r['rows']
    prefix=[]
    for r in cases:
        for a in range(1,r['start']):
            U=value(r['source']['exponential_coefficients'],a);c=U//2**v2(U);seed=band[c]
            for s in range(1,seed['v2U']+1):need(1+c*2**s in small,'literal finite-prefix head witness exists through9')
            prefix.append({'family':r['family'],'a':a,'c':c,'seed':seed,'heads':[{'s':s,'rows':small[1+c*2**s]} for s in range(1,seed['v2U']+1)]})
    return cases,prefix

def construct(case,prefix,a,s):
    need(a>=1 and s>=1,'positive parameters');src=case['source'];U=value(src['exponential_coefficients'],a);j=v2(U);c=U//2**j
    if a<case['start']:
        rec=next(r for r in prefix if r['family']==case['family'] and r['a']==a);r=rec['seed'];d=r['v2U']
        rows=tail(r['rows'],[v for v in range(len(r['rows'])) if r['selected']>>v&1],s-d) if s>d else next(h['rows'] for h in rec['heads'] if h['s']==s)
        method='finite-prefix coefficient-band witness'
    elif case['method']=='odd_seed':
        r=case['repair'];labs,rr=ib.inflate(r['rows'],r['port'],a-case['shift']);O=[v for v,q in enumerate(labs) if r['selected']>>q&1];rows=tail(rr,O,s);method='normalized odd slope from another core'
    elif s==1:
        r=case['repair'];_,rows=ib.inflate(r['rows'],r['port'],a-case['shift']);method='separate complete Hamilton profile'
    else:
        need(j==1,'declared eventual valuation');labs,rr=ib.inflate(src['rows'],src['port'],a);O=[v for v,q in enumerate(labs) if src['selected']>>q&1];rows=tail(rr,O,s-j);method='original unit tail'
    return {'family':case['family'],'a':a,'s':s,'c':c,'H':1+c*2**s,'order':len(rows),'method':method,'rows':rows}

def audit():
    cases,prefix=setup();controls=[]
    for r in cases:
        p=r['repair']
        if r['method']=='Hamilton_head':
            actual=[int(x) for x in skeleton(np.array(p['rows'],dtype=np.int64),p['port'])];need(trim(expcoeff(actual))==tuple(p['exponential_coefficients']),'full physical Hamilton kernel replay')
        else:
            from sparse_covers import kernel
            D,U=kernel(p['rows'],[p['port']],[v for v in range(len(p['rows'])) if p['selected']>>v&1]);need(D=={(1,):1},'full physical repair stability');need(trim(expcoeff([U.get((i,),0) for i in range(1,max(k[0] for k in U)+1)]))==tuple(p['exponential_coefficients']),'independent full odd-slope kernel')
        for a in range(1,7):
            for s in range(1,7):
                got=construct(r,prefix,a,s)
                if got['order']>18:continue
                h,_,_,_=tables(np.array(got['rows'],dtype=np.int64));need(int(h)==got['H'],'whole constructed graph HP count');controls.append({k:v for k,v in got.items() if k!='rows'})
    out={'scope':'All positive a,s for five original even higher-mode axes. Physical adjacency constructors, complete symbolic profile identities and finite exceptional parameter witnesses; no claim of all478 closures.','closed_families':5,'families':cases,'finite_prefix_certificates':prefix,'literal_whole_graph_controls':controls,'partially_closed_birth4':'The head s=1 is H of the6-vertex pointed profile[-1,6,4]; s=2 remains unpaid by this construction.'}
    save('repaired_rays.json.gz',out);print(json.dumps({'closed_families':5,'finite_prefixes':[{k:v for k,v in r.items() if k not in ['seed','heads']} for r in prefix],'literal_controls':len(controls)}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--family',type=int);p.add_argument('--parameter',type=int,default=1);p.add_argument('--exponent',type=int,default=1);a=p.parse_args()
    if a.family:
        cases,prefix=setup();print(json.dumps(construct(next(r for r in cases if r['family']==a.family),prefix,a.parameter,a.exponent),indent=2))
    else:audit()
