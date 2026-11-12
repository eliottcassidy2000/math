"""Compile the finite head obligations of all binary birth recurrences.

Search actual(D,oddpart(U)) interfaces on the complete <=9 catalogue and
the declared order10 source subbank. Keep a best actual graph for every
query and dyadic shift. No scalar residue is accepted as a graph witness.
"""
from pathlib import Path
from collections import Counter
from math import gcd
import gzip,json,time
import numpy as np
from numba import njit,types
from numba.typed import Dict
from census_nine import read,save,need,tables
from head_repair import graphs
HERE=Path(__file__).resolve().parent

def divisors(n):return [d for d in range(1,n+1) if n%d==0]

@njit(cache=True)
def scan(H,U,W,queries):
    n=len(U);N=1<<n;slopes=np.zeros(N,dtype=np.int64);taxes=np.zeros(N,dtype=np.int64);found=[]
    for m in range(1,N):
        b=m&-m;v=0
        while (1<<v)!=b:v+=1
        rest=m^b;slopes[m]=slopes[rest]+U[v];tax=taxes[rest]
        for w in range(v+1,n):
            if rest>>w&1:tax+=W[v,w]
        taxes[m]=tax;d=H-tax;u=slopes[m];odd=u;j=0
        while odd>0 and odd%2==0:odd//=2;j+=1
        pair=(d,odd)
        if pair in queries:found.append((m,queries[pair],u,j))
    return found

def run():
    families={}
    for r in read('affine_births.json.gz')['binary_families']:
        pair=(r['odd_coefficient_A'],r['odd_coefficient_B'])
        if pair not in families or r['valuation_after_start']<families[pair]['valuation_after_start']:families[pair]=r
    factors={1:[]}
    with gzip.open(HERE.parent/'metagraph_aspect_sequences_20261111/records.jsonl.gz','rt') as f:
        for line in f:
            r=json.loads(line)
            if r['H'] not in factors:factors[r['H']]=r['rows']
    query_set=set();obligations=[]
    for (A,B),r in sorted(families.items()):
        for h in r['unpaid_heads']:
            s=h['s'];D=h['offset'];options=[]
            for G in divisors(gcd(abs(D),A)):
                if G not in factors:continue
                a=A//G;v=(a&-a).bit_length()-1;pair=(D//G,a>>v);query_set.add(pair)
                options.append({'factor':G,'query':pair,'exponent_offset':s+v})
            obligations.append({'A':A,'B':B,'s':s,'offset':D,'coefficient':A*2**s,'options':options})
    pairs=sorted(query_set);qdict=Dict.empty(key_type=types.UniTuple(types.int64,2),value_type=types.int64)
    for i,pair in enumerate(pairs):qdict[pair]=i
    best={};stats=Counter();hits=0;start=time.time()
    for rows in graphs():
        n=len(rows);stats[n]+=1;h,u,w,_=tables(np.array(rows,dtype=np.int64))
        for m,qi,slope,j in scan(h,u,w,qdict):
            hits+=1;k=(qi,j);item={'query_index':qi,'rows':rows,'selected':m,'H':int(h),'defect':pairs[qi][0],'U':slope,'dyadic_shift':j}
            if k not in best or (n,rows,m)<(len(best[k]['rows']),best[k]['rows'],best[k]['selected']):best[k]=item
        if sum(stats.values())%100000==0:print(json.dumps({'source_counts':dict(stats),'matched_presentations':hits,'retained_interfaces':len(best),'seconds':round(time.time()-start,1)}),flush=True)
    need(stats==Counter({1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,10:465023}),'complete declared source universe')
    least={}
    for line in (HERE.parent/'path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt').read_text().splitlines():
        if line.startswith('#'):continue
        f=line.split();least[int(f[0])]=int(f[1]) if f[1].isdigit() else None
    resolved=[];unresolved=[];exceptions=[];used_factors=set();used_interfaces=set()
    for o in obligations:
        choices=[]
        for option in o['options']:
            qi=pairs.index(option['query']);G=option['factor'];off=option['exponent_offset']
            for (q,j),r in best.items():
                if q!=qi:continue
                startk=max(1,j-off+1);heads=[];unknown=[]
                for k in range(1,startk):
                    N=o['offset']+o['coefficient']*2**k
                    if N in [7,21]:heads.append({'k':k,'H':N,'status':'globally forbidden'})
                    elif least.get(N) is not None:heads.append({'k':k,'H':N,'status':'inherited complete order11','least_order':least[N]})
                    else:unknown.append({'k':k,'H':N})
                # No broad finite-prefix result is silently used here.
                cost=len(r['rows'])+1+off-j+len(factors[G])
                choices.append((len(unknown),startk,cost,qi,j,G,heads,unknown))
        if not choices:unresolved.append({**o,'reason':'no physical endpoint interface in the declared bank'});continue
        unknown_count,startk,cost,qi,j,G,heads,unknown=min(choices,key=lambda x:x[:6])
        answer={**o,'factor':G,'query_index':qi,'dyadic_shift':j,'tail_from_parameter':startk,'output_order_constant':cost,'finite_heads':heads}
        if unknown:unresolved.append({**answer,'reason':'finite heads outside the inherited order11 support','unverified_heads':unknown});continue
        resolved.append(answer);used_factors.add(G);used_interfaces.add((qi,j));exceptions.extend(h for h in heads if h['status']=='globally forbidden')
    unresolved_laws={(o['A'],o['B']) for o in unresolved};closed=[pair for pair in sorted(families) if pair not in unresolved_laws]
    # Recount every actual chosen graph, including factor witnesses.
    for qi,j in used_interfaces:
        r=best[qi,j];h,u,w,_=tables(np.array(r['rows'],dtype=np.int64));O=[i for i in range(len(r['rows'])) if r['selected']>>i&1]
        need(int(h)-sum(int(w[a,b]) for ia,a in enumerate(O) for b in O[ia+1:])==r['defect'] and sum(int(u[a]) for a in O)==r['U'],'retained physical interface')
    for G in used_factors:
        if G>1:need(int(tables(np.array(factors[G],dtype=np.int64))[0])==G,'literal factor witness')
    result={'scope':'All540 distinct binary coefficient laws from the complete birth bank through10, using the minimum available tail valuation for each law. Every lower-head obligation compiled into actual affine endpoint interfaces and actual ordered-product factors. The searched source universe is all ordinary orders<=9 plus the strong order10 sources with a unit selection. Finite head exceptions use only the inherited complete order11 spectrum. Missing compiler interfaces do not imply impossible targets.',
            'source_counts':dict(stats),'distinct_laws':len(families),'head_obligations':len(obligations),'queries':len(pairs),'matched_presentations':hits,'retained_interfaces':len(best),'resolved_head_obligations':len(resolved),'unresolved_head_obligations':len(unresolved),'fully_closed_laws':len(closed),'closed_laws':closed,'unresolved_laws':sorted(unresolved_laws),'resolved':resolved,'unresolved':unresolved,'queries_by_index':pairs,'interfaces':[best[k] for k in sorted(best)],'factor_rows':{str(g):factors[g] for g in sorted(used_factors)},'families':[families[p] for p in sorted(families)]}
    save('bulk_head_interfaces.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k in ['source_counts','distinct_laws','head_obligations','queries','matched_presentations','retained_interfaces','resolved_head_obligations','unresolved_head_obligations','fully_closed_laws']}),flush=True)
if __name__=='__main__':run()
