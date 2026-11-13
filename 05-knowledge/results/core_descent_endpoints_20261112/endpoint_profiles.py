"""Whole one-clock endpoint observers for the70 inherited birth families.

Queries are literal graphs; updates are the independently audited universal
two-live-port observer. Matching H alone is reported as such, never as seed
stability or closure of every lower head.
"""
from context import *
from fractions import Fraction
from collections import Counter,defaultdict
from math import gcd,isqrt
import time

OBS=load('two_port_envelope.json',HERE.parent/'adaptive_state_induction_20261103')
QUERIES=OBS['rank']['suffixes']
MATRICES=[[(i,j,int(2*Fraction(x))) for i,row in enumerate(M) for j,x in enumerate(row) if Fraction(x)] for M in OBS['rank']['query_matrices']]

def attach(rows,roots,mask,drop=2):
    n=len(rows);incoming={roots[i] for i in range(2) if mask>>i&1}
    out=[r|((1<<n) if i in incoming else 0) for i,r in enumerate(rows)]
    out.append(sum(1<<i for i in range(n) if i not in incoming))
    roots=list(roots)
    if drop<2:roots[drop]=n
    return out,roots

def observer(rows,roots,clock):
    profiles=[]
    for word in QUERIES:
        rr,pp=list(rows),list(roots)
        for mask in word:rr,pp=attach(rr,pp,mask)
        profiles.append(list(trim(expcoeff([int(x) for x in skeleton(np.array(rr,dtype=np.int64),clock)]))))
    width=max(map(len,profiles))
    return tuple(tuple(p+[0]*(width-len(p))) for p in profiles)

def advance(X,letter):
    width=len(X[0]);Y=[[0]*width for _ in range(17)]
    for i,j,c in MATRICES[letter]:
        for k,x in enumerate(X[i]):Y[j][k]+=c*x
    need(all(x%2==0 for p in Y for x in p),'integral physical exponential observer')
    return tuple(tuple(x//2 for x in p) for p in Y)

def sources():
    for i,r in enumerate(load('birth_sources_synchrony_12.json.gz')['records']):
        blocks=defaultdict(list)
        for v,label in enumerate(r['labels']):blocks[label].append(v)
        roots=[b[0] for b in blocks.values() if len(b)>1]
        need(len(roots)==2,'two synchronized blocks')
        for p in r['stable_ports']:yield {'family':f'even:{i}','rows':r['rows'],'roots':roots,'clock':p}
    for i,r in enumerate(load('mixed_grammar.json.gz')['records']):
        for f in r['free_ports']:yield {'family':f'odd:{i}','rows':r['rows'],'roots':r['coupled_ports'],'clock':f['port']}

def split_sources():
    # Name the two end vertices of a length(a+2) free module and keep its
    # middle TT_a as the clock. Future cuts can distinguish these endpoints.
    for src in sources():
        p=src['clock'];labels=[v for v in range(len(src['rows'])) for _ in range(3 if v==p else 1)]
        rr=[sum(1<<j for j,w in enumerate(labels) if (v==w and i<j) or (v!=w and src['rows'][v]>>w&1)) for i,v in enumerate(labels)]
        cc=[i for i,v in enumerate(labels) if v==p]
        yield {**src,'rows':rr,'roots':[cc[0],cc[2]],'clock':cc[1],'clock_offset':2}

def compile_sources(split=False):
    start=time.time();out=[];stats=Counter();checks=0
    for src in (split_sources() if split else sources()):
        need(src['clock'] not in src['roots'],'clock never named')
        X=observer(src['rows'],src['roots'],src['clock'])
        # Two physical updates, including a root replacement and two-in cut.
        for letter in [3,10]:
            rr,pp=attach(src['rows'],src['roots'],*divmod(letter,3))
            y=advance(X,letter)
            for k in [0,1,4]:
                rrr,ppp=rr,pp
                for mask in QUERIES[k]:rrr,ppp=attach(rrr,ppp,mask)
                e=trim(expcoeff([int(x) for x in skeleton(np.array(rrr,dtype=np.int64),src['clock'])]))
                need(trim(list(y[k]))==e,'literal whole-profile continuation control');checks+=1
        stats[len(trim(list(X[0])))]+=1
        out.append({**src,'observer':X})
        if len(out)%10==0:print('compiled',len(out),'seconds',round(time.time()-start,1),flush=True)
    result={'scope':'Every free clock of each of the70 inherited birth charts at its stored smallest source; two named roots chosen in the coupled blocks. Complete17-query exponential observers, all positive clock lengths; no arbitrary-size census.','source_observers':len(out),'H_dominant_base_distribution':dict(stats),'literal_profile_controls':checks,'records':out}
    result['split_clock_endpoints']=split
    save('split_endpoint_profiles.json.gz' if split else 'endpoint_profiles.json.gz',result)
    print({k:v for k,v in result.items() if k!='records'},flush=True)

def target_queries():
    factors={}
    for r in load('hamilton_clocks_9.json.gz')['records']:
        h=value(r['exponential_coefficients'],1)
        factors.setdefault(h,r['rows'])
    queries=defaultdict(list)
    for r in load('even_valuations.json.gz')['records']:
        for d in range(9):
            for shift in range(13):
                t=[Fraction(c*(i+1)**shift,2**d) for i,c in enumerate(r['exponential_coefficients'])];t[0]+=1
                if any(c.denominator!=1 for c in t):continue
                t=list(map(int,t));G=abs(gcd(*t));divs=set()
                for q in range(1,isqrt(G)+1):
                    if G%q==0:divs.update([q,G//q])
                for factor in divs:
                    if factor in factors:queries[trim([c//factor for c in t])].append({'birth_id':r['birth_id'],'port':r['port'],'d':d,'shift':shift,'factor':factor,'factor_rows':factors[factor]})
    return queries

def search(depth=4,split=False):
    queries=target_queries();maxh=max(value(p,1) for p in queries);maxbase=max(map(len,queries))
    records=load('split_endpoint_profiles.json.gz' if split else 'endpoint_profiles.json.gz',HERE)['records'];hits=[];stats=Counter();seen=set();start=time.time()
    # All12 commands, including neutral universal-source insertion. Depth is
    # explicit; source/word witnesses survive observer deduplication.
    frontier=[]
    for r in records:
        X=tuple(map(tuple,r['observer']))
        if len(trim(list(X[0])))>maxbase:stats['sources_pruned_by_growth']+=1;continue
        frontier.append((r,(),X))
    for level in range(depth+1):
        nxt=[]
        for src,word,X in frontier:
            key=X
            if key in seen:stats['duplicate_observers']+=1;continue
            seen.add(key);stats[f'states_depth{level}']+=1
            p=trim(list(X[0]));h=value(p,1)
            if p in queries:
                for q in queries[p]:hits.append({**q,'source':{k:v for k,v in src.items() if k!='observer'},'word':word,'H_profile':p,'full_observer':X})
            if level==depth:continue
            for a in range(12):
                Y=advance(X,a);hp=trim(list(Y[0]));hy=value(hp,1)
                need(hy>=h,'insertion monotonicity')
                if hy>maxh or len(hp)>maxbase:stats['pruned_children']+=1;continue
                nxt.append((src,word+(a,),Y))
        frontier=nxt
        print('depth',level,'seen',len(seen),'next',len(nxt),'hits',len(hits),'seconds',round(time.time()-start,1),flush=True)
    out={'scope':f'All12 universal two-live-port commands through word length{depth}, starting at the declared70-chart free-clock source observers. Matching exact entire H profiles with d0..8,shift0..12 and attained content factors. Observer equality is exact only for this continuation language.','statistics':dict(stats),'hits':hits,'matched_axes':sorted({(h['birth_id'],h['port']) for h in hits})}
    out['split_clock_endpoints']=split
    save(f'{"split_" if split else ""}endpoint_search_{depth}.json.gz',out)
    print({k:v for k,v in out.items() if k!='hits'},flush=True)

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='search':search(int(sys.argv[2]) if len(sys.argv)>2 else 4,'split' in sys.argv)
    else:compile_sources('split' in sys.argv)
