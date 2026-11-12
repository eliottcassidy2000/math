"""All-length intersection certificates for synchronized birth families."""
from common import *
from birth_cones import canonical_orders
from synchrony_probe import expand
from itertools import combinations
from collections import Counter

def run():
    bank=read('synchrony_kernels.json.gz')['records'];charts=[];pairs=[];literal=Counter()
    for idx,r in enumerate(bank):
        need(r['defect_exponential'] in [[[1,1,1,1],[1,2,1,-1],[2,1,1,1]],[[1,1,1,1],[1,2,1,1],[2,1,1,-1]]],'common exact synchrony gate up to swapping the two coupled ports')
        a,b,p=r['ports'];q,o,C=quotient(r['rows'],r['selected']);canon,orders=canonical_orders(q,o)
        for order in orders:
            weights=[(1,0) if a in C[i] or b in C[i] else (0,len(C[i])) for i in order]
            charts.append({'family':idx,'quotient':list(canon),'weights':weights})
        for f in r['existing_frozen_block_tests']:
            need(f['increment_defect_minus_one'] in [[[1,4],[4,4]],[[1,12],[2,-12],[4,8]]],'positive frozen-block polynomial')
    for i,j in combinations(range(len(charts)),2):
        a,b=charts[i],charts[j]
        if a['quotient']!=b['quotient']:continue
        # Total order is8+2m in every family, so equality forces the same m.
        constraints=[(x-u,v-y) for (x,y),(u,v) in zip(a['weights'],b['weights']) if (x,y)!=(u,v)]
        impossible=any(x==0 and y!=0 for x,y in constraints)
        forced={y//x for x,y in constraints if x and y%x==0}
        impossible|=any(x and y%x!=0 for x,y in constraints) or len(forced)>1
        need(impossible or (len(forced)==1 and next(iter(forced))==1),'distinct charts intersect only at the old base m=1')
        pairs.append({'first_chart':i,'second_chart':j,'constraints':constraints,'intersection':[] if impossible else [1]})
    extkeys={tuple(r['canonical_key']) for r in read('birth_extensions_12.json.gz')['records']};overlap=[]
    for m in [1,2,3,4]:
        keys=set()
        for i,r in enumerate(bank):
            rr,O,_=expand(r['rows'],r['selected'],{r['ports'][0]:m,r['ports'][1]:m});k=key(rr,O);keys.add(k)
            if m==2 and k in extkeys:overlap.append(i)
        literal[8+2*m]=len(keys)
        need(len(keys)==(12 if m==1 else 24),'literal canonical intersection control')
    out={'scope':'All-length ordinary marked-isomorphism separation of24 synchronized birth charts at free-port length1. Same weighted quotient implies equal m by total size. Exact linear coordinate equations give all intersections.','families':24,'charts':len(charts),'same_quotient_pair_tests':len(pairs),'birth_GF_including_old_bases':'12*z^10+24*z^12/(1-z^2)','new_birth_GF':'24*z^12/(1-z^2)','literal_birth_counts':dict(literal),'order12_extension_overlap':overlap,'charts_data':charts,'pair_certificates':pairs}
    save('synchrony_grammar.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['charts_data','pair_certificates']}),flush=True)
if __name__=='__main__':run()
