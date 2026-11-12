"""Certify and separate the new odd-order synchronized birth families."""
from common import *
from birth_cones import canonical_orders
from sparse_covers import kernel as covers
from joint_nine import expo
from synchrony_probe import expand
from itertools import combinations
from collections import Counter
from math import comb
import time

def run():
    source=read('mixed_family_probe.json.gz');need(source['complete'],'complete discovery input');found={};start=time.time()
    for r in source['records']:
        q,O,C=quotient(r['rows'],r['selected']);canon,orders=canonical_orders(q,O);need(len(orders)==1,'rigid marked quotient for this declared family bank');order=orders[0]
        weights=tuple((1,len(C[i])-1) if set(C[i])&set(r['coupled_ports']) else (0,len(C[i])) for i in order)
        k=(canon,weights)
        if k in found:continue
        signs=[]
        for f in r['frozen_tests']:
            E=dict(f['increment_defect_minus_one']);R=max(E);u=[sum(comb(b,j)*c for b,c in E.items() if b>=j) for j in range(1,R+1)]
            sign=1 if u[0]>0 else -1;need(all(sign*c>=0 for c in u) and sign*u[0]>0,'all-length positive Stirling nonvanishing certificate');signs.append({'block':f['block'],'sign':sign,'Stirling_coefficients':u})
        p=r['free_ports'][0]['port'];P=r['coupled_ports']+[p];D,U=covers(r['rows'],P,[v for v in range(len(r['rows'])) if r['selected']>>v&1])
        need([[*b,c] for b,c in sorted(expo(D).items())]==r['free_ports'][0]['triple_exponential'],'independent entire mixed-family free-port kernel')
        found[k]={**r,'quotient':list(canon),'weights':weights,'frozen_nonvanishing':signs}
        if len(found)%10==0:print(json.dumps({'families':len(found),'seconds':round(time.time()-start,2)}),flush=True)
    bank=list(found.values());intersections=[]
    for i,j in combinations(range(len(bank)),2):
        a,b=bank[i],bank[j]
        if a['quotient']!=b['quotient']:continue
        eq=[(x-u,v-y) for (x,y),(u,v) in zip(a['weights'],b['weights']) if (x,y)!=(u,v)]
        need(eq,'no identical charts after quotient deduplication')
        bad=any(x==0 and y!=0 or x and y%x for x,y in eq);roots={y//x for x,y in eq if x and y%x==0}
        bad|=len(roots)>1
        root=None if bad else next(iter(roots))
        need(root is None or root<=1,'all new positive-length overlaps confined to the original base')
        intersections.append({'first':i,'second':j,'equations':eq,'intersection_m':root})
    counts={};literal=[]
    for m in [1,2,3]:
        keys=set()
        for i,r in enumerate(bank):
            rr,O,_=expand(r['rows'],r['selected'],dict.fromkeys(r['coupled_ports'],m));keys.add(key(rr,O))
            if m==2:
                h,u,w,_=tables(np.array(rr,dtype=np.int64));need(defect(h,w,[v for v in range(len(rr)) if O>>v&1])==1,'literal mixed birth unit count');literal.append({'family':i,'m':m,'H':int(h),'order':len(rr)})
        counts[11+2*m]=len(keys)
        if m>=2:need(len(keys)==len(bank),'all-length chart count literal check')
    out={'scope':'All-length strong stable-clone birth proof from complete signed free-port kernels, positive Stirling certificates for every old frozen block, and rigid weighted-quotient intersections. This is a finite family grammar, not all births.','families':len(bank),'base_order13_births':counts[13],'counts':counts,'new_GF':f'{len(bank)}*z^15/(1-z^2)','whole_GF':f'{counts[13]}*z^13+{len(bank)}*z^15/(1-z^2)','same_quotient_pair_tests':len(intersections),'intersections':intersections,'literal_controls':literal,'records':bank}
    save('mixed_grammar.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k not in ['intersections','literal_controls','records']}),flush=True)
if __name__=='__main__':run()
