"""Whole seed tensors, valuation-aligned descent and cross-core relations."""
from context import *
from endpoint_profiles import sources
from sparse_covers import kernel
from joint_nine import expo
from collections import defaultdict,Counter
from math import gcd,comb
from fractions import Fraction
from itertools import combinations
import time

def diagonal(E):
    out=defaultdict(int)
    for (a,b,z),c in E.items():out[a*b,z]+=c
    return {k:v for k,v in out.items() if v}

def run():
    even=load('birth_sources_synchrony_12.json.gz')['records'];odd=load('mixed_grammar.json.gz')['records']
    records=[];unclassified=[];start=time.time();classes=Counter()
    for src in sources():
        typ,idx=src['family'].split(':');original=(even if typ=='even' else odd)[int(idx)];O=original['selected']
        ports=src['roots']+[src['clock']]
        need(strong(src['rows']),'physical strong source')
        need(O>0 and all(not O>>v&1 for v in ports),'nonempty selected gate and unselected variable ports')
        # Acyclic variable ports make the finite run polynomial exact.
        degrees=[sum(src['rows'][p]>>q&1 for q in ports) for p in ports]
        need(sorted(degrees)==[0,1,2],'acyclic three-port run kernel')
        D,U=kernel(src['rows'],ports,[v for v in range(len(src['rows'])) if O>>v&1])
        Ed=diagonal(expo(D));Eu=diagonal(expo(U))
        need(Ed=={(1,1):1},'entire synchronized unit defect including every free length')
        if not all(z in [1,2] for b,z in Eu):
            unclassified.append({**src,'selected':O,'reason':'higher free-clock mode','diagonal_U':[[*b,c] for b,c in sorted(Eu.items())],'defect_runs':[[*b,c] for b,c in sorted(D.items())],'slope_runs':[[*b,c] for b,c in sorted(U.items())]});continue
        A={b:c for (b,z),c in Eu.items() if z==1};B={b:c for (b,z),c in Eu.items() if z==2}
        if not(len(A)==1 and 4 in A and A[4]>0):
            unclassified.append({**src,'selected':O,'A':sorted(A.items()),'B':sorted(B.items()),'defect_runs':[[*b,c] for b,c in sorted(D.items())],'slope_runs':[[*b,c] for b,c in sorted(U.items())]});continue
        vA=v2(A[4]);anchor=A[4]>>vA
        vB=v2(abs(B[1]));need(all(v2(abs(c))+v2(b)>vB for b,c in B.items() if b!=1),'constant controls v2(B(m)) for m>=1')
        need(all(c%(1<<vB)==0 for c in B.values()),'coefficientwise exact normalization before binary shifts')
        need(3+vA-vB>0,'minimum positive free clock in the guarded chart')
        bnorm={b:c>>vB for b,c in B.items()};content=abs(gcd(*bnorm.values()));primitive={b:c//content for b,c in bnorm.items()}
        monotone=[sum(comb(b,j)*c for b,c in bnorm.items() if b>=j) for j in range(1,max(bnorm)+1)]
        need(monotone[0]>0 and all(c>=0 for c in monotone) and monotone[1]>0,'positive increasing odd B via Stirling basis')
        rec={**src,'selected':O,'defect_runs':[[*b,c] for b,c in sorted(D.items())],'slope_runs':[[*b,c] for b,c in sorted(U.items())],'A_coefficient':A[4],'A_valuation_offset':vA,'anchor':anchor,'B':sorted(B.items()),'B_valuation':vB,'B_odd':sorted(bnorm.items()),'B_odd_content':content,'primitive_B_odd':sorted(primitive.items()),'B_odd_Stirling':monotone,'z_offset':vA-vB}
        records.append(rec);classes[anchor,vA,vB,tuple(sorted(primitive.items()))]+=1
    links=[]
    # Match the inherited actual extension parent, not just proportional data.
    for child in records:
        if not child['family'].startswith('odd:'):continue
        orig=odd[int(child['family'].split(':')[1])];parent_name=f"even:{orig['parent_synchrony_family']}"
        for parent in records:
            if parent['family']!=parent_name or parent['clock']!=child['clock']:continue
            if parent['primitive_B_odd']!=child['primitive_B_odd'] or parent['anchor']!=child['anchor']:continue
            n=len(parent['rows']);mask=(1<<n)-1
            need(len(child['rows'])==n+1 and [r&mask for r in child['rows'][:n]]==parent['rows'],'literal core extension')
            need(child['selected']==parent['selected'],'selected gate retained by core extension')
            if parent['roots']!=child['roots']:continue
            ratio=Fraction(child['B_odd_content'],parent['B_odd_content'])
            if ratio<=1:continue
            a,b=ratio.numerator,ratio.denominator
            links.append({'parent_family':parent['family'],'child_family':child['family'],'clock':parent['clock'],'incoming_new': [v for v in range(n) if child['rows'][v]>>n&1],'outgoing_new':[v for v in range(n) if child['rows'][n]>>v&1],'child_z_minus_parent_z':child['z_offset']-parent['z_offset'],'ratio':[a,b],'anchor':parent['anchor'],'identity':f'{a}*parent_c-{b}*child_c={(a-b)*parent["anchor"]}'})
    out={'scope':'All78 free-port presentations of the70 stored infinite birth charts; exact three-port D,U run polynomials. A single base4 constant-mode subfamily receives the displayed valuation-aligned descent. Other profiles are retained separately. No assertion that these coefficient families cover every odd integer.','records':records,'other_profiles':unclassified,'cross_core_links':links,'statistics':{'presentations':len(records)+len(unclassified),'base4_descent_presentations':len(records),'other_presentations':len(unclassified),'distinct_tail_shapes':len(classes),'verified_parent_child_links':len(links)},'tail_shapes':[{'anchor':k[0],'A_valuation_offset':k[1],'B_valuation':k[2],'primitive_B_odd':k[3],'multiplicity':v} for k,v in classes.items()]}
    save('synchronized_descent.json.gz',out)
    print(out['statistics'],'seconds',round(time.time()-start,2),flush=True)
    print('anchors',Counter(r['anchor'] for r in records),'ratios',Counter(tuple(r['ratio']) for r in links),flush=True)
    print('links',links[:4],flush=True)
if __name__=='__main__':run()
