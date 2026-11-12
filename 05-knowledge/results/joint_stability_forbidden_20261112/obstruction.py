"""Independent checks of the first physical mixed stable-port obstruction."""
from itertools import permutations,combinations,product
from pathlib import Path
import json,sys
import numpy as np
from census_nine import read,save,need,tables,defect,strong
from joint_nine import cones,kernel,expo,evaluate,key
from c3_pair_family import exponentials
HERE=Path(__file__).resolve().parent

def permutation_counts(rows,O):
    h=w=0;n=len(rows)
    for p in permutations(range(n)):
        bad=[i for i in range(n-1) if not rows[p[i]]>>p[i+1]&1]
        h+=not bad
        if len(bad)>1:continue
        for i in range(n-1):
            if O>>p[i]&1 and O>>p[i+1]&1 and (not bad or bad==[i]):w+=1
    return h,w

def run():
    bank=read('stable_nine.json.gz')['records'];mixed=read('joint_nine.json.gz');fail=mixed['joint_failures'];need(len(fail)==2,'two first obstructions')
    out=[]
    for ix in fail:
        r=bank[ix];rows=r['rows'];O=r['selected'];P=r['stable_ports'];sel=[i for i in range(9) if O>>i&1]
        h,w=permutation_counts(rows,O);need((h,w)==(265,264),'independent 9! permutation/cut count')
        K,U=cones.kernel(rows,P,sel);need(K=={(1,1):1,(2,2):-1},'independent split-cover skeleton')
        need(K==kernel(rows,O,P),'signed-word skeleton agreement')
        E=expo(K);UE=exponentials(U);grid=[]
        for a,b in product(range(1,5),repeat=2):
            rr,oo,_=cones.expand(rows,sel,P,[a,b]);hh,uu,ww,_=tables(np.array(rr,dtype=np.int64));dd=defect(hh,ww,oo)
            need(dd==1-(2**a-2)*(2**b-2)==evaluate(E,[a,b]),'literal full grid defect')
            slope=sum(int(uu[v]) for v in oo);need(slope==evaluate(UE,[a,b]),'literal slope grid')
            grid.append({'lengths':[a,b],'H':int(hh),'W':int(hh)-dd,'defect':dd,'U':slope})
        reverse=[((1<<9)-1)^rows[i]^(1<<i) for i in range(9)]
        out.append({'record':ix,'rows':rows,'selected':O,'ports':P,'strong':strong(rows),'permutation_H':h,'permutation_W':w,'Stirling_defect':[[*k,v] for k,v in sorted(K.items())],'U_exponentials':[[*k,v] for k,v in sorted(UE.items())],'literal_grid':grid,'converse_key':key(reverse,O)})
    need(out[0]['converse_key']==key(out[1]['rows'],out[1]['selected']),'marked converse pair')
    result={'scope':'All-length physical refutation of unrestricted individual-implies-joint stability. Smallest order among strong marked unit seeds, conditional on the inherited complete through-eight bank and the new complete order-nine census. Does not refute contractible-block joint stability.','minimum_strong_order':9,'ordinary_marked_types':2,'converse_types':1,'defect_formula':'1-(2^a-2)*(2^b-2)','exact_positive_unit_locus':'a=1 or b=1','stable_face_ideal':'<x*y>','records':out}
    save('obstruction.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run()
