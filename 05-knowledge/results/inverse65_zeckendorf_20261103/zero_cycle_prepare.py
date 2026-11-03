from pathlib import Path
from fractions import Fraction
import json
import sympy as sp
from bank import *
HERE=Path(__file__).resolve().parent
units=[]
for l,m in MATS.items():
    det=int(sp.Matrix(m).det())
    if abs(det)!=1:continue
    units.append((l,det,m))
    need(m[0][0]%2==1 and m[1][0]%2==0 and m[2][0]%2==0,'mod2 triangular block')
    lower=tuple(tuple(m[i][j]%2 for j in (1,2)) for i in (1,2))
    need(lower in (((1,0),(0,1)),((0,1),(1,0))),'mod2 unipotent lower block')
need([l for l,d,m in units]==['a','b','q0','r0','r20','r30','r6'],'unimodular bank')
need(max(sum(row) for l,d,m in units for row in m)**7<2**63,'C++ matrix arithmetic range')
need(9**3>5*5**3 and 9**5>18*5**5 and 9**8>77*5**8,'exact length bounds')
(HERE/'units.txt').write_text('\n'.join(l+' '+str(d)+' '+' '.join(str(x) for row in m for x in row) for l,d,m in units)+'\n')
weights={l:min(Fraction(sum(m[i][j]*(6,9,5)[j] for j in range(3)),(6,9,5)[i]) for i in range(3)) for l,d,m in units}
certificates={}
for q,upper,trace,second,det in ((3,5,5,3,-1),(6,18,19,19,1),(9,77,77,75,-1)):
    queue=[((),Fraction(1),sp.eye(3))];records=[]
    for word,weight,matrix in queue:
        for l,sg,rows in units:
            nw=weight*weights[l]
            if nw>=upper:continue
            w=word+(l,);mm=matrix*sp.Matrix(rows);queue.append((w,nw,mm))
            coeff=mm.charpoly().all_coeffs();got=(-int(coeff[1]),int(coeff[2]),-int(coeff[3]))
            need(got!=(trace,second,det),('possible zero clock',q,w))
            records.append({'word':w,'norm_lower_bound':str(nw),'trace_second_det':got})
    certificates[str(q)]=records
result={'status':'PASS','unit_letters':[l for l,d,m in units],
        'norm_weights':{l:str(w) for l,w in weights.items()},
        'zero_period_candidates':[3,6,9],'maximum_emitted_lengths':[2,4,7],
        'complete_pruned_candidates':certificates}
(HERE/'zero_cycles.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','unit_letters':result['unit_letters'],'candidate_counts':{q:len(v) for q,v in certificates.items()}}))
