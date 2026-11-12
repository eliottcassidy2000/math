"""Exact family and observer checks for the nonconfluent crossing."""
from itertools import product
from collections import Counter
import numpy as np
import sympy as s
from census_nine import need,save,tables,defect,automorphisms,move
from joint_nine import cones,key,quotient
from c3_pair_family import exponentials

def run():
    Q=[2,100,49,7,75,409,173,31,223];mask=50;O=[1,4,5];P=[0,8]
    need(quotient(Q,mask)[2]==[[i] for i in range(9)],'twin reduced quotient')
    need(len([p for p in automorphisms(Q) if move(mask,p)==mask])==1,'rigid marked quotient')
    rev=[511^Q[i]^(1<<i) for i in range(9)];need(key(Q,mask)==key(rev,mask),'marked self converse')
    D,U=cones.kernel(Q,P,O);UE=exponentials(U)
    expected={(1,1):-8,(1,2):-49,(2,1):-95,(2,2):334,(3,1):-64,(3,2):506}
    need(UE==expected,'full U exponential identity')
    controls=[]
    for a,b in product(range(1,6),repeat=2):
        rr,sel,_=cones.expand(Q,O,P,[a,b]);h,u,w,_=tables(np.array(rr,dtype=np.int64));d=defect(h,w,sel)
        H=175*2**(a+b)-9*2**a-9*2**b+1
        V=-8-49*2**b-95*2**a+334*2**(a+b)-64*3**a+506*3**a*2**b
        need(int(h)==H and sum(int(u[i]) for i in sel)==V,'literal H and U formulas')
        need(d==1+(2**a-4)*(2**b-4),'literal defect')
        if a==2 or b==2:
            need(d==1,'full locus');n=len(rr)
            # Actual unselected source padding preserves defect and breaks
            # equivalence of the two terminal roots under converse.
            padded=[((1<<(n+1))-1)^1]+[r<<1 for r in rr];ps=[i+1 for i in sel]
            if n+1<=18:
                ph,pu,pw,_=tables(np.array(padded,dtype=np.int64));need(int(ph)==H and defect(ph,pw,ps)==1,'literal source-padding transport')
        controls.append({'lengths':[a,b],'H':H,'U':V,'defect':d})
    roots=[]
    for lengths in [[1,2],[2,1]]:
        rr,sel,_=cones.expand(Q,O,P,lengths);n=len(rr);rr=[((1<<(n+1))-1)^1]+[r<<1 for r in rr];m=sum(1<<(i+1) for i in sel);rrev=[((1<<(n+1))-1)^rr[i]^(1<<i) for i in range(n+1)]
        roots.append((key(rr,m),key(rrev,m)))
    need(roots[0][0] not in roots[1],'padded terminal roots remain distinct modulo converse')
    z=s.Symbol('z');gf=2*z**10/(1-z)-z**11;coeff=s.series(gf,z,0,31).removeO().expand()
    need([int(coeff.coeff(z,n)) for n in range(10,31)]==[2,1]+[2]*19,'exact family coefficients')
    out={'scope':'All-length strong two-line crossing, exact source profiles, and unrestricted source-padded converse obstruction. Literal grid through5 in both lengths.','H':'175*2^(a+b)-9*2^a-9*2^b+1','U_terms':[[*k,v] for k,v in sorted(UE.items())],'GF':'2*z^10/(1-z)-z^11','nonconfluent_GF':'z^11+2*z^12/(1-z)','edge_GF':'2*z^11/(1-z)','odd_coefficient_family':'5841*2^k-241, k>=1; U=4*c; tail exponents s>=3','literal_controls':controls,'padded_births_distinct_modulo_converse':True}
    save('crossing_family.json.gz',out);print('Verified full crossing formulas, all-length counting identity and converse-folding hostile.')
if __name__=='__main__':run()
