"""Exact interval kernels for two selected vertices in a C3 transitive block."""
from pathlib import Path
from itertools import product
from collections import defaultdict
from math import comb
import sys,json
import sympy as s
from tournament_tools import cyclic_blocks
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'seed_recursion_births_20261111'))
import cones

def need(x,msg):
    if not x:raise RuntimeError(msg)
def exponentials(poly):
    out=defaultdict(int)
    for visits,c in poly.items():
        for bases in product(*(range(1,j+1) for j in visits)):
            v=c
            for j,b in zip(visits,bases):v*=(-1)**(j-b)*comb(j,b)
            out[bases]+=v
    return {p:c for p,c in out.items() if c}
def evaluate(poly,args):
    ans=0
    for bases,c in poly.items():
        for p,a in zip(bases,args):c*=p**a
        ans+=c
    return ans
def symbolic(poly,names):
    ans=0
    for bases,c in poly.items():
        for p,name in zip(bases,names):
            if p>1:c*=s.Symbol(f'{p}^{name}')
        ans+=c
    return str(s.factor(ans))

def run():
    report=[];checks=0
    for orientation in [(1,2),(2,1)]:
        for present in product((0,1),repeat=3):
            ports=[];names=[];O=[];L=0
            for i,name in enumerate('abc'):
                if present[i]:ports.append(L);names.append(name);L+=1
                if i<2:O.append(L);L+=1
            rows=cyclic_blocks(L,*orientation)
            de,ue=cones.kernel(rows,ports,O);D=exponentials(de);U=exponentials(ue)
            for args in product(range(1,4),repeat=len(ports)):
                abc=[0,0,0]
                for name,a in zip(names,args):abc['abc'.index(name)]=a
                a,b,c=abc
                child,selected,_=cones.expand(rows,O,ports,args)
                h,u,w,*_=cones.pt.ib.tc.path_tables(child)
                delta=h-w[selected[0]][selected[1]];value=sum(u[i] for i in selected)
                need(evaluate(D,args)==delta and evaluate(U,args)==value,'literal interval kernel');checks+=1
                if b<=1:
                    coefficient=({0:30,1:78} if orientation==(1,2) else {0:50,1:130})[b]
                    constant=1 if orientation==(1,2) else 3
                    need(delta==1 and value==coefficient*2**a*3**c-constant*(a==0),'seed family formula')
            # For b=0 or b=1 substitute the gap variable in the complete
            # exponential kernel and compare the remaining exact polynomial.
            if 'b' not in names:gap_values=[0]
            else:gap_values=[1]
            exact=[]
            for gap in gap_values:
                collapsedD=defaultdict(int);collapsedU=defaultdict(int)
                for source,dest in [(D,collapsedD),(U,collapsedU)]:
                    for bases,value in source.items():
                        remaining=[]
                        for name,base in zip(names,bases):
                            if name=='b':value*=base**gap
                            else:remaining.append(base)
                        dest[tuple(remaining)]+=value
                collapsedD={p:c for p,c in collapsedD.items() if c}
                collapsedU={p:c for p,c in collapsedU.items() if c}
                kept=[q for q in names if q!='b']
                need(collapsedD=={(1,)*len(kept):1},'all-size unit identity')
                K=({0:30,1:78} if orientation==(1,2) else {0:50,1:130})[gap]
                expected=defaultdict(int)
                expected[tuple(2 if q=='a' else 3 for q in kept)]+=K
                if 'a' not in kept:expected[(1,)*len(kept)]-=1 if orientation==(1,2) else 3
                expected={p:c for p,c in expected.items() if c}
                need(collapsedU==expected,'all-size slope identity')
                exact.append({'gap':gap,'full_kernel_identity':True})
            report.append({'other_block_lengths':orientation,'positive_intervals':names,'defect':symbolic(D,names),'U':symbolic(U,names),'defect_terms':[[*p,c] for p,c in sorted(D.items())],'U_terms':[[*p,c] for p,c in sorted(U.items())],'unit_family_identities':exact})
            if orientation==(1,2):print(json.dumps({k:report[-1][k] for k in ('positive_intervals','defect','U')}),flush=True)
    out={'scope':'All nonnegative interval lengths a,b,c in C3[TT_(a+b+c+2),TT_1,TT_2] and its reversed small-block order; two selected vertices at positions a and a+b+1. Complete exponential kernels in every zero/positive boundary stratum.','literal_controls':checks,'records':report}
    dest=HERE/('c3_pair_family.optimized.json' if not __debug__ else 'c3_pair_family.json')
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
