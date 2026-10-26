"""Independent equations for nonsplit nodes and their exceptional conics."""
import itertools
import json
from fractions import Fraction as F


def need(ok,label):
    if not ok: raise RuntimeError(label)


def quadratic_square(a,D):
    x,y=a
    return (x*x+D*y*y,2*x*y)


def exact_controls():
    controls=[(2,3,(F(1),F(0)),(F(0),F(1))),
              (5,3,(F(1,2),F(1,2)),(F(1,2),F(-1,2))),
              (13,2,(F(1),F(0)),(F(1),F(0)))]
    for D,delta,a,b in controls:
        aa,bb=quadratic_square(a,D),quadratic_square(b,D)
        need((aa[0]+bb[0],aa[1]+bb[1])==(delta,0),'exact quadratic-field point')
    # All eight nodes: evaluate the seven original homogeneous equations in K.
    nodes=0
    for delta in [F(2),F(3),F(-1),F(5,7)]:
        d=[2/delta,1,1,1,-13,1,-1]
        L=[0,0,0,1,-1,1,-1]
        for eps,eta,theta in itertools.product([-1,1],repeat=3):
            s=[(0,0),(0,0),(0,0),(1,0),(0,eps),(eta,0),(theta,0)]
            need(all(quadratic_square(a,13)==(di*li,0) for a,di,li in zip(s,d,L)), 'original surface node')
            nodes+=1
        # Tangent equation follows from the first three rows, no fitting.
        need(d[1]==d[2]==1 and 2/d[0]==delta,'tangent coefficient')
    roots=[x for x in range(3) if (x*x-13)%3==0]
    hensel=[]
    modulus=3
    for e in range(1,10):
        need(len(roots)==2 and all(x%3 for x in roots),'two simple 3-adic roots')
        hensel.append({'exponent':e,'modulus':modulus,'roots':roots})
        roots=[x+t*modulus for x in roots for t in range(3)
               if ((x+t*modulus)**2-13)%(3*modulus)==0]
        modulus*=3
    primitive=[(a,b,c) for a,b,c in itertools.product(range(9),repeat=3)
               if (a*a+b*b-3*c*c)%9==0 and any(x%3 for x in (a,b,c))]
    need(not primitive,'primitive mod-nine obstruction')
    # Biquadratic upper-bound point 1^2+(sqrt2)^2=3.
    need(1+2==3 and all(a*a!=13 for a in range(4)),'degree-four positive control')
    return {'exact_quadratic_witnesses':[{'D':D,'delta':delta,'A':list(map(str,a)),'B':list(map(str,b))}
                for D,delta,a,b in controls], 'surface_nodes_checked':nodes,
            'hensel_roots_13':hensel,'primitive_mod9_conic3_solutions':len(primitive)}


class FF:
    # Tuple representation, independent of check_moments' encoded arithmetic.
    def __init__(self,p,degree):
        self.p,self.degree=p,degree
        self.nu=next(n for n in range(1,p) if all(x*x%p!=n for x in range(p)))
        self.el=list(itertools.product(range(p),range(p) if degree==2 else [0]))
        self.q=len(self.el)
        self.zero=(0,0); self.one=(1,0)
    def add(self,a,b): return tuple((x+y)%self.p for x,y in zip(a,b))
    def mul(self,a,b):
        x,y=a; z,w=b
        return ((x*z+self.nu*y*w)%self.p,(x*w+y*z)%self.p)
    def num(self,n): return (n%self.p,0)


def field_controls(p,degree):
    f=FF(p,degree)
    square={x:f.mul(x,x) for x in f.el}
    projective=[(f.one,a,b) for a in f.el for b in f.el]
    projective += [(f.zero,f.one,b) for b in f.el]+[(f.zero,f.zero,f.one)]
    conics={}
    for delta in f.el:
        if delta==f.zero: continue
        pts=[(a,b,c) for a,b,c in projective
             if f.add(square[a],square[b])==f.mul(delta,square[c])]
        need(len(pts)==f.q+1,'direct projective conic count')
        conics[delta]=pts
    fibre_cases=0
    for D in f.el:
        if D==f.zero: continue
        nodes=[(r,s,t) for r in f.el if square[r]==D
               for s in f.el if square[s]==f.one for t in f.el if square[t]==f.one]
        need(len(nodes) in [0,8],'actual root-node count')
        for delta,pts in conics.items():
            direct=sum(1 for node in nodes for point in pts)
            need(direct==(8*(f.q+1) if D in square.values() else 0),'complete resolved fibre')
            fibre_cases+=1
    D=f.num(13)
    need(D!=f.zero,'good prime')
    fixed_nodes=4*sum(s==D for s in square.values())
    fixed={str(delta):fixed_nodes*len(conics[f.num(delta)]) for delta in [2,3]}
    need(fixed['2']==fixed['3'],'two distinct arithmetic conics have same good-field counts')
    return {'p':p,'degree':degree,'q':f.q,'all_nonzero_D_delta_cases':fibre_cases,
            'all_conic_counts':f.q+1,'D13_node_count':fixed_nodes,'D13_resolved_counts':fixed}


if __name__=='__main__':
    result={'exact':exact_controls(),'fields':[field_controls(p,d) for p,d in
        [(5,1),(7,1),(11,1),(17,1),(19,1),(5,2),(7,2)]]}
    print(json.dumps(result,indent=2,sort_keys=True))
