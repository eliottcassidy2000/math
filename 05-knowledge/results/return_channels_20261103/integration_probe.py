"""A uniform 3-adic bulk elevator and its separate F9 phase clock."""
import json
import sympy as sp
from channels import LOOPS,SEEDS,T
from common import *

def level3(row):
    g=0
    from math import gcd
    for x in row:g=gcd(g,abs(int(x)))
    need(g>0,'nonzero bulk')
    e=0
    while g%3==0:g//=3;e+=1
    return e

def fieldmul(z,w):
    a,b=z;c,d=w
    return ((a*c+b*d)%3,(a*d+b*c+2*b*d)%3)

def power(z,k):
    out=(1,0)
    for _ in range(k):out=fieldmul(out,z)
    return out

def main():
    word=LOOPS['L2']+('b',)+LOOPS['L2']
    P=sp.Matrix(matrix(word));B=sp.Matrix([[0,1],[0,2],[1,-1]])
    Z=P*B
    Q=sp.Matrix([[Z[0,j]+Z[2,j] for j in range(2)],[Z[0,j] for j in range(2)]])
    need(B*Q==Z,'saturated quotient intertwining')
    R=Q/3
    need(all(x.q==1 for x in R) and R.det()==-100,'uniform exact level elevator')
    from math import gcd
    need(gcd(*[int(x) for x in R])==1,'normalized Smith factors one and one hundred')
    need(all((int(R[i,0])+50*int(R[i,1]))%100==0 for i in range(2)),'index-one-hundred lattice inclusion')
    return_gates=[]
    for seed in SEEDS:
        bulk=sp.Matrix([vector(seed)])*B
        residue=int(bulk[0]+50*bulk[1])%100
        need(residue!=0,'no exact D-based ternary return to these seeds')
        return_gates.append({'seed':seed,'bulk':[int(x) for x in bulk],'residue_mod100':residue})
    need(step(T,matrix(word))==tuple(-4*x for x in T),'separate error gain')
    need(sum(COSTS[l] for l in word)==19,'literal vertex clock')
    modR=R.applyfunc(lambda x:x%3)
    need(modR==sp.Matrix([[0,1],[1,2]]),'F9 multiplication matrix')
    order=next(k for k in range(1,9) if (modR**k).applyfunc(lambda x:x%3)==sp.eye(2))
    need(order==8,'primitive normalized bulk clock')
    for u in range(3):
        for v in range(3):
            if (u,v)==(0,0):continue
            orbit=[tuple(int(x)%3 for x in sp.Matrix([[u,v]])*(modR**k)) for k in range(8)]
            need(len(set(orbit))==8,'every nonzero normalized bulk has period eight')
    alpha=(1,1);eta=(0,1)
    need(alpha==power(eta,7),'actual U field phase')
    phase=[]
    for k in range(8):
        z=fieldmul(alpha,power(eta,k));f=power(z,3)
        need(f==fieldmul(alpha,power(eta,(3*k+6)%8)),'source-relative Frobenius phase')
        phase.append({'k':k,'alpha':z,'frobenius_index':(3*k+6)%8})
    X=sp.symbols('X');need(sp.expand(P.charpoly(X).as_expr()-(X+4)*(X**2-1093533*X-900))==0,'exact cubic characteristic')
    records=[]
    for k in range(4):
        sources=[]
        for i,seed in enumerate(SEEDS):
            w=seed+word*k;v=vector(w);rows=actual_target(w)
            h,width,peak=literal_frontier(rows)
            need(h==v[0],'literal new return-family count')
            beta=sp.Matrix([v])*B
            need(level3(beta)==k+(i>0),'exact bulk level at actual sources')
            sources.append({'seed':i,'word':w,'vertices':len(rows),'port':len(rows)-1,
                'state':v,'rows':rows,'bulk':[int(x) for x in beta],'level3':level3(beta),
                'backward_width':width,'frontier_peak':peak})
        m=sources[0]['state'][0]
        need(sources[1]['state'][0]==3*m+2*(-4)**k and sources[2]['state'][0]==3*m+4*(-4)**k,'transported ternary companions')
        records.append({'k':k,'sources':sources})
    h=[r['sources'][0]['state'][0] for r in records]
    need(h[3]==1093529*h[2]+4375032*h[1]+3600*h[0],'cubic H recurrence')
    result={'status':'PASS','return_word':word,'gain':-4,'vertex_cost':19,
        'matrix':[[int(x) for x in row] for row in P.tolist()],
        'bulk_basis':B.tolist(),'bulk_matrix':Q.tolist(),'normalized_bulk_matrix':R.tolist(),
        'normalized_determinant':int(R.det()),'normalized_mod3_order':order,
        'normalized_Smith_factors':[1,100],'exact_return_lattice_gates':return_gates,
        'phase':phase,'families':records,'H_values':h,
        'scope':'Uniform level increase for every nonzero bulk; an actual infinite family, not a closed ternary source bank.'}
    # Convert SymPy integers explicitly for portable JSON.
    (HERE/'integration.json').write_text(json.dumps(result,indent=2,default=int)+'\n',encoding='utf8')
    print(json.dumps({'status':'PASS','H_values':h,'gain':-4,'bulk_level_step':1,
        'normalized_clock_period':8,'literal_max_order':max(s['vertices'] for r in records for s in r['sources'])}))
if __name__=='__main__':main()
