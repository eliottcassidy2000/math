"""Integral pairing geometry, affine source guard, and dyadic local controller."""
from itertools import combinations,product
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from mechanisms import *
HERE=Path(__file__).resolve().parent

def radical(m):
    need(m%2==0 and m>=2,'even source length')
    l=m//2
    return [l-1]+[v for j in range(l) for v in (j-1,j+1)]

def peel(m):
    B=s.Matrix(wedge(m));N=m+1;p=s.eye(N)[:,0];q=s.eye(N)[:,-1]
    vectors=[p,q]
    for i in range(1,N-1):
        v=s.eye(N)[:,i]
        vectors.append(v+(q.T*B*v)[0]*p-(p.T*B*v)[0]*q)
    vectors=vectors[:2]+[vectors[-1]]+vectors[2:-1]
    return s.Matrix.hstack(*vectors)

def streamed(w,u,e):
    need(e>=1,'positive precision');M=1<<e;tail=[];acc=0;binary=0;last=0
    for a,b in zip(w,u):
        for distance,(c,d) in enumerate(reversed(tail),1):
            if distance>=2:acc=(acc+(1<<(distance-2))*(c*b-a*d))%M
        tail.append((a,b));tail=tail[-(e+1):]
        last=b-a;binary=(2*binary+last)%(2*M)
    need((binary+last)%2==0,'final exact half')
    return (acc+(binary+last)//2)%M

def main():
    invariants=[]
    for m in range(1,31):
        B=s.Matrix(wedge(m));need(B==-B.T,'alternating form')
        if m>=2:
            U=peel(m);target=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix(wedge(m-2)))
            need(abs(U.det())==1 and U.T*B*U==target,'integral symplectic peel')
        if m%2:
            need(B.det()==1,'odd-length unimodular form')
            S=B[1:,1:];diag=smith_normal_form(S,domain=s.ZZ)
            need([abs(int(diag[i,i])) for i in range(m)]==[1]*(m-1)+[0],'odd interaction cokernel is free')
            need(S.row_join(B[0,1:].T).rank()>S.rank(),'odd-length affine center obstruction')
        else:
            v=s.Matrix(radical(m));need(B*v==s.zeros(m+1,1),'primitive radical')
            S=B[1:,1:];d=m//2-1
            diag=smith_normal_form(S,domain=s.ZZ);factors=[abs(int(diag[i,i])) for i in range(m)]
            expected=([0,0] if m==2 else [1]*(m-2)+[d,d])
            need(factors==expected and S.det()==d*d,'paired torsion factors')
            if d:
                center=v[1:,0]/d
                need(S*center==B[0,1:].T,'exact affine centering')
            invariants.append({'base_length':m,'d':d,'smith':factors,'radical':list(map(int,v))})
    # Classical dimension-four geometry is attached to this actual response matrix.
    B=wedge(3);vec=lambda a:tuple(a>>i&1 for i in range(4))
    planes={tuple(sorted((a,b,a^b))) for a,b in combinations(range(1,16),2) if pair(B,vec(a),vec(b))%2==0}
    need(len(planes)==15,'F2 Lagrangians')
    inside=[p for p in sorted(planes) if not any(a&1 for a in p)]
    meets=[p for p in sorted(planes) if any(a&1 for a in p)]
    need(len(inside)==3 and len(meets)==12,'affine Boolean chart guard')
    for plane in meets:
        x,y=[vec(a) for a in plane if a&1]
        need(jump(x[1:],y[1:])%4==0,'realized isotropic edge jump')
    checked=0
    for m in range(1,6):
        for w in product((0,1),repeat=m):
            for u in product((0,1),repeat=m):
                for e in range(1,5):
                    need(streamed(w,u,e)==(jump(w,u)//2)%(1<<e),'dyadic stream');checked+=1
    # A genuine four-source hourglass curl in the strong chart.
    w0=(1,0,0,0,0);w1=(1,1,0,0,0);u0=w0;u1=(1,0,0,1,0)
    curl=jump(w1,u1)-jump(w1,u0)-jump(w0,u1)+jump(w0,u0)
    need(curl==2,'nonzero integer rectangle curvature')
    measured=0
    for m in range(1,9):
        base=(1,)+(0,)*m+(0,);S=s.Matrix(wedge(m))[1:,1:]
        for i in range(m):
            for j in range(m):
                wi=list(base);uj=list(base);wi[i+1]=1;uj[j+1]=1
                wi=tuple(wi);uj=tuple(uj)
                # Four actual strong source pairs, counted independently by gap partitions.
                def D(w,u):return gap_h(w,u)-gap_h(w,u,False)
                response=D(wi,uj)-D(wi,base)-D(base,uj)+D(base,base)
                need(response==2*S[i,j],'strong Boolean measurement of interaction lattice');measured+=1
    output={'status':'PASS','integral_forms_through':30,'even_length_invariants':invariants,
            'lagrangians_F2_dim4':{'total':15,'within_zero_source_hyperplane':inside,'meeting_source_chart':meets},
            'stream_checks':checked,'strong_interaction_entries_measured':measured,
            'strong_rectangle':{'words':[w0,w1,u0,u1],'curl':curl}}
    (HERE/'geometry.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf8')
    print('PASS: unimodular peel, paired torsion, affine Lagrangian guard, dyadic stream')
if __name__=='__main__':main()
