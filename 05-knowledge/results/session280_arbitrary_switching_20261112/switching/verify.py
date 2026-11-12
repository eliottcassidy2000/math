"""Exact rational switched-interface certificates; standard library only."""
from fractions import Fraction as F
from itertools import product,combinations
from functools import lru_cache
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
R=((1,1,0),(1,1,1),(0,1,0));I=((1,0,0),(0,1,0),(0,0,1))
def check(b,msg):
    if not b:raise AssertionError(msg)
def mv(A,x):return tuple(sum(a*b for a,b in zip(row,x)) for row in A)
def mm(A,B):return tuple(tuple(sum(a*b for a,b in zip(row,col)) for col in zip(*B)) for row in A)
@lru_cache(None)
def g(p):
    if p==0:return 1
    a=p//2
    return g(a) if p%2 else g(a)+a+a%2
def block(w):
    p=int(w,2)+1;q=1+sum(1<<j for j,b in enumerate(w) if b=='0');t=1<<len(w)
    return ((g(p-1),q,-1),(p,t,0),(1,0,0))
def active(N):return ['1'+''.join(w)+'0' for n in range(2,N+1) for w in product('01',repeat=n-2)]
def psd(Q,strict=False):
    a,b=Q[0];_,c=Q[1];det=a*c-b*b
    return (a>0 and c>0 and det>0) if strict else (a>=0 and c>=0 and det>=0)
def sub(A,B):return tuple(tuple(a-b for a,b in zip(x,y)) for x,y in zip(A,B))
def scale(c,A):return tuple(tuple(c*x for x in row) for row in A)
def quad(Q,z):return sum(a*b for a,b in zip(z,mv(Q,z)))
def char(P):
    tr=sum(P[i][i] for i in range(3));c=sum(P[i][i]*P[j][j]-P[i][j]*P[j][i] for i,j in combinations(range(3),2))
    d=P[0][0]*(P[1][1]*P[2][2]-P[1][2]*P[2][1])-P[0][1]*(P[1][0]*P[2][2]-P[1][2]*P[2][0])+P[0][2]*(P[1][0]*P[2][1]-P[1][1]*P[2][0])
    return (1,-tr,c,-d)
def polyval(p,x):
    ans=0
    for c in p:ans=ans*x+c
    return ans
def main():
    data=json.loads((ROOT/'metric56.json').read_text());words=data['words'];tau=F(*data['tau']);D=data['Q_den']
    expected=[w for w in active(7) if sum(a!=b for a,b in zip(w,w[1:]))<=3]
    check(words==expected,'precise56 alphabet');check(len(words)==56,'size')
    Q=[((F(a,D),F(b,D)),(F(b,D),F(c,D))) for a,b,c in data['Q_num']];M=list(map(block,words));ell=F(1,2)
    for q in Q:check(psd(sub(q,((ell,0),(0,ell))),True),'uniform norm lower bound')
    minimum=None;tests=0
    for u,U in enumerate(M):
        for v,V in enumerate(M):
            q=V[0][1];t=V[1][1]
            for r in R:
                a,w,d=mv(V,mv(U,r));rr=F(w,a);s=F(d,a)
                K=((t-q*rr,rr),(-q*s,s))
                loss=sub(scale(tau,Q[u]),mm(tuple(zip(*K)),mm(Q[v],K)))
                check(psd(loss,True),('metric inequality',words[u],words[v],r))
                det=loss[0][0]*loss[1][1]-loss[0][1]**2
                minimum=det if minimum is None else min(minimum,det);tests+=1
    initial=F(0)
    for u,U in enumerate(M):
        rays=[mv(U,r) for r in R]
        positive=[(a,F(w+d,2),F(w-d,2)) for a,w,d in rays]
        for j in range(3):check(max(r[j] for r in positive)<=2*min(r[j] for r in positive),'ratio2')
        for x,y in combinations(rays,2):
            z=tuple(F(y[j])-F(y[0],x[0])*x[j] for j in (1,2))
            initial=max(initial,quad(Q[u],z))
            # Retain both ordered pairs: the scalar normalization is asymmetric.
            z=tuple(F(x[j])-F(x[0],y[0])*y[j] for j in (1,2))
            initial=max(initial,quad(Q[u],z))
    B=96;k=1;bound=4*B*B*initial/ell
    while bound>=1:k+=1;bound*=tau
    check(bound<1,'all-depth threshold')
    hostile=('10','101010','1000');P=I
    for w in hostile:
        f=char(block(w));check(polyval(f,1)<0 and polyval(f,-1)<0,'individual Pisot')
        P=mm(block(w),P)
    f=char(P);check(f==(1,-15315,-19568,-4096),'hostile cubic');check(polyval(f,-1)==156,'expanding conjugate')
    check(polyval(f,-2)<0,'root below minus1 bracket')
    tested=0
    for ws in product(active(7),repeat=3):
        if sum(map(len,ws))>=12:continue
        S=I
        for w in ws:S=mm(block(w),S)
        p=char(S);check(polyval(p,1)<0 and polyval(p,-1)<0,'no smaller total triple');tested+=1
    out={'status':'PROVED rational certificate, not Chan-Pak coverage','alphabet_size':len(words),'metric_PSD_tests':tests,'tau':str(tau),'ell':str(ell),'initial_energy':str(initial),'B':B,'uniform_threshold':k,'threshold_squared_width_bound_decimal_display':float(bound),'minimum_determinant':str(minimum),'hostile_words':hostile,'hostile_matrix':P,'hostile_characteristic':f,'hostile_f_minus1':polyval(f,-1),'shorter_total_triples_tested':tested,'certificate_sha256':hashlib.sha256((ROOT/'metric56.json').read_bytes()).hexdigest()}
    (ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ['alphabet_size','metric_PSD_tests','tau','initial_energy','uniform_threshold','hostile_f_minus1','shorter_total_triples_tested']},indent=2))
if __name__=='__main__':main()
