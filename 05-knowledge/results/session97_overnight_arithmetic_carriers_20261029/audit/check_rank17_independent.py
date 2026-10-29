"""Stdlib coefficient-recursion section roots and rational height-matrix audit."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

def need(b,m):
    if not b:raise ArithmeticError(m)
def trim(a):
    while a and a[-1]==0:a.pop()
    return a

def add(*args):
    v=[F(0)]*max(map(len,args))
    for a in args:
        for i,x in enumerate(a):v[i]+=x
    return trim(v)
def scale(a,k):return trim([x*k for x in a])
def mul(a,b):
    v=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):v[i+j]+=x*y
    return trim(v)
def subtract(a,b):return add(a,scale(b,-1))
def rem(a,b):
    a=a[:]
    while len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]/b[-1]
        for j,v in enumerate(b):a[k+j]-=c*v
        trim(a)
    return a

def gcd_degree(a,b):
    while b:a,b=b,rem(a,b)
    return len(a)-1

def sqrtpoly(a,sgn):
    need((len(a)-1)%2==0,'even square degree');d=(len(a)-1)//2
    num,den=a[-1].numerator,a[-1].denominator;u,v=isqrt(num),isqrt(den)
    need(u*u==num and v*v==den,'leading square');out=[F(0)]*(d+1)
    out[d]=F(u if sgn=='+' else -u,v)
    for k in range(d-1,-1,-1):
        old=sum(out[i]*out[d+k-i] for i in range(k+1,d) if k+1<=d+k-i<=d)
        out[k]=(a[d+k]-old)/(2*out[d])
    need(mul(out,out)==a,'square polynomial identity');return out

def det(M):
    A=[list(map(F,r)) for r in M];ans=F(1)
    for j in range(len(A)):
        p=next((i for i in range(j,len(A)) if A[i][j]),None)
        if p is None:return F(0)
        if p!=j:A[j],A[p]=A[p],A[j];ans=-ans
        z=A[j][j];ans*=z
        for i in range(j+1,len(A)):
            ratio=A[i][j]/z
            for k in range(j,len(A)):A[i][k]-=ratio*A[j][k]
    return ans

def main():
    path=Path(__file__).resolve().parents[1]/'ranks/phase2/elkies_948_mathdata.json'
    data=json.loads(path.read_text());poly=lambda a:trim(list(reversed(list(map(F,a)))))
    A=scale(poly(data['S']),-27);B=scale(poly(data['T']),F(27,4));X=list(map(poly,data['x_sections']));Y=[]
    for x,s in zip(X,data['y_leading_signs']):Y.append(sqrtpoly(add(mul(mul(x,x),x),mul(A,x),B),s))
    need(Y[0]==poly(data['y1']),'named y1')
    gram=[[4 if i==j else 0 for j in range(17)] for i in range(17)];inf_total=0
    for i in range(17):
        for j in range(i):
            dx=subtract(X[i],X[j]);dy=subtract(Y[i],Y[j]);infinity=min(5-len(dx),7-len(dy));inf_total+=infinity
            gram[i][j]=gram[j][i]=2-gcd_degree(dx,dy)-infinity
    need(gram==data['published_gram'],'recomputed complete Gram')
    need(det(gram)==948,'rank17 determinant')
    minors=[det([r[:k] for r in gram[:k]]) for k in range(1,18)]
    need(all(v>0 for v in minors),'positive definite')
    for j,m in enumerate(data['secant_slopes'],1):need(subtract(Y[j],Y[0])==mul(poly(m),subtract(X[j],X[0])),'named secant')
    print(json.dumps({'independent_sections':17,'height_determinant':948,
        'leading_minors':[int(v) for v in minors],'infinity_intersections':inf_total,
        'secant_identities':15},sort_keys=True))
if __name__=='__main__':main()
