#!/usr/bin/env python3
from fractions import Fraction as F
from itertools import combinations
import json

def need(x,label):
    if not x: raise RuntimeError(label)
def eye(d): return [[F(i==j) for j in range(d)] for i in range(d)]
def det(a):
    a=[r[:] for r in a]; out=F(1)
    for j in range(len(a)):
        p=next((i for i in range(j,len(a)) if a[i][j]),None)
        if p is None:return F(0)
        if p!=j:a[p],a[j]=a[j],a[p];out=-out
        z=a[j][j];out*=z
        for i in range(j+1,len(a)):
            t=a[i][j]/z
            for h in range(j+1,len(a)):a[i][h]-=t*a[j][h]
    return out

def inv(a):
    d=len(a);b=[row[:]+e for row,e in zip(a,eye(d))]
    for j in range(d):
        p=next(i for i in range(j,d) if b[i][j]);b[p],b[j]=b[j],b[p]
        t=b[j][j];b[j]=[x/t for x in b[j]]
        for i in range(d):
            if i!=j:
                t=b[i][j];b[i]=[x-t*y for x,y in zip(b[i],b[j])]
    return [r[d:] for r in b]
def bvec(k,i,j):return [F(h==i)-F(h==j) for h in range(k-1)]
def quad(a,b):return sum(b[i]*a[i][j]*b[j] for i in range(len(b)) for j in range(len(b)))
def network(v,corr=F(0),pair=None):
    k=len(v);d=k-1;r=[F(1,x) for x in v];cs=[r[i]+r[-1] for i in range(d)]
    s=eye(d)
    if pair:
        i,j=pair;s[i][j]=s[j][i]=corr
    sigma=[[cs[i]*s[i][j]*cs[j] for j in range(d)] for i in range(d)]
    L=inv(sigma); edges={}
    for i,j in combinations(range(k),2):
        edges[i,j]=sum(L[i]) if j==k-1 else -L[i][j]
    need(all(x>=0 for x in edges.values()),'nonnegative conductance')
    rebuilt=[[sum(w*bvec(k,i,j)[a]*bvec(k,i,j)[b] for (i,j),w in edges.items()) for b in range(d)] for a in range(d)]
    need(rebuilt==L,'Laplacian reconstruction')
    resist={}
    for (i,j),w in edges.items():
        z=quad(sigma,bvec(k,i,j));cap=(r[i]+r[j])**2
        need(z<=cap,'all pair resistance bounds')
        need(not w or z==cap,'active contact equality')
        resist[i,j]=z
    need(sum(w*(r[i]+r[j])**2 for (i,j),w in edges.items())==d,'contact trace')
    tree_sum=F(0);support=[e for e,w in edges.items() if w]
    for T in combinations(support,d):
        B=[[bvec(k,*e)[h] for e in T] for h in range(d)]
        minor=det(B);need(minor in (-1,0,1),'incidence minor')
        p=F(1)
        for e in T:p*=edges[e]
        tree_sum+=minor**2*p
    need(tree_sum==det(L),'matrix tree identity')
    for den in range(1,15):
        for num in range(1,den+1):
            t=F(num,den);T=[[t*z for z in row] for row in sigma]
            need(det(T)<=det(sigma),'feasible scaled determinant')
            trace=sum(L[i][j]*T[j][i] for i in range(d) for j in range(d))
            need(det(T)*det(L)<=(trace/d)**d,'AMGM certificate')
    if pair:
        e=tuple(sorted(pair));w=edges[e];need(w>0,'extra contact edge')
        b=bvec(k,*e);cut=[[L[i][j]-w*b[i]*b[j] for j in range(d)] for i in range(d)]
        need(det(cut)>0,'deleted edge graph stays connected')
        cutcov=inv(cut)
        need(any(quad(cutcov,bvec(k,i,j))>(r[i]+r[j])**2 for i,j in combinations(range(k),2)),'tree deletion is not feasible')
    return {'speeds':v,'det_covariance_in_grounded_potential_coordinates':str(det(sigma)),
            'tree_sum':str(tree_sum),'conductances':{str(e):str(w) for e,w in edges.items()},
            'resistances':{str(e):str(z) for e,z in resist.items()}}
print(json.dumps({'status':'PASS','examples':[network([1,2,4]),network([1,2,3,5],F(3,28),(1,2))],
                  'scope':'Exact electrical contact certificates; no LRC orbit-hit claim.'},indent=2))
