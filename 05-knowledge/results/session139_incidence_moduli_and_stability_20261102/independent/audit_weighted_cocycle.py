#!/usr/bin/env python3
"""Independent exact binary rounding of weighted four-face second moments."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from audit_bridges import P,matmul,transpose,det,check
from audit_integral_cocycle import eta,polynomial_eta
import argparse,json

def poly_family():
    a=P({(1,0):1});b=P({(0,1):1});n=5
    K=[[P() for _ in range(n)] for _ in range(n)]
    K[0][3]=a;K[3][0]=-a;K[1][4]=b;K[4][1]=-b
    w=[1,1,-2,0,0];Q=[[P(int(i==j)-F(w[i]*w[j],3)) for j in range(n)] for i in range(n)]
    return K,matmul(matmul(Q,K),transpose(Q))
def q_data(K):
    n=len(K);edges=list(combinations(range(n),2))
    q={e:(1-K[e[0]][e[1]]*K[e[0]][e[1]])*F(1,4) for e in edges}
    V=sum(q.values());M22=sum(q[e]*q[f] for e,f in combinations(edges,2) if not set(e)&set(f))
    e4=sum(det([[K[i][j] for j in s] for i in s]) for s in combinations(range(n),4))
    return V,M22,e4

def signed_mean_formula(K,s):
    i,j,k,l=s;a,b,c,d,e,f=K[i][j],K[i][k],K[i][l],K[j][k],K[j][l],K[k][l]
    return (a*b*d-a*c*e+b*c*f-d*e*f)/4

def enumerate_rounding(K,mutate=False):
    n=len(K);edges=list(combinations(range(n),2));faces=list(combinations(range(n),4))
    probs={(i,j):(1+K[i][j])/2 for i,j in edges}
    means={s:F(0) for s in faces};seconds=means.copy();pfseconds=means.copy();mass=F(0)
    for mask in range(1<<len(edges)):
        prob=F(1);S=[[0]*n for _ in range(n)]
        for bit,(i,j) in enumerate(edges):
            z=(mask>>bit)&1;prob*=probs[i,j] if z else 1-probs[i,j]
            S[i][j]=2*z-1;S[j][i]=-S[i][j]
        mass+=prob
        for s in faces:
            e=eta(S,s);i,j,k,l=s;pf=S[i][j]*S[k][l]-S[i][k]*S[j][l]+S[i][l]*S[j][k]
            means[s]+=prob*e;seconds[s]+=prob*e*e;pfseconds[s]+=prob*pf*pf
    check(mass==1,'rounding mass1')
    for s in faces:
        i,j,k,l=s
        pf=K[i][j]*K[k][l]-K[i][k]*K[j][l]+K[i][l]*K[j][k]
        squares=(K[i][j]*K[k][l])**2+(K[i][k]*K[j][l])**2+(K[i][l]*K[j][k])**2
        check(pfseconds[s]==pf*pf+3-squares,'Pf second moment')
        check(seconds[s]==(pf*pf+(1 if mutate else 2)-squares)/8,'eta second moment')
        check(means[s]==signed_mean_formula(K,s),'eta mean cubic')
    V,M22,e4=q_data(K)
    # q_data determinant returns P because independent determinant is polynomial-generic.
    e4=e4.at(0,0) if isinstance(e4,P) else e4
    predicted=(e4-comb(n,4)+4*comb(n-2,2)*V-16*M22)/8
    check(sum(seconds.values())==predicted,'global matching-incidence formula')
    return {'total_bad_probability':str(sum(seconds.values())),'mean_charges':[str(means[s]) for s in faces], 'V':str(V),'M22':str(M22),'e4':str(e4)}

def run(mutate):
    K,L=poly_family();V,M,e4=q_data(K);V2,M2,e42=q_data(L)
    check(V==V2 and e4==e42,'family fixes V/e4')
    check(M2-M==P({(2,2):F(-1,54)}),'symbolic M22 difference')
    point=F(1,4),F(1,2)
    A=[[F(p.at(*point)) for p in row] for row in K];B=[[F(p.at(*point)) for p in row] for row in L]
    ra=enumerate_rounding(A,mutate);rb=enumerate_rounding(B,mutate)
    check(F(rb['M22'])-F(ra['M22'])==F(-1,3456),'M22 rational hostile difference')
    check(F(rb['total_bad_probability'])-F(ra['total_bad_probability'])==F(1,1728),'expected-support rational hostile difference')
    check(set(ra['mean_charges']+rb['mean_charges'])=={'0'},'same zero mean cocycle')
    uniform=enumerate_rounding([[F(0) for _ in range(5)] for _ in range(5)],mutate)
    check(uniform['total_bad_probability']=='5/4','uniform expected bad count')
    # Unstructured rational input ensures identities are not only tested on the bipartite family.
    K4=[[F(0) for _ in range(4)] for _ in range(4)]
    for (i,j),v in zip(combinations(range(4),2),[F(1,2),F(-1,3),F(3,4),F(-2,5),F(1,7),F(4,9)]):K4[i][j]=v;K4[j][i]=-v
    general=enumerate_rounding(K4,mutate)
    return {'status':'PASS','binary_roundings':3*1024+64,'symbolic_M22_difference':'-a^2*b^2/54','symbolic_expected_bad_difference':'a^2*b^2/27','A':ra,'B':rb,'uniform':uniform,'unstructured_four':general}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mutate-second-constant',action='store_true');a=ap.parse_args()
    print(json.dumps(run(a.mutate_second_constant),sort_keys=True,indent=2))
