#!/usr/bin/env python3
"""Independent rational polynomial audit of regular cospectral critical-line hostile."""
from fractions import Fraction as F
from math import comb
import argparse,json
from audit_bridges import P,matmul,transpose,inv,check

def zprod(f,g):
    h=[P() for _ in range(len(f)+len(g)-1)]
    for i,x in enumerate(f):
        for j,y in enumerate(g):h[i+j]+=x*y
    return h
def regular_family(mutate):
    t=P({(1,0):1});n=5
    K=[[P() for _ in range(n)] for _ in range(n)]
    for i,j in [(0,3),(3,4),(4,0)]:K[i][j]=t;K[j][i]=-t
    w=[1,1,-2,0,0];Q=[[P(int(i==j)-F(w[i]*w[j],3)) for j in range(n)] for i in range(n)]
    L=matmul(matmul(Q,K),transpose(Q))
    matrices=[]
    for J in [K,L]:
        check(all(sum(row)==0 for row in J),'regular skew row sums')
        A=[[(int(i!=j)+J[i][j])*F(1,2) for j in range(n)] for i in range(n)]
        for row in A:
            check(sum(row)==2,'regular weighted row sum')
            for p in row:
                c=p.d.get((0,0),0);r=sum(abs(v) for k,v in p.d.items() if k!=(0,0))
                check(0<=c-r<=c+r<=1,'whole parameter interval weighted')
        matrices.append(A)
    HA,VA,ta,MA,ua,va,pa=inv(matrices[0]);HB,VB,tb,MB,ub,vb,pb=inv(matrices[1])
    check(pa==pb and VA==VB and ta==tb and ua==ub and va==vb,'regular spectral invariants')
    check(HA==P(F(15,2))+t*t*F(9,4),'source H regular formula')
    check(HB-HA==t*t*t*t*F(-1,28 if mutate else 27),'regular H difference')
    check(MB-MA==t*t*t*t*F(-1,54),'regular M23 difference')
    expected=zprod(zprod([P(1),P(-2)],zprod([P(1),P(F(1,2))],[P(1),P(F(1,2))])),[P(1),P(1),P(F(1,4))+t*t*F(3,4)])
    check(pa==expected,'critical line spectral factorization')
    # Divide P by the Perron factor1-2z, then transform Q by z->-z/(1+z).
    q=[pa[0]]
    for j in range(1,5):q.append(pa[j]+2*q[j-1])
    check(pa[5]==-2*q[4],'Perron quotient exact')
    transformed=[P() for _ in range(5)]
    for j,c in enumerate(q):
        for k in range(5-j):transformed[j+k]+=c*((-1)**j)*comb(4-j,k)
    check(q==transformed,'reduced functional equation')
    vertex_delta=matmul(matrices[1],matrices[1])[0][0]-matmul(matrices[0],matrices[0])[0][0]
    check(vertex_delta==t*t*F(5,18),'individual vertex-root response is not preserved')
    point=F(1,4),F(0)
    check((HB-HA).at(*point)==F(-1,6912),'regular rational hostile')
    return {'status':'PASS','symbolic_H_A':str(HA),'symbolic_H_B_minus_H_A':str(HB-HA),'rational_t':'1/4','H_A':str(HA.at(*point)),'H_B':str(HB.at(*point)),'difference':'-1/6912','V':str(VA.at(*point)),'tau':str(ta.at(*point)),'u':str(ua.at(*point)),'v':str(va.at(*point)),'reduced_functional_equation':'Q(z)=(1+z)^4 Q(-z/(1+z))','weighted_interval':'-1<=t<=1','vertex_root_second_moment_difference':'5t^2/18'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mutate-h-coefficient',action='store_true');a=ap.parse_args()
    print(json.dumps(regular_family(a.mutate_h_coefficient),indent=2,sort_keys=True))
