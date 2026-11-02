#!/usr/bin/env python3
"""Exact rational-parameter circle certificate; custom bivariate polynomial engine."""
from fractions import Fraction as F
from audit_bridges import P,inv,check
import json

def power(p,n):
    out=P(1)
    for _ in range(n):out*=p
    return out

def run():
    t=P({(1,0):1});u=P({(0,1):1});d=u*u+u+1
    xn=[-u,u*(u+1),u+1]
    check(sum(xn)==d,'rational circle sum')
    check(sum(x*x for x in xn)==d*d,'rational circle square norm')
    Knum=[[P() for _ in range(5)] for _ in range(5)]
    for i in range(3):
        Knum[i][3]=t*xn[i];Knum[3][i]=-t*xn[i]
        Knum[i][4]=-t*xn[i];Knum[4][i]=t*xn[i]
    Knum[3][4]=t*d;Knum[4][3]=-t*d
    Anum=[[(d*int(i!=j)+Knum[i][j])*F(1,2) for j in range(5)] for i in range(5)]
    H,_,_,M,_,_,_=inv(Anum)
    # Actual A=Anum/d; H degree4 and M23 degree5 in adjacency entries.
    product_num=xn[0]*xn[1]*xn[2]
    expectedH=(P(F(15,2))+t*t*F(9,4))*power(d,4)+power(t,4)*product_num*d*F(1,4)
    expectedM=power(d,5)*F(5,8)+power(t,4)*product_num*power(d,2)*F(1,8)
    check(H==expectedH,'exact circle H identity after denominator clearing')
    check(M==expectedM,'exact circle M23 identity after denominator clearing')
    # At t=1, u=1 is the lower endpoint; u=-1 is upper endpoint.
    for v,r in [(F(1),F(-4,27)),(F(-1),F(0))]:
        den=d.at(F(1),v);prod=product_num.at(F(1),v)/den**3
        check(prod==r,'circle extremal product')
    # Cubic discriminant and endpoint factorizations are stated analytically.
    return {'status':'PASS','parameter':'s=(-u,u(u+1),u+1)/(u^2+u+1)','denominator_positive':'u^2+u+1=(u+1/2)^2+3/4','H':'15/2+9t^2/4+t*x1*x2*x3/4','M23':'5/8+t*x1*x2*x3/8','range':'[15/2+9t^2/4-t^4/27,15/2+9t^2/4]','scope':'3+2 regular rank-two slice only','missed_parameter_point':'(0,1,0), recovered by u→infinity (or direct zero product)'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
