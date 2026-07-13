import math
from math import gcd
def lcm(xs):
    r=1
    for x in xs:
        if x: r=r*x//gcd(r,x)
    return r
def moments(E,Ng):
    s1=s2=0
    for k in range(1,Ng):
        x=k/Ng; occ=0
        for e in E: occ|=1<<(int((e*x%1.0)*7.0)%7)
        N=7-bin(occ).count("1"); s1+=N; s2+=N*(N-1)
    n=Ng-1; return s1/n,s2/n
def J(E,Ng=90000): m1,m2=moments(E,Ng); return 6*m1-m2
def J_inf(C,Ng=90000): m1,m2=moments(C,Ng); return 6*(6/7)*m1-(5/7)*m2
# k=9: compact 8-cluster C, far element added -> 9-core
print("k=9 J two-scale constant (controlled Ng=400*w):")
for C in [[0,1,2,3,4,5,6,7],[0,1,2,3,4,5,6,8],[0,1,2,3,4,5,6,9]]:
    L=lcm(C); worst=0
    for tag,w in [("prime",101),("prime",211),("prime",503),("lcmM",L),("lcmM",2*L),("lcmM",5*L)]:
        Ng=max(90000,400*w)
        e=abs(J(C+[w],Ng)-J_inf(C,Ng)); worst=max(worst,e*w)
    print(f"  C={C} lcm={L} Jinf={J_inf(C):.4f}  worst err*w={worst:.3f}")
# J threshold margin: J_inf(consec8)=? vs 432/91=4.7473
print(f"  J_inf(consec8)={J_inf([0,1,2,3,4,5,6,7]):.4f}, threshold=432/91={432/91:.4f}, margin={J_inf([0,1,2,3,4,5,6,7])-432/91:.4f}")
