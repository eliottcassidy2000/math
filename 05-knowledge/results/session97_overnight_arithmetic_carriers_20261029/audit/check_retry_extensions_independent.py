"""Independent finite guards and clocks audit for cheap-controller retries."""
from fractions import Fraction as F
from itertools import product
import json

def need(b,m):
    if not b:raise ArithmeticError(m)
# Native lengths J, division clocks D, carry C, source residue and modulus.
C=[(1,1,1,3,4),(4,6,179,125,128),(7,11,16177,1357,4096),
   (7,11,12289,2269,4096),(7,11,9697,2877,4096)]
def edge(n):
    t=3*n+1;a=0
    while t%2==0:t//=2;a+=1
    return t,a

def main():
    p=[F(1,1<<d) for j,d,c,r,m in C]
    rho=sum(p)
    need(rho==F(1059,2048),'failure measure')
    need(sum(w*c[0] for w,c in zip(p,C))/(1-rho)==F(1173,989),'odd clock mean')
    need(sum(w*c[1] for w,c in zip(p,C))/(1-rho)==F(1249,989),'division clock mean')
    alpha=[-F(c,3**j-2**d) for j,d,c,r,m in C]
    need(min(alpha)==-F(16177,139) and max(alpha)==-1,'fixedpoint interval')
    edges={}
    for n in range(-115,0,2):
        hits=[(j,d,c) for j,d,c,r,m in C if n%m==r]
        if not hits:continue
        need(len(hits)==1,'disjoint guards')
        j,d,c=hits[0];y=F(3**j*n+c,2**d)
        need(y.denominator==1,'guarded integer')
        edges[n]=int(y)
        need(y<n or (n,y) in [(-3,-1),(-1,-1)],'cycle exclusion')
    need(edges[-3]==-1 and edges[-1]==-1,'two exceptional ports')
    rays=0
    for j,d,c,r,m in C[1:]:
        gap=2**(d+2)-3**(j+1);carry=3*c+2**d
        need(gap>0 and gap*r>carry,'uniform original-source repair')
        for k in range(100):
            n=r+2*m*k;y=(3**j*n+c)//2**d;z,a=edge(y)
            need(a>=2 and z<n,'repaired even k half-ray');rays+=1
    for t in range(2,301):
        n=2**t-1;z=n;D=0
        for i in range(t):z,a=edge(z);D+=a
        c=1 if t%2 else 2+(t&-t).bit_length()-1
        need(D==t+c and z==F(3**t-1,2**c),'compressed clock')
        need((z<n)==(t in [2,4,8]),'Mersenne first-return classifier')
    # Independently compose every pair through its original source pullback.
    for left,right in product(C,repeat=2):
        j,d,c,r,m=left;k,e,b,s,l=right
        A=3**(j+k);D=d+e;B=3**k*c+2**d*b;M=2**(D+1)
        R=((2**D-B)*pow(A,-1,M))%M
        for shift in [0,1,7]:
            n=R+shift*M;y=F(3**j*n+c,2**d);z=F(3**k*y+b,2**e)
            need(n%m==r and y.denominator==1 and int(y)%l==s,'pulled-back guards')
            need(z==F(A*n+B,2**D) and n<y<z,'named intermediate port')
    print(json.dumps({'pair_charts':25,'pair_sources':75,'repaired_sources':rays,
        'rho':str(rho),'periodic_interval':[-115,-1],'only_cycle':[-1],
        'Mersenne_t_checked':[2,300],'clock_means':['1173/989','1249/989']},sort_keys=True))
if __name__=='__main__':main()
