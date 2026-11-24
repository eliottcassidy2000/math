from pathlib import Path
from fractions import Fraction as F
from math import factorial, comb
import sys, subprocess

p=Path(__file__).resolve().parent
sys.path[:0]=[str(p.parent/'moore_chanpak_counts_20261124'),
             str(p.parent/'moore57_compact_blue_routes_20261124')]
import counts as a
import compact as b

def upper_factor(d):
    t=factorial(d)<<(40*d); l,h=0,1<<47
    while h-l>1:
        mid=(h+l)//2
        if mid**d>=t:h=mid
        else:l=mid
    assert (h-1)**d<t<=h**d
    return h

def unhex(s):
    a,b=s.lower().split('p'); hg,lo=a[2:].split('.')
    return F(int(hg+lo,16))*F(2)**(int(b)-4*len(lo))

T,_=b.field_host(5)
P=a.hs_spine()
Ts=[T,a.tournament(50,lambda i,j:j==i+1 or not a.hs_edge(P[i],P[j])),
      a.tournament(50,lambda i,j:j==i+1 or a.hs_edge(P[i],P[j]))]
n=50; mu=F(factorial(n),2**(n-1))
v=sum((F(2**r*comb(j-1,r-1)*comb(n-2*j,r)*factorial(n-2*j),factorial(n))
      for j in range(1,(n-1)//2+1) for r in range(1,min(j,n-2*j)+1)),F())
print("relative variance",float(v))
for k,t in enumerate(Ts):
    name=['compact','native sparse','native dense'][k]
    low=a.fixed_block(t,[list(range(i,i+5)) for i in range(0,50,5)])
    txt=' '.join(map(str,[upper_factor(i) for i in range(1,52)]+
      [x|(1<<50) for x in t]+[(1<<50)-1]))+'\n'
    log,ceiling=subprocess.check_output([str(p/'hs-bracket'),str([9,4,4][k])],
      input=txt,text=True).split()
    r=unhex(ceiling)/mu
    assert r<[F(80,100),F(51,100),F(53,100)][k]
    greater=(1-r)**2/(v+(1-r)**2)
    assert greater>[F(50,100),F(85,100),F(84,100)][k]
    print(name,low,"binary log upper",log,"absolute hex upper",ceiling,
          "labelled above fraction >",float(greater))
