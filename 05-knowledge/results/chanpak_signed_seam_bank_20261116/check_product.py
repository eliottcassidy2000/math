from oracle import lift,path_dp,tournament,block
from itertools import combinations
def need(q):
 if not q:raise AssertionError
def TT(n):return [[int(i<j) for j in range(n)] for i in range(n)]
def rev(w,x):return sum(((x>>i)&1)<<(w-1-i) for i in range(w))
def J(V):
 return sum((1<<(v-u-1))+(v==u+1) for u,v in combinations(sorted(V),2))
def digits(w,z):return {i for i in range(w) if z>>i&1}
def count(k,r,t,a,b):
 A=digits(k,rev(k,b));B=digits(r,a)
 V=sorted(A|{v+k+t for v in B})
 N=k+r+t
 # Direct predecessor/successor expansion on the whole tape,
 # without grouping the spacer coefficient.
 E=sum(1<<v for v in V); L=sum(1<<(N-1-v) for v in V)
 NTT=(1<<N)+1+E-L-J(V)
 P=((1<<k)-b)*((1<<r)+a)
 C=1+rev(k,b)-rev(r,a)-J(digits(k,b))-J(digits(r,a))
 need(NTT==P*2**t+C)
 return P,C,NTT,(A,B)
def literal_tt_head(N, marked):
 # Exact status of sentinel and head; old vertices in either part have
 # a forced increasing order. Distinct from the pair-coefficient proof.
 F=(1<<N)-1
 H=2 # x first, and the empty second-slot case with t immediately before x
 for after in range(1,F+1):
  if (after&-after)&marked:continue
  H+=1 # sentinel just before x
  before=F^after
  if before and marked>>(before.bit_length()-1)&1:
   H+=1 # largest before vertex precedes x; sentinel is last
 return H

def main():
 for k in range(1,6):
  for r in range(1,6):
   for a in range(1<<(r-1)):
    for b in range(1<<(k-1)):
     for t in (1,2,7):count(k,r,t,a,b)
 for a,b,C in [(5,13,3),(3,11,5)]:
  P,c,H,(A,B)=count(5,4,1,a,b)
  need((P,c)==(399,C))
  W=lift([TT(5),TT(1),TT(4)],[A,set(),B])
  need(path_dp(W)==H)
 for a,b,C in [(153,334,-391),(494,489,-389),(241,382,-387)]:
  for t in (1,2,7,19):
   P,c,H,(A,B)=count(10,10,t,a,b)
   need((P,c)==(812130,C))
   if t==1:
    marked=sum(1<<v for v in A)+sum(1<<(v+11) for v in B)
    need(literal_tt_head(21,marked)==H)
 A=tournament(4,4);B=tournament(3,2)
 for mask in [{0},{0,2},{1,2,3}]:
  for V in [set(),{0},{0,1,2}]:
   h,q,E,LJ,L,_,_=block(A,mask)
   hh,qq,EE,JJ,LL,_,_=block(B,V)
   for t in (1,2,3):
    X,Y=2,4
    H=(q*X-L)*(qq*Y+EE)*2**t+h*hh+hh*(X*E-LJ)-h*(Y*LL+JJ)
    need(path_dp(lift([TT(1),A,TT(t),B,TT(2)],[set(),mask,set(),V,set()]))==H)
 print('All small products, two literal small tournaments, three growing mask states by literal assignment, and 27 nontransitive two-core literal counts verified.')
if __name__=='__main__':main()
