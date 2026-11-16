import json, sys
from pathlib import Path
from itertools import permutations
from collections import Counter
from fractions import Fraction

def need(check):
 if not check: raise ArithmeticError("certificate failure")
def det(a):
 n=len(a); out=0
 for p in permutations(range(n)):
  v=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
  for i in range(n):v*=a[i][p[i]]
  out+=v
 return out
def mm(A,B,m):
 return [[sum(a*b for a,b in zip(x,y))%m for y in zip(*B)] for x in A]
def mv(A,v,m):
 return tuple(sum(a*b for a,b in zip(row,v))%m for row in A)
I=[[int(i==j) for j in range(4)] for i in range(4)]
def power(M,n,m):
 A=I
 while n:
  if n&1: A=mm(A,M,m)
  M=mm(M,M,m);n//=2
 return A
def main(file):
 d=json.loads(Path(file).read_text()); j=d['words'].index([0,6,5])
 q=d['matrices'][j]; M=[q[i:i+4] for i in range(0,16,4)]
 need(q==[5632,2235,1222,4221,672,267,149,506,-680,-272,-147,-511,-6896,-2736,-1498,-5169])
 s=tuple(d['start']);f=d['readout'];need(s==(-43,-5,4,53) and f==[32,14,10,27])
 def H(v,m):return sum(x*y for x,y in zip(v,f))%m
 for t in range(5):
  pol=t**4-583*t**3+1025*t**2-1047*t-424
  need(pol==det([[t*I[i][j]-M[i][j] for j in range(4)] for i in range(4)]))
  need((pol-(t-3)**2*(t*t-2*t-1))%5==0)
 N=[[M[i][j]-3*I[i][j] for j in range(4)] for i in range(4)]
 u=mv(N,(0,0,0,1),5);need(u==(1,1,4,3) and mv(N,u,5)==(0,0,0,0))
 need(power(M,60,5)==I)
 need(all(mv(power(M,a,5),s,5)!=tuple(t%5 for t in s) for a in (12,20,30)))
 O=[mm([f],power(M,a,5),5)[0] for a in range(4)]
 need(det(O)%5==1)
 need(tuple((x-y)%25 for x,y in zip(mv(power(M,60,25),s,25),s))==(5,5,20,15))
 phases=[Counter() for _ in range(4)]
 for j in range(60):
  phases[j%4][H(mv(power(M,j,5),s,5),5)]+=1
  need(H(mv(power(M,j,5),u,5),5)==(2,1,3,4)[j%4])
 need(all(q==Counter({i:3 for i in range(5)}) for q in phases))
 need(mv(power(M,17,32),s,32)==mv(M,s,32))
 need(all(H(mv(power(M,j,32),s,32),32)==(25,29,17,1)[j%4] for j in range(17)))
 A=168;D=8**8;B=13460883;z=Fraction(A,100*D)
 need(B>A**3 and (A+D*z)**300*z**(-3)<B**97)
 print("Golden block characteristic, state/observability, all 60 five-adic derivative phases, mod-32 readout and exact sparse-growth inequalities PASS")
if __name__=='__main__':main(sys.argv[1])
