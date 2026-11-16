from itertools import combinations
from math import gcd
from oracle import tournament,lift,path_dp
def need(x):
 if not x:raise AssertionError
def energy(x):
 bits=[i for i in range(x.bit_length()) if x>>i&1]
 return x*x+(x&1)+2*sum((1<<(j-i-1))+(j==i+1) for i,j in combinations(bits,2))
def generate(s,m):
 cursor=s+3; gap=s+3; dd=1
 if m>1:
  while dd<s+2 or pow(2,dd,m)!=1:dd+=1
 base=0; groups=[]
 for _ in range(m-1):
  x=cursor;cursor+=dd+gap+3
  y=cursor;cursor+=dd+gap+4
  base|=(1<<x)|(1<<(x+1))|(1<<y)|(1<<(y+2+dd))
  groups.append((1<<(x+1))|(1<<(x+1+dd))|(1<<(y+2))|(1<<(y+2+dd)))
 dy=[] if s==1 else [1]
 for j in range(1,s-1):
  dy.append((1<<cursor)|(1<<(cursor+j+1)))
  cursor+=j+1+gap
 def visit():
  for w in range(1<<len(dy)):
   z=base
   for i,bit in enumerate(dy):
    if w>>i&1:z|=bit
   yield z
   for mark in groups:
    z^=mark;yield z
 return cursor,visit()
def main():
 for s,m in [(1,1),(2,1),(7,1),(1,3),(2,5),(3,7),(4,9),(5,11)]:
  n,X=generate(s,m);mod=(1<<s)*m
  H={((1<<(2*n+2))+1-energy(x))%mod for x in X}
  need(H==set(range(1,mod,2)))
 for n in range(1,5):
  a=tournament(n,(1<<(n*(n-1)//2))-1)
  for x in range(1<<n):
   marks={i for i in range(n) if x>>(n-1-i)&1}
   W=lift([[[0]],a,a,[[0]]],[set(),marks,{n-1-i for i in marks},set()])
   need(path_dp(W)==(1<<(2*n+2))+1-energy(x))
 print("All CRT packet residue banks and 30 literal strong one-apex counts passed.")
if __name__=='__main__':main()
