from itertools import product
from math import isqrt
from slot_profiles import profile,append,ordinal,external,strong,end_profiles
A=((0,41,11,1,15,29),12)
B=((28,17,18,22,0,31),43)
Q=[profile(*R) for R in (A,B)]
if Q!=[(5,128,34,41,0,1,-3),(5,128,106,41,0,5,-27)]:
 raise ValueError(Q)

def state(bits):
 s=append((1,1,1,1),Q[0])
 for e in bits:s=append(s,Q[e])
 return s

def graph(bits):
 a=(0,);q=1
 for u in (0,*bits):
  b,tag=(A,B)[u]
  q|=tag<<len(a);a=ordinal(a,b)
 return external(a,q,2)

def count(bits):
 j,w,k,e=state(bits)
 n=len(bits)
 if (j,e,3*k+w)!=(5**(n+1),0,5**n*(27+10*n)):raise ValueError(bits)
 return 4*w+k

def length(n):
 if n<594:return -1
 m=0;lower=594
 while 128*lower<=n:m+=1;lower*=128
 return m

def decode(n):
 m=length(n)
 if m<0 or n>=625*128**m:return None
 c=5**m*(27+10*m)
 w,r=divmod(3*n-c,11)
 if r or w<0:return None
 ans=[]
 for i in range(m,0,-1):
  rec=[b for b in (0,1) if (w-5**i*(34+72*b))%128==0]
  if not rec:return None
  ans+=rec;w=(w-5**i*(34+72*rec[0]))//128
 bits=tuple(reversed(ans))
 return bits if w==162 and count(bits)==n else None

def pack(m):
 w0=162*128**m+170*(128**m-5**m)//123
 delta=[72*5**(i+1)*128**(m-1-i) for i in range(m)]
 return w0,delta,5**m*(27+10*m)

def greedy(w,weights):
 bits=[]
 for x in weights:
  b=w>=x;bits.append(int(b))
  if b:w-=x
 return tuple(bits)

def predecessor(n):
 m=length(n)
 while m>=0:
  w0,v,c=pack(m)
  budget=(3*n-c)//11-w0
  if budget>=0:return count(greedy(budget,v))
  m-=1
 return None

def successor(n):
 m=max(0,length(n))
 while True:
  w0,v,c=pack(m)
  all_budget=3*(count((1,)*m)-n)//11
  if all_budget>=0:return count(tuple(1-i for i in greedy(all_budget,v)))
  m+=1

for m in range(7):
 found=[]
 for bits in product((0,1),repeat=m):
  t=count(bits);found.append(t)
  if decode(t)!=bits or not 594*128**m<=t<625*128**m:raise ValueError((m,bits,t))
  if m in (0,1):
   a=graph(bits)
   if not strong(a) or end_profiles(a)[0][-1]!=t:raise ValueError((m,bits,t))
 if found!=sorted(set(found)):raise ValueError('not strict lex')
 for lo,hi in zip(found,found[1:]):
  n=(lo+hi)//2
  if predecessor(n)!=lo or successor(n)!=hi:raise ValueError((lo,hi,n))
 print(m,len(found),found[0],found[-1])
for word,p in (((0,),76717),((0,1,1,0),160923394019)):
 if count(word)!=p or not all(p%d for d in range(2,isqrt(p)+1)) or not strong(graph(word)):
  raise ValueError(p)
 print('prime',word,p)
