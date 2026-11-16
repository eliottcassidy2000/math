from math import isqrt
from itertools import product

def h(s):
 x=0; c=1;p=0
 for b in s:
  if b=='0':c+=x+p
  p=int(b);x=2*x+p
 return c

def rc(s):
 return s[::-1].translate(str.maketrans('01','10'))
def addr(x):return bin(x)[2:]
def isprime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def fixed(m):
 if m>=2:
  yield '0'+'1'*(m-2)+'0'
 if m==3:
  yield '001'
  yield '100'

def bank(B):
 # Sufficient construction; all active exteriors are by positive x,z
 for x in range(1,isqrt(B//16)+1,2):
  u,v=addr(x),rc(addr(x+1))
  a,b=u+'0110'+v,u+'1001'+v
  n=h(a); assert h(b)==n+2
  if n<=B:yield n,x,x+1,'',a,b
 for m in range(2,B.bit_length()-4):
  for w in fixed(m):
   for x in range(1,isqrt(B//(1<<(m+4)))+1):
    u,v=addr(x),rc(addr(x))
    a,b=u+'01'+w+'10'+v,u+'10'+w+'01'+v
    n=h(a); assert h(b)==n+2
    if n<=B:yield n,x,x,w,a,b

def select(P):
 # If canonical A->B has negative sign, reverse-complement interchanges
 # the canonical roles, preserves both counts, and normalizes to +2.
 for target in (P-2,P):
  if target<1: continue
  for m in [0]+list(range(2,target.bit_length()-4)):
   q=1<<(m+4)
   xmax=isqrt(target//q)
   for w in [''] if m==0 else fixed(m):
    def out(j):
     x = 2*j+1 if m==0 else j+1
     z = x+(m==0)
     u,v=addr(x),rc(addr(z))
     a=u+'01'+w+'10'+v
     b=u+'10'+w+'01'+v
     return h(a),a,b
    lo,hi=0,(xmax+1)//2 if m==0 else xmax
    while lo<hi:
     j=(lo+hi)//2
     N,a,b=out(j)
     if N<target:lo=j+1
     elif N>target:hi=j
     else:
      assert h(b)==N+2
      yield N,w,a,b
      break

def main():
 B=100000
 for p in (29989,30011,30013):
  nmax=2+(p-1).bit_length()-1
  assert not any(h('1'+format(i,'0'+str(n-2)+'b')+'0')==p
   for n in range(3,nmax+1) for i in range(1<<(n-2)))
  assert not list(select(p))
 assert {row[0] for row in select(10007)}=={10007}
 assert {row[0] for row in select(12541)}=={12539}

 b=sorted(bank(B))
 p=[x for x in b if isprime(x[0]) and isprime(x[0]+2)]
 print('closest to 30,000:', sorted(b,key=lambda v:abs(v[0]-30000))[:12])
 print('closest twin pairs:', sorted(p,key=lambda v:abs(v[0]-30000))[:12])
 print('windows at 30011 or 29989',[v for v in b if v[0] in (30009,30011,29987,29989)])
 for m in range(1,16):
  cand=set()
  for w in product('01',repeat=m):
   w=''.join(w); z=int(rc(w),2); x=int(w,2)
   if z-x+int(w[0])+int(w[-1])-1==2:cand.add(w)
  assert cand==set(fixed(m)),m

if __name__=='__main__':main()
