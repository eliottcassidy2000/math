"""Exact target-first decoder for a transitive one-head two-block seam."""
from math import isqrt
from itertools import combinations

def rev(w,x):return sum(((x>>i)&1)<<(w-1-i) for i in range(w))
def J(x):
 bits=[i for i in range(x.bit_length()) if x>>i&1]
 return sum((1<<(j-i-1))+(j==i+1) for i,j in combinations(bits,2))
def C(k,r,a,b):return 1+rev(k,b)-rev(r,a)-J(b)-J(a)
def pair(P,k,r):
 # Enumerate by shorter of left band, right band, and full divisor search.
 kl,kh,rl,rh=(1<<(k-1))+1,1<<k,1<<r,(3<<(r-1))-1
 sizes=(kh-kl+1,rh-rl+1,isqrt(P)); which=min(range(3),key=sizes.__getitem__)
 if which==0:
  for d in range(kl,kh+1):
   q,z=divmod(P,d)
   if z==0 and rl<=q<=rh:yield d,q
 elif which==1:
  for e in range(rl,rh+1):
   q,z=divmod(P,e)
   if z==0 and kl<=q<=kh:yield q,e
 else:
  for d in range(1,isqrt(P)+1):
   e,z=divmod(P,d)
   if z==0:
    if kl<=d<=kh and rl<=e<=rh:yield d,e
    if d!=e and kl<=e<=kh and rl<=d<=rh:yield e,d
def decode(N,k,r,t):
 if N<1 or N%2==0:return
 S=1<<t
 cmin=5-(3<<(r-1))-(1<<(k-1));cmax=(1<<k)-1
 low=max(1,(N-cmax+S-1)//S);high=(N-cmin)//S
 for P in range(low,high+1):
  for d,e in pair(P,k,r):
   a=e-(1<<r);b=(1<<k)-d
   if C(k,r,a,b)==N-P*S:yield a,b
def decode_all(N):
 # One spacer suffices, by relocating any excess spacer bits to a width.
 # H≥(3/8)2^n means 2^n≤floor(8N/3).
 lim=(8*N//3).bit_length()-1
 for k in range(1,lim-1):
  for r in range(1,lim-k):
   for a,b in decode(N,k,r,1):yield k,r,a,b

def need(x):
 if not x:raise AssertionError
def main():
 for a,b,n in [(153,334,1623869),(494,489,1623871),(241,382,1623873)]:
  need((a,b) in list(decode(n,10,10,1)))
  need(not list(decode(2**32+3,10,10,12)))
 for k in range(1,5):
  for r in range(1,5):
   for t in range(1,5):
    tab={}
    for a in range(1<<(r-1)):
     for b in range(1<<(k-1)):
      h=((1<<k)-b)*((1<<r)+a)*(1<<t)+C(k,r,a,b)
      tab.setdefault(h,set()).add((a,b))
    for h,rows in tab.items():
     need(set(decode(h,k,r,t))==rows)
     need(bool(list(decode_all(h))))
 for n in range(3,9):
  need(bool(list(decode_all(3*2**(n-2)+3))))
 need(not list(decode_all(7)) and not list(decode_all(21)))
 print("Physical bank targets and 64 exhaustive small seams invert exactly.")
if __name__=="__main__":main()
