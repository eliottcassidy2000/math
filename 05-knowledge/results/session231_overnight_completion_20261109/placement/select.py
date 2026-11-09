"""Exact finite inverse for the twelve-bit interval + arbitrary singleton contexts.
No claim that every target is accepted. python3 select.py 9909 9907 404037039
"""
import json,sys

def f(bits):
 return sum((1<<(j-i-1))+(j==i+1) for i in range(len(bits)) if bits[i] for j in range(i+1,len(bits)) if bits[j])
def digits(p,n):return [(p>>(n-i-1))&1 for i in range(n)]
def data(k):
 r=4096**k; j=(r-1)//4095
 q=11520*r; A=7516*r+343260*j+45;G=2353+9532800*j+5760*r
 D=(245051264*r+1602510*k+1118023)//24843
 return q,A,G,2*A+45-D,D

def baseline(k,P,R,b):
 q,A,G,T,D=data(k);pb=digits(P,P.bit_length());rb=digits(R,b)
 rev=sum(x<<i for i,x in enumerate(rb));e=rb[0] if rb else 0
 return ((q*P+A)*((1<<(b+1))-rev)+45*(2*R+1-f(rb))-P*G-45*f(pb)-T-45*(P%2)-45*e)

def select(N):
 if N<1 or N%2==0:return None
 k=0
 while data(k)[4]-2*k<=N:
  cap=(N-45)//(5760*4096**k)
  if cap<1:break
  for m in range(cap.bit_length()):
   for a in range(m+1):
    b=m-a
    for P in (range(1<<(a-1),1<<a) if a else [0]):
     for R in range(1<<b):
      delta=baseline(k,P,R,b)-N
      if 0<=delta<=2*k and delta%2==0:
       return dict(k=k,j=delta//2,prefix=P,suffix=R,suffix_length=b,vertices=12*k+13+m)
  k+=1
 return None
if __name__=='__main__':
 for s in sys.argv[1:]:
  n=int(s);print(json.dumps(dict(target=n,witness=select(n))))
