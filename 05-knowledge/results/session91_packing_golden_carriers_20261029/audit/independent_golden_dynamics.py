#!/usr/bin/env python3
"""Independent one-operation separation and nonlumpability of golden fibres."""
from independent_golden_quotients import words,read,key,modulus
import json
N=0
def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)
def toggle(w,j):
 a=list(w);a[j]^=1
 return tuple(a) if all(not(a[i] and a[(i+1)%len(a)]) for i in range(len(a))) else w

def main():
 rows=[]
 for n in range(1,19):
  g=modulus(n);W=list(words(n));fib={}
  for w in W:fib.setdefault(key(read(w),g),[]).append(w)
  if n%2:
   ck(len(fib)==len(W),'odd no loss');continue
  Z=fib[(0,0)];ck(len(Z)==3,'zero triple')
  probabilities=[]
  for w in Z:
   probabilities.append(sum(key(read(toggle(w,j)),g)==(0,0) for j in range(n)))
  for j in range(n):
   profiles={(key(read(w),g),key(read(toggle(w,j)),g)) for w in W}
   ck(len(profiles)==(len(W) if n>=4 else 1),'single operation exact residual separation')
  ck(sorted(probabilities)==([0,n//2,n//2] if n>=4 else[2,2,2]),'zero transition probability')
  rows.append(dict(n=n,words=len(W),quotient=len(fib),one_step_states=len(profiles),zero_transition_numerators=sorted(probabilities),denominator=n))
 print(json.dumps(dict(status='PASS',checks=N,rows=rows),indent=2))
if __name__=='__main__':main()
