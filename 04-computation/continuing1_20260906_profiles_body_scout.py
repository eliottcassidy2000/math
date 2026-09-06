"""Exact bounded search for a physical gain from the rank-four tail gate.

Universe: every eleven-subset of 1..24, with body dilations 1 and 2.
This is a bounded negative search, not a universal redundancy theorem.
"""
from itertools import combinations
from math import lcm, comb
from fractions import Fraction as F
import sys
sys.stdout.reconfigure(newline="\n")

def need(ok,why):
 if not ok:raise ArithmeticError(why)

def main():
 N=24;D=14*lcm(*range(1,N+1))
 safe={h:[((14*k+1)*(D//(14*h)),(14*k+13)*(D//(14*h))) for k in range(h)] for h in range(1,N+1)}
 def lengths(H):
  A=safe[H[0]]
  for h in H[1:]:
   B=safe[h];out=[];i=j=0
   while i<len(A) and j<len(B):
    x=max(A[i][0],B[j][0]);y=min(A[i][1],B[j][1])
    if x<=y:out.append((x,y))
    if A[i][1]<B[j][1]:i+=1
    else:j+=1
   A=out
  return sorted((y-x for x,y in A if y>x),reverse=True)
 count=0;near=[];hits1=[];hits2=[]
 for H in combinations(range(1,N+1),11):
  count+=1;W=lengths(H);M=sum(W);L=max(W,default=0)
  # Necessity for a gain solely at the (19,25) intermediate rank:
  # (7,13) must be paid by mass, while the (19,25) mass gate fails.
  if M*49<2*D or M*3325>=138*D:continue
  if L*3325<37*D and L*287>=2*D:
   near.append((H,F(M,D),F(L,D),F(sum(W[:4]),D)))
   if sum(W[:4])*3325>=120*D:hits1.append(H)
  # Doubling every body speed repeats each component at half its length.
  # A mixed-parity body is paired already; an all-odd body's singleton
  # central component is handled by the literal replicated list here.
  W2=sorted([w for w in W for _ in range(2)],reverse=True)
  if L*3325<74*D and L*287>=4*D and sum(W2[:4])*3325>=240*D:hits2.append(H)
 need(count==comb(24,11)==2496144,'entire body universe')
 need(not hits1 and not hits2,'no rank-four improvement in the declared body/dilation universe')
 print('FINITE-EXACT BOUNDED NEGATIVE; no universal redundancy assertion')
 print('UNIVERSE: all 2496144 eleven-subsets of1..24; body dilations1,2; strict rank-four gain over five mass/width gates')
 print('INTEGER_GRID_DENOMINATOR:',D)
 print('HITS_DILATION1:',len(hits1),'HITS_DILATION2:',len(hits2))
 print('NEAR_DILATION1:',len(near))
 for H,M,L,S in near:print('BODY',H,'mass',M,'max',L,'top4',S)
 print('STOPPING_REASON: no physical gain found; retain the proved sufficient criterion and seek body geometry rather than claiming family closure')
 print('PASS')

if __name__=='__main__':main()
