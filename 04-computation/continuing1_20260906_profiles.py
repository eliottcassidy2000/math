"""Ranked all-height component envelope: exact finite head and hostiles."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0

def need(ok,why):
 global GATES
 GATES+=1
 if not ok:raise ArithmeticError(why)

def eligible(p,q):
 if not(0<p<q and gcd(p,q)==1 and gcd(p*q,6)==1):return False
 n=p+q;d=2
 while d*d<=n:
  e=0
  while n%d==0:n//=d;e+=1
  if e and (d%3!=2 or e>2):return False
  d+=1
 return n==1 or n%3==2

def intervals(p,q):
 # Literal strict intersections, no merge at excluded endpoints.
 A=[(max(F(0),F(7*k-1,7*p)),min(F(1),F(7*k+1,7*p)),k%2) for k in range(p+1)]
 B=[(max(F(0),F(7*k-1,7*q)),min(F(1),F(7*k+1,7*q)),k%2) for k in range(q+1)]
 out=[];i=j=0
 while i<len(A) and j<len(B):
  x=max(A[i][0],B[j][0]);y=min(A[i][1],B[j][1])
  if x<y and A[i][2]!=B[j][2]:out.append((x,y))
  if A[i][1]<B[j][1]:i+=1
  else:j+=1
 return out

def widths(p,q):return sorted([b-a for a,b in intervals(p,q)],reverse=True)
def top(w,k):return sum(w[:k],F(0))
def weak(u,v):return all(top(u,k)<=top(v,k) for k in range(1,max(len(u),len(v))+1))
def split(w,g):return [a/g for a in w for _ in range(g)]
def escapes(u,v):return any(0<top(u,k)>=top(v,k) for k in range(1,len(u)+1))

FRONT={
 (7,13):[F(1,49)]*2,
 (19,25):[F(37,3325)]*2+[F(23,3325)]*2+[F(9,3325)]*2,
 (5,41):[F(2,287)]*6,
 (5,53):[F(2,371)]*6+[F(9,1855)]*2,
 (1,67):[F(2,469)]*10,
}

def main():
 head={(p,q):widths(p,q) for q in range(2,68) for p in range(1,q) if eligible(p,q)}
 need(len(head)==46,'complete q<=67 head')
 maximal={a for a,u in head.items() if not any(a!=b and weak(u,v) and u!=v for b,v in head.items())}
 need(maximal==set(FRONT),'exact ranked frontier')
 for pair,v in FRONT.items():need(head[pair]==v,'whole component vector')
 for pair,u in head.items():
  witnesses=[v for v in FRONT.values() if weak(u,v)]
  need(bool(witnesses),'head domination')
  for g in (1,3,5):
   ug=split(u,g)
   need(weak(ug,u),'dilation refines lengths')
   for v in witnesses:need(weak(ug,split(v,g)),'dilation preserves dominance')
   for k in range(1,len(ug)+2):
    r,s=divmod(k,g)
    need(top(ug,k)==top(u,r)+(F(s,g)*u[r] if r<len(u) else 0),'interpolation formula')
 # The all-height tail proof uses the inherited mass cap, plus the width bound.
 need(F(2,7*69)<F(2,469),'odd q>67 width tail')
 for q in range(69,126,2):
  for p in range(1,q,2):
   if eligible(p,q):
    u=widths(p,q)
    need(sum(u)<F(20,469) and max(u,default=0)<F(2,469),'finite tail controls')
    need(weak(u,FRONT[1,67]),'tail domination control')
 # Abstract gain, explicitly not a physical eleven-body assertion.
 u=[F(1,100)]*4+[F(1,2000)]*2
 need(sum(u)==F(41,1000),'abstract mass')
 need(sum(u)<sum(FRONT[19,25]) and max(u)<max(FRONT[19,25]),'old coupled profile misses')
 need(top(u,4)==F(1,25)>F(120,3325),'new rank-four escape')
 need(all(escapes(u,v) for v in FRONT.values()),'all ranked profile gates')
 # Rank one and final rank alone really lose the intermediate coordinate.
 need(not weak(u,FRONT[19,25]),'endpoint-only hostile')
 # A single component can contain several closed body components. Do not
 # impose an injection between the two component sets.
 need(weak([F(2,5),F(2,5)],[F(9,10)]),'many body pieces may share one carrier')
 need(escapes([F(9,10)],[F(9,10)]),'closed/open equality boundary')
 need(not escapes([],[F(9,10)]),'isolated-only set has no length certificate')
 # Exhaust every paired list on a declared rational alphabet, padding zeros.
 values=[F(0),F(1,2000),F(1,200),F(7,1000),F(1,100),F(11,1000),F(1,50)]
 paired=0
 from itertools import combinations_with_replacement
 for half in combinations_with_replacement(values,5):
  z=sorted([a for a in half for _ in range(2) if a],reverse=True)
  M=sum(z,F(0));L=max(z,default=F(0))
  for pair,v in FRONT.items():
   reduced=(M>=sum(v) or L>=v[0] or (pair==(19,25) and top(z,4)>=F(120,3325)))
   need(escapes(z,v)==reduced,'paired breakpoint characterization')
  paired+=1
 # Validate the abstract strict containment inequality for all finite
 # subsets of explicit closed subpieces inside three open carrier arcs.
 capacities=[F(1,5),F(1,7),F(1,11)]
 pieces=[(j,c/8) for j,c in enumerate(capacities) for _ in range(3)]
 for mask in range(1,1<<len(pieces)):
  z=sorted([pieces[i][1] for i in range(len(pieces)) if mask>>i&1],reverse=True)
  for k in range(1,len(z)+1):need(top(z,k)<top(capacities,k),'strict packing control')
 print('STATUS: ranked all-height inert envelope; five complete maximal profiles')
 print('FINITE_HEAD: 46 eligible pairs q<=67; complete partial-sum comparisons')
 for pair,v in FRONT.items():print('FRONTIER',pair,'widths',','.join(map(str,v)))
 print('TAIL_PROOF: q>67 width<2/469 and inherited mass<20/469')
 print('DILATIONS: g=1,3,5; exact interpolation and dominance controls')
 print('ABSTRACT_GAIN: mass41/1000 max1/100 top4=1/25; physical realization NOT claimed')
 print('PAIRED_LIST_CONTROLS:',paired,'declared seven-value alphabet, ten entries with zeros')
 print('HOSTILES: endpoint-only loss; many-to-one containment; equality; isolated-only geometry')
 print('SCOPE: sufficient physical-body certificate; no new eleven-body family established')
 print('ACTIVE_GATES:',GATES)

if __name__=='__main__':main()
