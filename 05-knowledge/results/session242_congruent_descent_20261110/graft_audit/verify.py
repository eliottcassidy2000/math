#!/usr/bin/env python3
"""Independent recursive first-vertex counts and literal Hamilton permutations."""
import functools,itertools,json,math,sys

def need(x,s):
 if not x:raise RuntimeError(s)

def stats(a,u=0):
 n=len(a);top=(1<<n)-1
 need(all(0<=r<=top and not(r>>i&1) for i,r in enumerate(a)), 'shape')
 need(all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i in range(n) for j in range(i)), 'orientation')
 @functools.lru_cache(None)
 def start(m,v):
  rest=m^(1<<v)
  if not rest:return 1
  return sum(start(rest,w) for w in range(n) if rest>>w&1 and a[v]>>w&1)
 def h(m):return sum(start(m,v) for v in range(n) if m>>v&1) if m else 1
 hs=[h(m) for m in range(top+1)]
 q=sum(hs[m]*hs[top^m] for m in range(top+1))
 U=sum(start(m,u)*hs[top^m] for m in range(top+1) if m>>u&1)
 # End at u: the preceding path starts anywhere and ends in a predecessor;
 # reversing edges turns this into a first-vertex count.
 reverse=[top^(1<<v)^r for v,r in enumerate(a)]
 @functools.lru_cache(None)
 def end(m,v):
  rest=m^(1<<v)
  return 1 if not rest else sum(end(rest,w) for w in range(n) if rest>>w&1 and reverse[v]>>w&1)
 return hs[-1],q,[hs[-1],end(top,u),U,hs[top^(1<<u)]]

def graft(a,u,b):
 n=len(a);out=list(a)+[r<<n for r in b]
 for i in range(n):
  for j in range(len(b)):
   tail,head=(n+j,i) if i==u else(i,n+j)
   out[tail]|=1<<head
 return out

def brute(a):
 return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)) for p in itertools.permutations(range(len(a))))

def strong(a):
 for v in range(len(a)):
  seen={v};todo=[v]
  while todo:
   x=todo.pop()
   for w in range(len(a)):
    if a[x]>>w&1 and w not in seen:seen.add(w);todo.append(w)
  if len(seen)<len(a):return False
 return True

G=[22,60,8,81,4,29,55]; Z=[126,180,120,2,168,72,154,45]
W=graft([2,4,1],0,G)
S=[r|(1<<len(Z)) for r in Z]+[0]
V=graft(S,3,[0])
need(stats(G)[:2]==(49,812),'gadget pair')
need(stats(graft([0],0,G))[0]//7==7,'improper singleton-root counterexample')
need(stats(Z,3)[2]==[259,70,266,77],'root state')
receipts=[]
for label,a,target in [('gadget',G,49),('root',Z,259),('C3_graft',W,1673),('SB',V,791)]:
 counted=stats(a)[0]; literal=brute(a)
 need(counted==literal==target,label+' count')
 need(strong(a),label+' strong')
 receipts.append(dict(label=label,rows=a,vertices=len(a),recursive_count=counted,permutation_count=literal,strong=True))
# Complete labelled source universe through 4 vertices, gadgets through 3.
def allrows(n):
 pairs=list(itertools.combinations(range(n),2))
 for code in range(1<<len(pairs)):
  a=[0]*n
  for k,(i,j) in enumerate(pairs):
   x,y=(i,j) if code>>k&1 else(j,i);a[x]|=1<<y
  yield a
bank=[a for n in range(1,4) for a in allrows(n)]
checks=0;sinks=0
for n in range(1,5):
 for a in allrows(n):
  for u in range(n):
   D,E,U,d=stats(a,u)[2]
   for b in bank:
    H,Q,_=stats(b);N=stats(graft(a,u,b))[0]
    need(N==H*(D-E+U)+(Q-2*H)*(U-d),'scalar row')
    if a[u]:
     need(N>=H+Q and N//math.gcd(H,Q)!=7,'proper bound/exclusion');checks+=1
    else:
     need(N==H*d,'sink boundary');sinks+=1
# Singleton gadget is the boundary omitted by divisibility argument.
singleton_pairs=[stats([2,4,1],0)[:2],stats([2,0],0)[:2],stats([0],0)[:2]]
need(singleton_pairs==[(3,12),(1,4),(1,2)],'boundary controls')
# All integer normalized pairs allowed by necessary bounds at quotient7.
forced=[(h,q) for h in range(1,8,2) for q in range(4,8,4) if q>=2*h and h+q<=7]
need(forced==[(1,4)],'forced pair')
if '--fail' in sys.argv:need(False,'deliberate hostile control')
print(json.dumps(dict(status='PASS',witnesses=receipts,proper_literal_checks=checks,sink_literal_checks=sinks,normalized_pairs_at_7=forced,root_state_D_E_U_d=stats(Z,3)[2],gadget_gcd=math.gcd(*stats(G)[:2])),indent=2))
