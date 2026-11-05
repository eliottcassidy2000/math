"""Independent U-ordered profile/cost-bitmask census and height ledger."""
import heapq,json
from math import isqrt
from itertools import product

def need(q,m):
 if not q:raise RuntimeError(m)
T=(1,0,2,1)
def step(s,a):
 D,E,U,d=s
 if a=='O':return D,0,2*U,d
 if a=='Z':return D-E+U,d,2*U-d,d
 if a=='J':return D-E+U,D-E,D+U,D
 if a=='I':return s
 raise ValueError(a)
def G(s):return s[0]+2*s[2]
def word(w):
 s=T
 for a in w:s=step(s,a)
 return s
cap=50001;width=32;masklimit=(1<<(4*width))-1
costs={T:1};heap=[(T[2],T)];support=set();unrestricted=set();triples=0;hostile=[];ledger=0
while heap:
 _,s=heapq.heappop(heap);bits=costs[s];D,E,U,d=s;n=G(s);W=U-d
 unrestricted.add(n);tri=bits
 while tri:
  lsb=tri&-tri;idx=lsb.bit_length()-1;tri-=lsb;r,k=divmod(idx,width)
  need(W>=2**k and k<=(n-3).bit_length()-2,'charged ledger')
  need(k+r+4<=(n-3).bit_length()-1+r+3,'graph order')
  if k in (n.bit_length()-3,n.bit_length()-2):support.add(n)
  if n==49937:hostile.append({'profile':s,'J_count':r,'retained_count':k})
  triples+=1;ledger+=1
 for a in 'OZJ':
  q=step(s,a);newbits=(bits<<(width if a=='J' else 1))&masklimit
  if not newbits or G(q)>cap:continue
  need(q[2]>U and G(q)>n,'prefix pruning')
  old=costs.get(q,0)
  costs[q]=old|newbits
  if old==0:heapq.heappush(heap,(q[2],q))
missing=[n for n in range(3,cap+1,2) if n not in support]
need(triples==641253 and len(support)==23502 and len(missing)==1498,'full census totals')
need(max(missing)==49937 and 49937 in unrestricted,'largest hostile')
need(all(x['retained_count'] not in (13,14) for x in hostile),'no two-top lengths')
need(word('OOOZOZOJOOZOZJ')==(16593,6588,16672,6667),'explicit hostile source')
need(G(word('OOOZOZOJOOZOZJ'))==49937,'positive hostile target')
need(not [n for n in range(27783,50002,2) if n not in unrestricted],'full three-root finite coverage')
# Exact ledger maps on an unrestricted signed basis.
for s in [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-3,5,7,-1)]:
 D,E,U,d=s;W=U-d
 for a in 'OZJ':
  DD,EE,UU,dd=step(s,a);expected={'O':(2*W+d,d),'Z':(2*W,d),'J':(W+d,D)}[a]
  need((UU-dd,dd)==expected,'ledger map')
for k in range(1,13):
 s=word('Z'*k);n=G(s)
 need(s==(2**k+1,1,2**k+1,1),'sharp profile')
 need(k+4==(n-3).bit_length()-1+3,'sharp graph order')
need(G(T)==5 and 4==(5-3).bit_length()-1+3,'empty sharp boundary')
# Maximal-letter bounds, root-budget field predicate, and large-prime scale.
word_controls=0
for length in range(8):
 for letters in product('OZJ',repeat=length):
  s=word(letters)
  need(s[2]<=2**(length+1) and G(s)<=3*2**(length+1),'finite word height')
  word_controls+=1
for n in range(5,1001):
 j=(n-3).bit_length()-1;L=j+2;bound=24*(n-3)
 need(3*2**(L+1)<=bound and n<bound and L<bound,'prime threshold')
 # For any prime p>bound, the possible ordinary J counts are below p.
 for c in range(L+1):need((c*(c-1)*(c-2)*(c-3)==0)==(c<=3),'root count polynomial over integers')
print(json.dumps({'status':'PASS','cap':cap,'cost_triples':triples,'physical_profiles':len(costs),'two_length_support':len(support),'two_length_missing':len(missing),'largest_two_length_gap':max(missing),'hostile49937':sorted(hostile,key=lambda z:(z['profile'],z['J_count'],z['retained_count'])),'ledger_checks':ledger,'short_word_upper_bound_checks':word_controls,'prime_scale_n_through':1000,'full_three_root_coverage_above27781_to50001':True},sort_keys=True,indent=2))
