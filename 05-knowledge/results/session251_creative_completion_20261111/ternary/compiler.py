#!/usr/bin/env python3
"""Fixed-core attachment compiler; actual tournaments, exact integer arithmetic."""
import itertools,json,functools,hashlib
from pathlib import Path

def hp_subsets(rows):
 n=len(rows); full=(1<<n)-1
 dp=[[0]*n for _ in range(full+1)]
 for v in range(n): dp[1<<v][v]=1
 for mask in range(1,full+1):
  for v in range(n):
   x=dp[mask][v]
   if not x:continue
   avail=rows[v]&~mask&full
   while avail:
    b=avail&-avail;avail-=b;dp[mask|b][b.bit_length()-1]+=x
 h=[sum(x) for x in dp];h[0]=1
 return h

def signature(rows,O):
 n=len(rows);full=(1<<n)-1;h=hp_subsets(rows);p=[0]*8
 def walk(mask,start,end):
  w=h[full^mask];sgn=(-1)**mask.bit_count()
  for cat,yes in enumerate((True,not(O>>end&1),bool(O>>start&1),bool(O>>start&1) and not(O>>end&1))):
   if yes:p[2*cat]+=w;p[2*cat+1]+=sgn*w
  avail=rows[end]&~mask&full
  while avail:
   b=avail&-avail;avail-=b;walk(mask|b,start,b.bit_length()-1)
 for v in range(n):walk(1<<v,v,v)
 return (h[full],*p)

@functools.lru_cache(None)
def g(a):
 if a<2:return 1
 q=a//2
 return g(q) if a&1 else g(q)+q+(q&1)

def count(sig,a,b):
 d,kp,km,ip,im,op,om,jp,jm=sig
 x,y=a&1,b&1
 return d*(g(a)+g(b)-1)+(d+kp)*a*b+(d+km)*x*y+ip*a-im*x+op*b-om*y+jp+jm

def construct(rows,O,a,b):
 left=list(map(int,bin(a)[2:])) if a else []
 # B is complement of right word read in reverse.
 right=[1-int(x) for x in bin(b)[2:][::-1]] if b else []
 l=len(left);n=len(rows);r=len(right);old=l+n+r
 out=[0]*(old+1)
 for i in range(old):
  for j in range(i+1,old):
   if l<=i<j<l+n:win=bool(rows[i-l]>>(j-l)&1)
   else:win=True
   u,v=(i,j) if win else (j,i);out[u]|=1<<v
 bits=left+[(O>>v)&1 for v in range(n)]+right
 for v,bit in enumerate(bits):out[old if bit else v]|=1<<(v if bit else old)
 return out

def generate(sig,bound):
 d,kp,km,ip,im,op,om,jp,jm=sig
 m=min((d+km)*x*y-im*x-om*y+jp+jm for x in (0,1) for y in (0,1))
 cap=bound-d-m
 if cap<0:return
 # Uniform contact has a degenerate axis d*g(t); reduce t to even 2q.
 if not op:
  for q in range(bound//d):
   v=d*g(2*q)
   if v<=bound:yield v,0,2*q
 if not ip:
  for q in range(bound//d):
   v=d*g(2*q)
   if v<=bound:yield v,2*q,0
 maxa=cap//ip if ip else cap//(d+kp)
 for a in range(maxa+1):
  denom=(d+kp)*a+op
  if not denom:continue
  for b in range((cap-ip*a)//denom+1):
   if (not ip and b==0) or (not op and a==0):continue
   value=count(sig,a,b)
   if value<=bound:yield value,a,b

def tournament(n,code):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   u,v=(i,j) if code>>k&1 else (j,i);r[u]|=1<<v;k+=1
 return r

def canonical(rows):
 n=len(rows)
 return min(tuple(sum(((rows[p[i]]>>p[j])&1)<<j for j in range(n)) for i in range(n)) for p in itertools.permutations(range(n)))

def bank(maxn=5):
 # Full labelled enumeration then isomorphism quotient. No converse merge.
 ans=[]
 for n in range(1,maxn+1):
  classes=sorted({canonical(tournament(n,c)) for c in range(1<<(n*(n-1)//2))})
  ans.extend(classes)
 return ans

def main():
 cores=bank();sigs={};tests=0
 for rows in cores:
  for O in range(1<<len(rows)):
   sig=signature(rows,O);sigs.setdefault(sig,(rows,O))
   for a,b in itertools.product(range(4),repeat=2):
    child=construct(rows,O,a,b)
    actual=hp_subsets(child)[-1]
    expected=count(sig,a,b)
    if actual!=expected:raise AssertionError((rows,O,a,b,actual,expected))
    tests+=1
 bound=4095;witness={};by_parent={}
 for sig,(rows,O) in sorted(sigs.items()):
  d=sig[0]
  for val,a,b in generate(sig,bound):
   w={'rows':list(rows),'out_set':O,'A':a,'B':b,'signature':list(sig)}
   witness.setdefault(val,w);by_parent.setdefault(d,set()).add(val)
 # Include complete one-apex family up to target. Every g-value has a reduced even input 2q <=2(bound-1).
 apex={g(2*q) for q in range(bound)}
 groups={}
 cumulative=set(apex)
 for d in sorted(by_parent):
  cumulative|=by_parent[d]
  missing=[h for h in range(1,bound+1,2) if h not in cumulative]
  groups[str(d)]={'first_misses':missing[:30],'missing_count':len(missing),'initial_allowed_interval':next((h-2 for h in missing if h not in (7,21)),bound)}
 selected={str(v):witness[v] for v in (55,95,107,109,127,265,281,287,289) if v in witness}
 for v,w in selected.items():
  child=construct(w['rows'],w['out_set'],w['A'],w['B']);w['child_rows']=child
  if hp_subsets(child)[-1]!=int(v):raise AssertionError(v)
 report={'scope':'All cuts of every core class of orders <=5, arbitrary transitive flanks, plus one-apex graphs, target <=4095. Complete target inversion, not whole tournament spectrum.', 'core_classes_by_order':{str(n):sum(len(c)==n for c in cores) for n in range(1,6)},'distinct_signatures':len(sigs),'direct_dp_tests':tests,'cumulative_parent_banks':groups,'selected_witnesses':selected}
 Path(__file__).with_name('results.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
