from functools import lru_cache
from fractions import Fraction
from itertools import combinations,permutations
import json,sys

def verify(p):
 n=len(p);full=(1<<n)-1
 def require(x,msg):
  if not x:raise RuntimeError(msg)
 require(len(p)==n,'size')
 for b in range(n):
  require(p[b]<(1<<b),'natural labeling')
  for a in range(n):
   if p[b]>>a&1:require(not p[a]&~p[b],'transitivity')
 u=[sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n)]
 deg=[n-1-(p[i]|u[i]).bit_count() for i in range(n)]
 h=[]
 for i in range(n):h.append(1+max([h[j] for j in range(i) if p[i]>>j&1],default=0))
 require(max(deg)<=7 and max(deg)>0,'thin non-total');require(max(h)<=3,'height')
 def chain(m,d):return all(d[a]>>b&1 or d[b]>>a&1 for a,b in combinations([i for i in range(n) if m>>i&1],2))
 very=[];eligible=[]
 for dual,(d,v) in enumerate(((p,u),(u,p))):
  for a,b in permutations(range(n),2):
   if (d[a]|v[a])>>b&1:continue
   if d[a]==d[b] and chain(v[a]&~v[b],d) and chain(v[b]&~v[a],d):very.append([dual,a,b])
   if not d[a]&~d[b] and chain(v[b]&~v[a],d):eligible.append((dual,a,b))
 require(not very,'very good pair found')
 @lru_cache(None)
 def f(mask):
  if mask==full:return 1
  return sum(f(mask|1<<i) for i in range(n) if not mask>>i&1 and not p[i]&~mask)
 E=f(0);N=[[0]*n for _ in range(n)];forward={0:1}
 for mask in range(full+1):
  z=forward.get(mask,0)
  if not z:continue
  for b in range(n):
   if mask>>b&1 or p[b]&~mask:continue
   nxt=mask|1<<b;forward[nxt]=forward.get(nxt,0)+z
   for a in range(n):
    if mask>>a&1:N[a][b]+=z*f(nxt)
 require(forward[full]==E,'DP count')
 entries=[{'dual':d,'a':a,'b':b,'numerator':N[b][a] if d else N[a][b],'probability':str(Fraction(N[b][a] if d else N[a][b],E))} for d,a,b in eligible]
 good=[r for r in entries if 2*r['numerator']<=E]
 half=[[a,b] for a,b in combinations(range(n),2) if 2*N[a][b]==E]
 require(not good,'Good pair found')
 require(len(set((p[i].bit_count(),u[i].bit_count()) for i in range(n)))==n or n!=10,'Rigidity signature')
 out={'status':'EXACT_HEIGHT_THREE_NO_VERY_GOOD_WITNESS_VERIFIED','pred':p,'height':max(h),'incomparability_degrees':deg,'extensions':E,'very_good_pairs':0,'good_pairs':good,'eligible':entries,'half_pairs':half}
 return out

def enumerate_gold(p, first, second):
 n=len(p);full=(1<<n)-1;cells=[0,0,0,0];total=0
 def visit(word,mask):
  nonlocal total
  if mask==full:
   total+=1;r={v:i for i,v in enumerate(word)}
   a=r[first[0]]<r[first[1]];b=r[second[0]]<r[second[1]]
   cells[(0 if a else 2)+(0 if b else 1)]+=1;return
  for v in range(n):
   if not mask>>v&1 and not p[v]&~mask:visit(word+[v],mask|1<<v)
 visit([],0);t=cells[0]+cells[1]
 margins=[total-t-cells[0],total-t-cells[1],t-cells[2],t-cells[3]]
 if min(margins)<0:raise RuntimeError('Gold failed')
 return {'extensions':total,'first_comparison':first,'second_comparison':second,'four_cells_11_10_01_00':cells,'gold_margins':margins}

def main():
 examples=[(0,0,0,6,0,19,2,18,86,151),(0,0,0,0,7,13,14,12,206,5,685,535)]
 results=[]
 for p in examples:
  out=verify(p)
  gold=enumerate_gold(p,[3,5],[2,4]) if len(p)==10 else enumerate_gold(p,[0,3],[0,3])
  if gold['extensions']!=out['extensions']:raise RuntimeError('Enumeration/DP mismatch')
  if len(p)==10:
   if out['extensions']!=11971 or gold['four_cells_11_10_01_00']!=[3281,3281,1899,3510]:raise RuntimeError('Ten-point regression')
   if len(out['eligible'])!=48 or min(r['numerator'] for r in out['eligible'])!=6547:raise RuntimeError('Eligibility regression')
  elif out['extensions']!=61180:raise RuntimeError('Twelve-point regression')
  out['gold']=gold;results.append(out)
 print(json.dumps({'status':'BOTH_HEIGHT_THREE_NO_GOOD_WITNESSES_AND_GOLD_VERIFIED','examples':results},indent=2))
if __name__=='__main__':main()
