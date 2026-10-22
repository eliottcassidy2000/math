import itertools as it, collections as co, json, math
from pathlib import Path
OUT=Path(__file__).parent

def require(c,m):
 if not c: raise RuntimeError(m)

def edges(n): return list(it.combinations(range(n),2))
def arc(mask,n,u,v):
 if u>v: return 1-arc(mask,n,v,u)
 return (mask>>edges(n).index((u,v)))&1

def induce(mask,n,vs):
 return sum(arc(mask,n,vs[i],vs[j])<<b for b,(i,j) in enumerate(edges(len(vs))))

def classes(n):
 E=edges(n); ix={e:i for i,e in enumerate(E)}; lookup={}; reps=[]
 perms=list(it.permutations(range(n)))
 maps=[[(ix[tuple(sorted((p[i],p[j])))],int(p[i]>p[j])) for i,j in E] for p in perms]
 for m in range(1<<len(E)):
  if m in lookup: continue
  reps.append(m)
  for mp in maps:
   x=sum((((m>>a)&1)^f)<<b for b,(a,f) in enumerate(mp))
   lookup[x]=m
 return reps,lookup
C={}; L={}
for n in range(1,7): C[n],L[n]=classes(n)
def profile(m,n,k): return tuple(sorted(co.Counter(L[k][induce(m,n,vs)] for vs in it.combinations(range(n),k)).items()))
def score(m,n): return tuple(sorted(sum(arc(m,n,u,v) for v in range(n) if u!=v) for u in range(n)))
def hp(m,n):
 dp={(1<<v,v):1 for v in range(n)}
 for s in range(1,1<<n):
  for v in range(n):
   x=dp.get((s,v),0)
   if not x:continue
   for w in range(n):
    if not (s>>w)&1 and arc(m,n,v,w):dp[s|1<<w,w]=dp.get((s|1<<w,w),0)+x
 return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
def cycsets(m,n): return [frozenset(vs) for vs in it.combinations(range(n),3) if sum(arc(m,n,vs[i],vs[(i+1)%3]) for i in range(3)) in (0,3)]
def overlap(m,n,r):
 cs=cycsets(m,n)
 return sum(len(a&b)==r for a,b in it.combinations(cs,2))
result={'class_counts':{n:len(C[n]) for n in C},'levels':{},'witnesses':{}}
require(list(result['class_counts'].values())==[1,1,2,4,12,56],'A000568 n<=6')
for n in range(3,7):
 rows=[]
 for m in C[n]:
  s=score(m,n); c=len(cycsets(m,n)); require(c==math.comb(n,3)-sum(math.comb(d,2) for d in s),'triangle identity')
  for j in range(n):
   star=sum(math.comb(d,j) for d in s)
   rhs=sum(mult*sum(sum(arc(h,j+1,u,v) for v in range(j+1) if u!=v)==j for u in range(j+1)) for h,mult in profile(m,n,j+1))
   require(star==rhs,'star moments')
  rows.append({'mask':m,'scores':s,'H':hp(m,n),'D3':profile(m,n,3),'proper':profile(m,n,n-1),'overlap1':overlap(m,n,1),'overlap2':overlap(m,n,2)})
 result['levels'][n]={'score_classes':len(set(r['scores'] for r in rows)),'Dk_classes':{k:len(set(profile(m,n,k) for m in C[n])) for k in range(1,n+1)},'scores_plus_Dk_classes':{k:len(set((score(m,n),profile(m,n,k)) for m in C[n])) for k in range(1,n+1)}}
 for a,b in it.combinations(rows,2):
  tests={'same_scores_different_H':a['scores']==b['scores'] and a['H']!=b['H'],'same_proper_deck_not_converse':a['proper']==b['proper'] and L[n][a['mask']^((1<<math.comb(n,2))-1)]!=b['mask'],'same_D3_different_overlap2':a['D3']==b['D3'] and a['overlap2']!=b['overlap2'],'same_scores_proper_deck_different_H':a['scores']==b['scores'] and a['proper']==b['proper'] and a['H']!=b['H']}
  for key,ok in tests.items():
   if ok and key not in result['witnesses']:result['witnesses'][key]={'n':n,'a':a,'b':b}
 # Verify overlap reconstruction from union decks for r1,2.
 for m in C[n]:
  for r in (1,2):
   k=6-r
   if n>=k: require(overlap(m,n,r)==sum(mult*overlap(h,k,r) for h,mult in profile(m,n,k)),'overlap union identity')
print(json.dumps(result,indent=2)); (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
