import importlib.util,itertools,collections,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
spec=importlib.util.spec_from_file_location('compiler',ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
def canonical(rows):
 groups=[]
 for _, vs in itertools.groupby(sorted(range(len(rows)),key=lambda v:rows[v].bit_count()),key=lambda v:rows[v].bit_count()):groups.append(list(vs))
 return min(tuple(sum(((rows[p[i]]>>p[j])&1)<<j for j in range(len(rows))) for i in range(len(rows))) for gs in itertools.product(*[itertools.permutations(g) for g in groups]) for p in [sum(gs,())])
def expand(old):
 ans=set()
 for row in old:
  n=len(row)
  for out in range(1<<n):
   r=list(row)+[out]
   for v in range(n):
    if not(out>>v&1):r[v]|=1<<n
   ans.add(canonical(r))
 return sorted(ans)
def cuts(rows):
 n=len(rows);full=(1<<n)-1;h=c.hp_subsets(rows);wp=[[0]*n for _ in rows];wm=[[0]*n for _ in rows]
 def walk(mask,first,last):
  w=h[full^mask];wp[first][last]+=w;wm[first][last]+=(-1)**mask.bit_count()*w
  todo=rows[last]&~mask&full
  while todo:
   b=todo&-todo;todo-=b;walk(mask|b,first,b.bit_length()-1)
 for v in range(n):walk(1<<v,v,v)
 kp=sum(map(sum,wp));km=sum(map(sum,wm))
 for O in range(1<<n):
  ip=sum(wp[a][b] for a in range(n) for b in range(n) if not(O>>b&1));im=sum(wm[a][b] for a in range(n) for b in range(n) if not(O>>b&1))
  op=sum(wp[a][b] for a in range(n) for b in range(n) if O>>a&1);om=sum(wm[a][b] for a in range(n) for b in range(n) if O>>a&1)
  jp=sum(wp[a][b] for a in range(n) for b in range(n) if (O>>a&1) and not(O>>b&1));jm=sum(wm[a][b] for a in range(n) for b in range(n) if (O>>a&1) and not(O>>b&1))
  yield O,(h[-1],kp,km,ip,im,op,om,jp,jm)
def main(maxn=7):
 old=[(0,)];groups=collections.defaultdict(dict);out=[]
 for n in range(1,maxn+1):
  if n>1:old=expand(old)
  found=[]
  for rows in old:
   for O,s in cuts(rows):
    key=s[:7];q=s[7]+s[8]
    for q2,w in groups[key].items():
     if q2!=q:found.append({'difference':q-q2,'first':w,'second':{'rows':rows,'O':O,'signature':s}})
    groups[key].setdefault(q,{'rows':rows,'O':O,'signature':s})
  units=[x for x in found if abs(x['difference'])==2]
  out.append({'n':n,'classes':len(old),'distinct_central_collisions':len(found),'unit_collisions':len(units),'first_units':units[:3],'first_any':found[:1]})
  print(json.dumps(out[-1]),flush=True)
  # Continue through the declared order.
 (HERE/'central_search.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
