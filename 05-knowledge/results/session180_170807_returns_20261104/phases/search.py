import itertools,json
from pathlib import Path

def response(rows):
 n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)];hs=[0]*N;hs[0]=1
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    dp[mask][v]=1 if rest==0 else sum(dp[rest][u] for u in range(n) if rows[u]>>v&1)
  hs[mask]=sum(dp[mask])
 F=[sum(dp[s][w]*hs[(N-1)^s] for s in range(N)) for w in range(n)]
 return hs[-1],F

def strong(rows):
 n=len(rows)
 for reverse in [False,True]:
  seen=1
  while True:
   new=seen
   for v in range(n):
    if seen>>v&1:
     new|=sum(1<<u for u in range(n) if (rows[u]>>v&1)) if reverse else rows[v]
   if new==seen:break
   seen=new
  if seen!=(1<<n)-1:return False
 return True

def main():
 bank={};orders={}
 for n in range(1,7):
  pairs=list(itertools.combinations(range(n),2));strongn=0
  for code in range(1<<len(pairs)):
   rows=[0]*n
   for k,(u,v) in enumerate(pairs):
    if code>>k&1:rows[u]|=1<<v
    else: rows[v]|=1<<u
   if not strong(rows):continue
   strongn+=1;h,fs=response(rows)
   for w,f in enumerate(fs):bank.setdefault((h,f),{'rows':rows[:],'port':w})
  orders[n]=strongn;print(n,strongn,len(bank),flush=True)
 matches={str(d):[] for d in [2,4]}
 for (h,f),w in bank.items():
  for d in [2,4]:
   if (3*h+d,3*f) in bank:matches[str(d)].append({'input':[h,f], 'output':[3*h+d,3*f], 'source':w, 'target':bank[(3*h+d,3*f)]})
 Path(__file__).with_name('small_suffix_search.json').write_text(json.dumps({'orders':orders,'bank':[{'h':h,'f':f,**w} for (h,f),w in sorted(bank.items())],'matches':matches},indent=2))
 print('matches', {d:len(ms) for d,ms in matches.items()});print(json.dumps(matches,indent=2))
if __name__=='__main__':main()
