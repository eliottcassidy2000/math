from pathlib import Path
rows=[tuple(map(int,l.split())) for l in (Path(__file__).resolve().parent.parent / 'certificates_interval' / 'witnesses.txt').read_text().splitlines()]
if [r[0] for r in rows]!=list(range(80407,80875,2)): raise ValueError("interval coverage")
for h,code,cut in rows:
 n=12; out=[0]*n; k=0
 for i in range(11):
  for j in range(i+1,11):
   if code>>k&1:out[i]|=1<<j
   else:out[j]|=1<<i
   k+=1
 for i in range(11):
  if cut>>i&1:out[i]|=1<<11
  else:out[11]|=1<<i
 if code>>55 or cut>>11: raise ValueError("out-of-range code")
 dp=[[0]*12 for _ in range(1<<12)]
 for i in range(12):dp[1<<i][i]=1
 for s in range(1,1<<12):
  for v in range(12):
   z=out[v]&~s
   while z:
    b=z&-z;w=b.bit_length()-1;z-=b
    dp[s|b][w]+=dp[s][v]
 got=sum(dp[-1])
 if got!=h: raise ValueError((h,got))
 for v in range(12):
  seen=1<<v;before=0
  while seen!=before:
   before=seen
   for w in range(12):
    if seen>>w&1:seen|=out[w]
  if seen!=(1<<12)-1: raise ValueError((h,v))
print('PASS',len(rows),'order-12 strong, direct subset-DP and reachability; inclusive',rows[0][0],rows[-1][0])
