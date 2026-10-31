import itertools,json
from check import direct

def one_ear(word):
 n=len(word)+1;ad=[0]*n
 for i,c in enumerate(word,1):
  if c==2:ad[0]|=1<<i
  else:ad[i]|=1
  for j in range(1,i):ad[j]|=1<<i
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   nxt=ad[i]&~mask
   while nxt:
    bit=nxt&-nxt;nxt-=bit;j=bit.bit_length()-1
    dp[mask|bit][j]+=dp[mask][i]
 return sum(dp[-1]),dp[-1][0]
checks=0
for n in range(7):
 for word in itertools.product((0,2),repeat=n):
  s,t=one_ear(word);v=1
  for mask in range(1,1<<n):
   first=(mask&-mask).bit_length()-1
   if word[first]==2:v+=1
  if direct(word)!=s-t+v:raise RuntimeError(('decomposition',word))
  checks+=1
  for k in range(1,4):
   if direct(word+(2,)*k)!=s+(1<<k)*v:raise RuntimeError(('iteration',word,k))
   checks+=1
  if one_ear(word+(2,0))[0]<3:raise RuntimeError(('triangle',word))
  checks+=1
for k in range(2,6):
 if direct((1,2)*k)<3**(2*k-4):raise RuntimeError(('mixedgrowth',k))
 checks+=1
print(json.dumps(dict(identity_checks=checks,max_seed_length=6,optimized_compatible=True)))
