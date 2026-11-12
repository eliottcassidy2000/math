"""Independent backward subset recurrence for the finite receiver-repair seeds."""
import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def join(L,R):
 n=len(L);m=len(R);rows=[0]*(n+m-1);maps=(list(range(n)),[0]+list(range(n,n+m-1)))
 for T,ids in zip((L,R),maps):
  for i,j in ((i,j) for i in range(len(T)) for j in range(len(T))):
   if T[i]&(1<<j):rows[ids[i]]|=1<<ids[j]
 for i in range(1,n):
  for j in range(n,n+m-1):rows[i]|=1<<j
 return tuple(rows)
def apex(word):
 rows=[0]*(len(word)+1)
 for i,c in enumerate(word,1):
  if c=='1':rows[0]|=1<<i
  else:rows[i]|=1
  for j in range(i+1,len(rows)):rows[i]|=1<<j
 return tuple(rows)
def profile(T,full_profile=True):
 n=len(T);size=1<<n;pred=[sum(1<<j for j in range(n) if T[j]&(1<<i)) for i in range(n)]
 need(all(not(T[i]&(1<<i)) and all(bool(T[i]&(1<<j))!=bool(T[j]&(1<<i)) for j in range(i)) for i in range(n)),'tournament')
 dp=[[0]*n for _ in range(size)];H=[0]*size;H[0]=1
 for S in range(1,size):
  bits=S
  while bits:
   bit=bits&-bits;bits-=bit;v=bit.bit_length()-1;rest=S^bit
   if not rest:dp[S][v]=1
   else:
    incoming=pred[v]&rest;total=0
    while incoming:
     p=incoming&-incoming;incoming-=p;total+=dp[rest][p.bit_length()-1]
    dp[S][v]=total
  H[S]=sum(dp[S])
 full=size-1;V=e=0
 for S in range(1,size,2):
  value=dp[S][0]*H[full^S];V+=value;e+=(-1 if S.bit_count()%2==0 else 1)*value
 return {'H':H[full],'h':H[full^1],'V':V,'e':e}
T41=(56,13,25,48,34,6);T43=(52,13,24,49,34,6);T45=(52,13,24,33,42,6)
R=(386,508,489,497,485,193,257,320,32)
bank=(join(apex('100'),T41),join(apex('100'),T43),R)
profiles=[profile(T) for T in bank]
need([p['H'] for p in profiles]==[269,271,273],'base counts')
need(all((p['h'],p['V'],p['e'])==(15,441,15) for p in profiles),'rearmed coordinates')
counts=[];receivers=[]
for word in ('111000','111001'):
 for T in bank:
  p=profile(join(apex(word),T));counts.append(p['H']);receivers.append(p['V'])
need(counts==[24875,24877,24879,24881,24883,24885],'six counts')
need(receivers==[28119]*3+[27639]*3,'paid receiver defect')
print(json.dumps({'status':'PASS','method':'backward incoming-edge subset recurrence; no imports from candidate verifiers','base_profiles':profiles,'six_counts':counts,'six_receivers':receivers},indent=2))
