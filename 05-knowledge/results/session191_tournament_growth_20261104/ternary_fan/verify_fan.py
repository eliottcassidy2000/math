#!/usr/bin/env python3
"""Exact self-contained checks for the two-port ternary fan. No external input."""
import hashlib, itertools, json
from pathlib import Path

CORES={
 'P':[10,4,9,2],
 'Q':[26,20,25,2,8],
 'A':[58,52,57,34,8,16],
 'C':[58,52,57,2,8,24],
 'B':[26,52,57,2,8,25],
}
EXPECTED={'P':(5,7),'Q':(9,11),'A':(17,21),'C':(15,19),'B':(19,21)}
R,V=1,3

def require(b,msg):
 if not b: raise ValueError(msg)

def valid(rows):
 n=len(rows)
 require(all(not (rows[i]>>i&1) and rows[i]>>n==0 for i in range(n)), 'loop or stray bit')
 require(all(((rows[i]>>j&1)+(rows[j]>>i&1))==1 for i in range(n) for j in range(i)), 'not tournament')

def strong(rows):
 n=len(rows);allmask=(1<<n)-1
 for v in range(n):
  seen=1<<v
  while True:
   nxt=seen
   for u in range(n):
    if seen>>u&1:nxt|=rows[u]
   if nxt==seen:break
   seen=nxt
  if seen!=allmask:return False
 return True

def profile(rows,root):
 valid(rows);n=len(rows);N=1<<n;end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for u in range(n):end[1<<u][u]=start[1<<u][u]=1
 for mask in range(1,N):
  for u in range(n):
   if not(mask>>u&1):continue
   old=mask^(1<<u)
   if old:
    end[mask][u]=sum(end[old][z] for z in range(n) if rows[z]>>u&1)
    start[mask][u]=sum(start[old][z] for z in range(n) if rows[u]>>z&1)
  h[mask]=sum(end[mask])
 F=sum(end[mask][root]*h[(N-1)^mask] for mask in range(1,N))
 J=[]
 for v in range(n):
  d=[[0]*n for _ in range(N)];d[1<<v][v]=1
  for mask in range(1,N):
   for u in range(n):
    if not(mask>>u&1):continue
    old=mask^(1<<u)
    if old:d[mask][u]=sum(d[old][z] for z in range(n) if rows[z]>>u&1)
  J.append(sum(d[mask][root]*h[(N-1)^mask] for mask in range(1,N)))
 return {'H':h[-1],'F':F,'s':start[-1],'e':end[-1],'J':J,'strong':strong(rows)}

def brute(rows,root):
 n=len(rows);p={};ends={};j=[0]*n
 for mask in range(1<<n):
  names=[u for u in range(n) if mask>>u&1];p[mask]=[]
  for word in itertools.permutations(names):
   if all(rows[u]>>v&1 for u,v in zip(word,word[1:])):p[mask].append(word)
 full=(1<<n)-1;s=[0]*n;e=[0]*n;F=0
 for word in p[full]:s[word[0]]+=1;e[word[-1]]+=1
 for mask,words in p.items():
  for word in words:
   if word and word[-1]==root:
    z=len(p[full^mask]);F+=z;j[word[0]]+=z
 return {'H':len(p[full]),'F':F,'s':s,'e':e,'J':j,'strong':strong(rows)}

def wrapper(core,k):
 require(k>=1,'positive clock');n=len(core);w=n+k;rows=core+[0]*(k+1)
 # Core->w only at R; every z beats every core name; w beats every z.
 rows[R]|=1<<w
 rows[w]=sum(1<<u for u in range(n+k) if u!=R)
 for i in range(k):rows[n+i]=((1<<n)-1)|sum(1<<(n+j) for j in range(i+1,k))
 return rows,w

def path_only(rows):
 n=len(rows);N=1<<n;end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)]
 for u in range(n):end[1<<u][u]=start[1<<u][u]=1
 for mask in range(1,N):
  for u in range(n):
   if not(mask>>u&1):continue
   old=mask^(1<<u)
   if old:
    end[mask][u]=sum(end[old][z] for z in range(n) if rows[z]>>u&1)
    start[mask][u]=sum(start[old][z] for z in range(n) if rows[u]>>z&1)
 return sum(end[-1]),start[-1]

def main():
 receipt={'cores':{},'wrappers':[]}
 for name,rows in CORES.items():
  p=profile(rows,R);require(p==brute(rows,R),'independent core mismatch '+name)
  require((p['H'],p['F'])==EXPECTED[name] and p['strong'] and p['J'][V]==1,'core profile '+name)
  receipt['cores'][name]={'rows':rows,**p}
 require(CORES['Q'][:4]==[r|(16 if u!=V else 0) for u,r in enumerate(CORES['P'])], 'first birth')
 require(CORES['Q'][4]==1<<V, 'first birth exit')
 require(CORES['A'][:5]==[r|(32 if u!=4 else 0) for u,r in enumerate(CORES['Q'])], 'second birth')
 require(CORES['A'][5]==16,'second birth exit')
 for a,b,edge in [('A','C',(3,5)),('C','B',(0,5))]:
  changed=[(i,j) for i in range(6) for j in range(i+1,6) if ((CORES[a][i]^CORES[b][i])>>j)&1]
  require(changed==[edge],'specified physical flip')
 for k in range(1,7):
  x=1<<k;N=5+7*x
  for name,core in CORES.items():
   rows,w=wrapper(core,k);h,s=path_only(rows);p=receipt['cores'][name]
   require(strong(rows),'lift strong')
   require(h==p['H']+x*p['F'],'lift count')
   require(s[:len(core)]==p['J'],'old start vector')
   require(s[len(core):w]==[(1<<(k-i-1))*p['F'] for i in range(k)],'clock start vector')
   require(s[w]==p['H'] and s[V]==1,'outer and retained starts')
   if name=='A':require(h==3*N+2,'lower ternary branch')
   if name=='B':require(h==3*N+4,'upper ternary branch')
   if name=='C':require(h<3*N+2,'nonovershoot')
   if k<=2:
    lifted=profile(rows,w)
    require(lifted['F']==h,'outer zero response')
   receipt['wrappers'].append({'name':name,'k':k,'rows':rows,'H':h,'start':s})
 out=Path(__file__).with_name('fan_receipt.json')
 raw=(json.dumps(receipt,indent=2,sort_keys=True)+'\n').encode();out.write_bytes(raw)
 print(json.dumps({'status':'PASS','independent_core_profiles':5,'wrapped_checks':30,'clocks':[1,6],'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
