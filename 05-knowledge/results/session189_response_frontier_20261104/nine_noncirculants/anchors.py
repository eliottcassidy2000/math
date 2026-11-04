from pathlib import Path
import json,subprocess,importlib.util
D=Path(__file__).resolve().parent
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def mods(Q):
 n=len(Q);out=[]
 for s in range(1,(1<<n)-1):
  if s.bit_count()<2:continue
  if all(not(s>>i&1) and (Q[i]&s) in (0,s) or s>>i&1 for i in range(n)):out.append(s)
 return out
def endpoints(Q):
 n=len(Q);joint=[[0]*n for _ in range(n)]
 for start in range(n):
  dp=[[0]*n for _ in range(1<<n)];dp[1<<start][start]=1
  for m in range(1,1<<n):
   for last in range(n):
    value=dp[m][last]
    if not value:continue
    todo=Q[last]&~m
    while todo:
     bit=todo&-todo;todo-=bit;j=bit.bit_length()-1;dp[m|bit][j]+=value
  joint[start]=dp[-1]
 return joint


rows=[[192,289,99,71,15,153,306,94,189],[320,193,83,71,43,15,432,61,190],[448,209,195,71,45,15,432,312,62],[384,177,147,135,105,77,15,368,126]]
subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'anchors.cpp'),'-o',str(D/'anchors')],check=True)
raw=subprocess.check_output([str(D/'anchors')],input=('4\n'+''.join(' '.join(map(str,Q))+'\n' for Q in rows)).encode()).decode();out=[]
(D/'anchor_values.txt').write_text(raw)
for Q,line in zip(rows,raw.splitlines()):
 J=endpoints(Q);S=list(map(sum,J));E=[sum(J[i][j] for i in range(9)) for j in range(9)];p=[(x+y)%2 for x,y in zip(S,E)];vals=list(map(int,line.split()));start=[vals[9*i:9*i+9] for i in range(9)];d0=vals[81:90];d2=[vals[90+9*i:99+9*i] for i in range(9)];patterns=[]
 for V in range(512):
  d=[(d0[i]+2*sum(d2[i][j] for j in range(9) if V>>j&1))%4 for i in range(9)];t=[(S[i]-E[i]+2*sum(start[i][j] for j in range(9) if V>>j&1))%4 for i in range(9)];a=[(2*t[i]+2*d[i])%8 for i in range(9)];patterns.append(dict(V=V,first=d,second=a))
 out.append(dict(rows=Q,joint_endpoints=J,nontrivial_modules=mods(Q),S=S,E=E,p=p,sigma=[(S[i]-E[i]+2*sum(start[i]))%4 for i in range(9)],start=start,d0=d0,d2=d2,patterns=patterns,blind=[r['V'] for r in patterns if not any(r['first']) and not any(r['second'])]))
(D/'anchor_data.json').write_text(json.dumps(out,separators=(',',':'))+'\n');print(json.dumps([{k:v for k,v in q.items() if k not in ['patterns','start','d0','d2']} for q in out],indent=2))
