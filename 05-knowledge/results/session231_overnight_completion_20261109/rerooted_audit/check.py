"""Independent literal interfaces, exact algebra, and physical full-subset DP.
Build first: g++ -O3 -std=c++17 physical.cpp -o physical
No producer imports. Explicit gates remain effective under python -O.
"""
import itertools, json, subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
C=[4,1,18,7,75,31,175,63]; X='101000111000';Y='111011100110'
def need(x,m):
 if not x: raise RuntimeError(m)
def compose(a,b):
 if a is None:return b
 if b is None:return a
 h,q,A,G,T,S,E=a;v,w,B,F,U,V,W=b
 return h*v,q*w,A*w+h*B,F*q+v*G,A*F+T*v+h*U+E*V,S*v,h*W
def word(w):
 s=None
 for bit in w: s=compose(s,(1,2,int(bit),int(bit),0,int(bit),int(bit)))
 return s
def literal_state(rows,cut):
 n=len(rows);paths={0:[()]}
 for size in range(1,n+1):
  for ss in itertools.combinations(range(n),size):
   mask=sum(1<<i for i in ss)
   paths[mask]=[p for p in itertools.permutations(ss) if all(rows[u]>>v&1 for u,v in zip(p,p[1:]))]
 full=(1<<n)-1;h=len(paths[full]);q=A=G=T=0
 for m,p in paths.items():
  a=paths[full^m]; q+=len(p)*len(a)
  g=sum(bool(t) and bool(cut>>t[-1]&1) for t in p)
  b=sum(bool(t) and bool(cut>>t[0]&1) for t in a)
  A+=len(p)*b;G+=g*len(a);T+=g*b
 return h,q,A,G,T,sum(cut>>p[0]&1 for p in paths[full]),sum(cut>>p[-1]&1 for p in paths[full])
core=literal_state(C,124)
need(core==(45,1440,878,577,332,25,14),'literal core state')
need(word(Y)==(1,4096,3814,1655,1540,1,0),'Y state')
a=compose(compose(compose(word('0'),core),word(Y)),word('1'))
b=compose(compose(compose(word('0'),word(X)),core),word('1'))
need(tuple(v-u for u,v in zip(a,b))==(0,0,0,0,2,0,0),'all-order move certificate')
# These exact identities, h=1 contexts, X[-1]=0 and Y[0]=1 prove every move.
need(X[-1]=='0' and Y[0]=='1','context boundary bits')
records=[]
for k in range(31):
 R=4096**k; J=(R-1)//4095
 num=13617152*R+1602510*k+717259
 need(num%24843==0,'integrality')
 z=num//24843
 for j in range(k+1):
  s=compose(compose(word(X*j),core),word(Y*(k-j)))
  h,q,A,G,T,S,E=s
  need((h,q,A,G)==(45,1440*R,878*R+171630*J,577+2383200*J),'bare interface')
  need(A-T+h-E==z-2*j,'bare scalar')
  need(E==(14 if j==k else 0),'terminal endpoint')
 records.append({'k':k,'upper':z,'order':12*k+9})
def graph(k,j,full=False):
 bits=list(map(int,('10' if full else '')+X*j))+[(124>>i)&1 for i in range(8)]+list(map(int,Y*(k-j)+('10' if full else '')))
 # full suffix 10 consists of selected trailing vertex and unselected v.
 n=len(bits);off=(2 if full else 0)+12*j
 rows=[sum(1<<v for v in range(u+1,n)) for u in range(n)]
 for u in range(8):rows[off+u]=(rows[off+u]&~(255<<off))|(C[u]<<off)
 rows=[a|((1<<n) if not bits[u] else 0) for u,a in enumerate(rows)]+[sum(b<<u for u,b in enumerate(bits))]
 return rows

def strong(rows):
 for u in range(len(rows)):
  mask=1<<u
  while True:
   nxt=mask
   for v,a in enumerate(rows):
    if mask>>v&1:nxt|=a
   if nxt==mask:break
   mask=nxt
  if mask!=(1<<len(rows))-1:return False
 return True
for k in range(8):
 for j in range(k+1):need(strong(graph(k,j)),'physical strongness')
cases=[(0,0,False),(1,0,False),(1,1,False),(0,0,True),(1,0,True),(1,1,True)]
inputs=''.join(str(len(g))+'\n'+' '.join(map(str,g))+'\n' for g in [graph(*c) for c in cases])
p=subprocess.run([str(HERE/'physical')],input=inputs,text=True,capture_output=True,check=True)
answers=p.stdout.strip().split('END');physical=[]
for (k,j,full),answer in zip(cases,answers):
 states=[tuple(map(int,line.split())) for line in answer.strip().splitlines()]
 if not full:
  R=4096**k;J=(R-1)//4095;z=(13617152*R+1602510*k+717259)//24843
  expected=(z-2*j,31 if j==k else 45,878*R+171630*J+45,45)
  need(states[-1]==expected,'literal full graph bare state')
 elif k==0:
  need(states[-1]==(9909,45,15167,45),'literal original r')
  need(states[-2]==(9909,2303,9909,2303),'literal original v')
 physical.append({'case':(k,j,full),'apex':states[-1],'states':states})
f0,f1=physical[-2]['states'],physical[-1]['states']
indices0={'lead':0,'zero':1,**{'C'+str(i):2+i for i in range(8)},'trail':22,'v':23,'r':24}
indices1={'lead':0,'zero':1,**{'C'+str(i):14+i for i in range(8)},'trail':22,'v':23,'r':24}
diffs={}
for name,i in indices0.items():
 delta=tuple(y-x for x,y in zip(f0[i],f1[indices1[name]])); diffs[name]=delta
 D,E,U,d=delta; J=D-E+U; V=U-d
 if name in ('r','v'):need((J,V)==(-2,0),'closing gap')
 elif name=='trail':need((J,V)==(8,-2),'trailing parity obstruction')
 else:need(J < -2 and V < 0,'named port universal gadget inequality')
# Literal deletion controls, including the global sink that changes root E.
g=graph(0,0,True)
def induced(rows,keep):
 return [sum(((rows[u]>>v)&1)<<i for i,v in enumerate(keep)) for u in keep]
delgraphs=[induced(g,[u for u in range(13) if u!=11]),induced(g,[u for u in range(13) if u not in (10,11)])]
raw=''.join(str(len(g))+'\n'+' '.join(map(str,g))+'\n' for g in delgraphs)
res=subprocess.run([str(HERE/'physical')],input=raw,text=True,capture_output=True,check=True)
ds=[[tuple(map(int,l.split())) for l in block.strip().splitlines()] for block in res.stdout.strip().split('END') if block.strip()]
need(ds[0][-1]==(2303,0,7606,45),'literal deleted-v sidecar')
need(ds[1][-1]==(2303,31,3803,45),'literal remove-trailing-sink sidecar')
# Exact scalar recurrence proves the displayed rational baseline integral and odd.
need(13617152+717259==24843*577,'baseline zero')
need(4095*1602510==24843*264150,'baseline recurrence k coefficient')
need(1602510-4095*717259==-24843*118165,'baseline recurrence constant')
# Terminal endpoint defect: reordered seven-value exception, duplicates thereafter.
for k in range(1,15):
 offsets=[-2*j for j in range(k)]+[-2*k+14]
 unit=len(set(offsets))==k+1 and max(offsets)-min(offsets)==2*k
 need(unit==(k==6),'unordered bare graft exception')
print(json.dumps({'status':'PASS','literal_deletions':[ds[0][-1],ds[1][-1]],'k1_named_port_differences':diffs,'all_order_certificate':{'core':core,'flanked_left':a,'flanked_right':b},'finite_algebra_through_k':30,'finite_strongness_through_k':7,'intervals':records[:4],'physical_subset_controls':physical},indent=2))
