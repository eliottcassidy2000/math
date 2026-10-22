import itertools,json,collections,math
from pathlib import Path

def step(s,b,M=None):
 h,a,p=s
 t=(h,2*a+1,1) if b else (h+a+p,2*a,0)
 return t if M is None else (t[0]%M,t[1]%M,t[2])
def state(w,M=None):
 s=(1 if M is None else 1%M,0,0)
 for b in w:s=step(s,int(b),M)
 return s
def h(w):return state(w)[0]
def dp(w):
 n=len(w)+1; out=[0]*n
 for i in range(n-1):
  for j in range(i+1,n-1):out[i]|=1<<j
  if w[i]=='1':out[-1]|=1<<i
  else:out[i]|=1<<(n-1)
 d={(1<<i,i):1 for i in range(n)}
 for k in range(1,n):
  e=collections.defaultdict(int)
  for (mask,v),c in d.items():
   avail=out[v]&~mask
   while avail:
    q=avail&-avail;avail-=q;e[mask|q,q.bit_length()-1]+=c
  d=e
 return sum(d.values())
def check(x,msg):
 if not x:raise RuntimeError(msg)
checks=0
for n in range(9):
 for t in itertools.product('01',repeat=n):
  w=''.join(t);check(h(w)==dp(w),w);checks+=1
rows=[]
for M in range(1,81):
 start=state('',M);seen={start};q=collections.deque([start])
 while q:
  s=q.popleft()
  for b in [0,1]:
   t=step(s,b,M)
   if t not in seen:seen.add(t);q.append(t)
 expected=2*M*M if M%2 else M*M//2
 check(len(seen)==expected,f'reach {M}')
 # independent partition refinement to Moore equivalence
 ss=sorted(seen);idx={s:i for i,s in enumerate(ss)}; part=[s[0] for s in ss]
 while True:
  keys=[(part[i],part[idx[step(s,0,M)]],part[idx[step(s,1,M)]]) for i,s in enumerate(ss)]
  dic={k:i for i,k in enumerate(sorted(set(keys)))};new=[dic[k] for k in keys]
  if len(set(new))==len(set(part)):break
  part=new
 minimal=len(set(new));expected_min=1 if M<=2 else expected
 check(minimal==expected_min,f'min {M}')
 probes={(s[0],(s[0]+s[1]+s[2])%M,(s[0]+2*s[1]+2)%M) for s in ss}
 check(len(probes)==expected_min,f'probes {M}')
 rows.append([M,len(seen),minimal]);checks+=3
# Hankel minor and exact determinant
U=['','1','10','11']; V=['','0','10','00']; mat=[[h(u+v) for v in V] for u in U]
def det(a):
 if not a:return 1
 return sum((-1)**j*a[0][j]*det([r[:j]+r[j+1:] for r in a[1:]]) for j in range(len(a)))
print(json.dumps({'word_dp_controls':511,'modular_states':[{'M':r[0],'reachable':r[1],'minimal':r[2]} for r in rows], 'hankel_rows':U,'hankel_columns':V,'hankel':mat,'determinant':det(mat),'checks':checks},indent=2))
