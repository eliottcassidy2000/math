#!/usr/bin/env python3
"""Independent direct-subset path DP for selected fiber claims."""
import itertools as it,json
import sympy as s

def totals(rows,start=None):
 n=len(rows);N=1<<n; d=[[0]*n for _ in range(N)];out=[0]*N;out[0]=1
 for i in range(n):
  if start is None or start==i:d[1<<i][i]=1
 for mask in range(1,N):
  for v,c in enumerate(d[mask]):
   if not c:continue
   todo=rows[v]&~mask
   while todo:
    bit=todo&-todo;todo-=bit;w=bit.bit_length()-1;d[mask|bit][w]+=c
  out[mask]=sum(d[mask])
 return out

def transpose(r):return [sum(1<<j for j in range(len(r)) if r[j]>>i&1) for i in range(len(r))]
def packet(r,x):
 n=len(r);full=(1<<n)-1;sub=totals(r);t=[row|((1<<n) if not x>>i&1 else 0) for i,row in enumerate(r)]+[x]
 st=totals(t,n);en=totals(transpose(t),n);allp=totals(t)
 vals=[0]*4
 for S in range(1<<n):
  weight=sub[full^S];sgn=(-1)**S.bit_count();mask=S|(1<<n)
  vals[0]+=weight*st[mask];vals[1]+=weight*en[mask]
  vals[2]+=weight*sgn*st[mask];vals[3]+=weight*sgn*en[mask]
 return (allp[-1],*vals)

def circ(n,steps):return [sum(1<<((i+d)%n) for d in steps) for i in range(n)]
def spectrum(r,k):
 ps=[packet(r,sum(1<<i for i in c)) for c in it.combinations(range(len(r)),k)]
 if len({p[1:] for p in ps})!=1:raise RuntimeError('responses differ')
 return sorted({p[0] for p in ps})
r=[14,12,8,0];p0=packet(r,0)
M=s.Matrix([[a-b for a,b in zip(packet(r,x),p0)] for x in [1,2,3,4,8]])
if M.det()!=24:raise RuntimeError('rank-five witness')
checks={'rank_five_determinant':int(M.det()),'examples':[]}
for n,steps,k,want in [(3,[1],1,[5]),(5,[1,2],2,[41,43]),(7,[1,2,4],3,[657]),(7,[1,2,3],2,[559,561,575]),(7,[1,2,3],3,[629,641,643,659])]:
 got=spectrum(circ(n,steps),k)
 if got!=want:raise RuntimeError((n,steps,k,got,want))
 checks['examples'].append(dict(n=n,steps=steps,k=k,counts=got))
# Independent fixed-deletion switch nonexistence in labelled cores through4.
small=[]
for n in range(1,5):
 pairs=list(it.combinations(range(n),2));cores=0
 for bits in range(1<<len(pairs)):
  r=[0]*n
  for j,(a,b) in enumerate(pairs):
   if bits>>j&1:r[a]|=1<<b
   else:r[b]|=1<<a
  seen={}
  for x in range(1<<n):
   p=packet(r,x)
   if p[1:] in seen and seen[p[1:]]!=p[0]:raise RuntimeError('early central switch')
   seen[p[1:]]=p[0]
  cores+=1
 small.append([n,cores])
checks['no_central_switch_labelled_core_universe']=small
r=[226,452,393,275,39,78,156,312,113]
pa,pb=packet(r,21),packet(r,11)
if pa!=(14207,18906,34455,2238,1119) or pb!=(14209,18906,34455,2238,1119):raise RuntimeError(('n9 pair',pa,pb))
hs=totals(r);full=(1<<9)-1;Q=sum(hs[S]*hs[full^S] for S in range(1<<9))
if hs[-1]!=3357 or Q!=50004:raise RuntimeError('n9 fixed core')
checks['n9_central_pair']={'rows':r,'cuts':[21,11],'H_U_V_d_e':[pa,pb],'h':hs[-1],'Q':Q}
r=circ(11,[1,3,4,5,9]);hs=totals(r);full=(1<<11)-1;Q=sum(hs[S]*hs[full^S] for S in range(1<<11));p=packet(r,3)
if (hs[-1],Q,p[0])!=(95095,1776500,354673):raise RuntimeError('Paley11 control')
checks['Paley11_h_Q_H2']=[hs[-1],Q,p[0]]
print(json.dumps({'status':'PASS','checks':checks},sort_keys=True,indent=2))
