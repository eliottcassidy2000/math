import itertools, functools,json
C=[4,1,18,7,75,31,175,63];X='101000111000';Y='111011100110'
def family(k,j):
 blocks=[];bits=[];names=[]
 for name,word in [('lead','1'),('zero','0'),('X',X*j)]:
  for i,c in enumerate(word): blocks.append([len(bits)]);bits.append(int(c));names.append(name if len(word)==1 else name+str(i))
 off=len(bits); blocks.append(list(range(off,off+8)));bits.extend([(124>>i)&1 for i in range(8)]);names.extend('C'+str(i) for i in range(8))
 for name,word in [('Y',Y*(k-j)),('trail','1')]:
  for i,c in enumerate(word):blocks.append([len(bits)]);bits.append(int(c));names.append(name if len(word)==1 else name+str(i))
 n=len(bits);rows=[sum(1<<v for v in range(u+1,n)) for u in range(n)]
 for u in range(8): rows[off+u]=(rows[off+u]&~(((1<<8)-1)<<off))|(C[u]<<off)
 v=n;r=n+1
 rows=[x|(1<<v)|(0 if bits[u] else 1<<r) for u,x in enumerate(rows)]+[1<<r,sum(b<<u for u,b in enumerate(bits))]
 blocks.append([v]);names+=['v','r']
 return rows,blocks,names

def count(rows,blocks,gaps):
 # A gap (left special, right special), None indicates free endpoint.
 b=len(gaps);dp={(0,)*b:1}
 for block in blocks:
  n=len(block)
  if not n:continue
  # endpoint matrix for paths on each local subset
  ends={}
  for s in range(1,1<<n):
   mat={}
   for first in range(n):
    if not s>>first&1:continue
    if s==1<<first:mat[first,first]=1;continue
    for last in range(n):
     if last==first or not s>>last&1:continue
     t=s^(1<<last)
     z=sum(ends[t].get((first,p),0) for p in range(n) if rows[block[p]]>>block[last]&1)
     if z:mat[first,last]=z
   ends[s]=mat
  # transitions given incoming state and subset assigned to a gap
  tran={}
  for g,(left,right) in enumerate(gaps):
   for status in range(3):
    tran[g,status,0]=[(status,1)]
    for s,mat in ends.items():
     vals=[0,0,0]
     for (first,last),z in mat.items():
      if status==0 and left is not None and not rows[left]>>block[first]&1:continue
      good=right is None or bool(rows[block[last]]>>right&1)
      vals[2 if good else 1]+=z
     tran[g,status,s]=[(i,z) for i,z in enumerate(vals) if z]
  assignments=[]
  for labels in itertools.product(range(b),repeat=n):
   masks=[0]*b
   for i,g in enumerate(labels):masks[g]|=1<<i
   assignments.append(masks)
  nd={}
  for st,w in dp.items():
   for masks in assignments:
    options=[tran[g,st[g],masks[g]] for g in range(b)]
    for comb in itertools.product(*options):
     dest=tuple(z[0] for z in comb);weight=w
     for z in comb:weight*=z[1]
     nd[dest]=nd.get(dest,0)+weight
  dp=nd
 total=0
 for st,w in dp.items():
  ok=True
  for s,(left,right) in zip(st,gaps):
   if s==1 or (s==0 and left is not None and right is not None and not rows[left]>>right&1):ok=False;break
  if ok:total+=w
 return total

def states(k,j):
 rows,blocks,names=family(k,j);r=len(rows)-1
 D=count(rows,blocks,[(None,r),(r,None)])
 result={}
 for a,name in enumerate(names):
  if name not in ['lead','zero','trail','v','r'] and not name.startswith('C'):continue
  bs=[[v for v in B if v!=a] for B in blocks]
  if a==r:
   d=count(rows,blocks,[(None,None)])
   E=count(rows,blocks,[(None,r)])
   U=count(rows,blocks,[(r,None),(None,None)])
  else:
   d=count(rows,bs,[(None,r),(r,None)])
   E=count(rows,bs,[(None,r),(r,a)])
   U=count(rows,bs,[(a,r),(r,None),(None,None)])+count(rows,bs,[(a,None),(None,r),(r,None)])
  result[name]=(D,E,U,d)
 return result
if __name__=='__main__':
 import sys
 k=int(sys.argv[1]) if len(sys.argv)>1 else 0
 out={j:states(k,j) for j in range(k+1)}
 print(json.dumps(out,indent=2))
