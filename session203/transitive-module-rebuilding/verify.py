#!/usr/bin/env python3
"""Exact transitive-module path/cycle skeletons and physical C1-to-C3 replacement."""
import itertools,json

def need(b,s):
 if not b:raise RuntimeError(s)

FAMILIES=[
 {'name':'A','P':[2,24,3,37,13,23],'T':[2,8,67,101,15,87,19],
  'flip':[1,4],'mask':15,'extra':1,
  'hpP':[25,59,41,8],'hpT':[51,119,82,16],
  'cycP':[2,3,1],'cycT':[3,5,2]},
 {'name':'B','P':[32,5,1,39,47,70,31],'T':[0,5,129,39,47,199,159,27],
  'flip':[0,5],'mask':5,'extra':0,
  'hpP':[25,45,20],'hpT':[51,91,40],
  'cycP':[1,1],'cycT':[0], 'cycReduced':[3,3]},
 {'name':'C','P':[2,4,41,3,15,91,31],'T':[0,5,41,131,15,219,159,23],
  'flip':[0,1],'mask':27,'extra':1,
  'hpP':[25,53,27],'hpT':[51,107,54],
  'cycP':[1,1],'cycT':[0], 'cycReduced':[3,3]}
]

def valid(rows):
 n=len(rows)
 for i,a in enumerate(rows):
  need(0<=a<1<<n and not(a>>i&1),'bad row')
  for j in range(i):need(((a>>j)&1)+((rows[j]>>i)&1)==1,'bad pair')

def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a

def skeleton(rows,mask,cyclic=False):
 n=len(rows);ans=[0]*(n+2)
 orders=((0,)+p for p in itertools.permutations(range(1,n))) if cyclic else itertools.permutations(range(n))
 for p in orders:
  a=[1]
  gaps=[(p[i],p[(i+1)%n]) for i in range(n)] if cyclic else [(None,p[0])]+list(zip(p,p[1:]))+[(p[-1],None)]
  for left,right in gaps:
   empty=left is None or right is None or bool(rows[left]>>right&1)
   filled=(left is None or not(mask>>left&1)) and (right is None or bool(mask>>right&1))
   b=[0]*(len(a)+1)
   for j,x in enumerate(a):
    if empty:b[j]+=x
    if filled:b[j+1]+=x
   a=b
  for j,x in enumerate(a):ans[j]+=x
 return trim(ans)

def inflate(rows,mask,k):
 n=len(rows);a=rows.copy()+[0]*k
 for j in range(n,n+k):
  for i in range(n):
   if mask>>i&1:a[j]|=1<<i
   else:a[i]|=1<<j
  for i in range(n,j):a[i]|=1<<j
 valid(a);return a

def delete(rows,v):
 ix=[i for i in range(len(rows)) if i!=v]
 return [sum(1<<j for j,w in enumerate(ix) if rows[u]>>w&1) for u in ix]

def literal_dp(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)];c=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 c[1][0]=1
 for s in range(1,1<<n):
  for v in range(n):
   q=rows[v]&~s
   while q:
    b=q&-q;q-=b;w=b.bit_length()-1
    d[s|b][w]+=d[s][v];c[s|b][w]+=c[s][v]
 return sum(d[-1]),sum(c[-1][v] for v in range(n) if rows[v]&1)

def evaluate(coeff,k):
 # r! S(k,r) counts surjections from k ordered module vertices to r slots.
 s=[0]*(len(coeff));s[0]=1
 for _ in range(k):s=[0]+[r*(s[r]+s[r-1]) for r in range(1,len(s))]
 return sum(a*b for a,b in zip(coeff,s))

def main():
 out=[];graphs=0
 for f in FAMILIES:
  p,t=f['P'],f['T'];n=len(p);valid(p);valid(t);need(len(t)==n+1,'wrong birth')
  i,j=f['flip'];need(p[i]>>j&1,'named arc absent')
  q=p.copy();q[i]^=1<<j;q[j]^=1<<i
  need([a&((1<<n)-1) for a in t[:-1]]==q,'not exactly named reversal and birth')
  tm=f['mask']|(f['extra']<<n)
  for key,rows,mask,cyc in [('hpP',p,f['mask'],False),('hpT',t,tm,False),('cycP',p,f['mask'],True),('cycT',t,tm,True)]:
   need(skeleton(rows,mask,cyc)==f[key],f['name']+' wrong '+key)
  if f['name'] in ['B','C']:
   need(t[0]==0 and tm&1,'vertex0 not universal sink for all module sizes')
   red=delete(t,0);rm=tm>>1;need(skeleton(red,rm,True)==f['cycReduced'],'wrong reduced cycles')
   need(skeleton(red,rm)==f['hpT'],'deleting sink changes path polynomial')
   # Deleting old0 first leaves all other old edges unchanged; new last core vertex replaces it.
   need(delete(p,0)==delete(t[:-1],0),'replacement changes retained old edge')
  samples=[]
  for k in range(5):
   pp=inflate(p,f['mask'],k);tt=inflate(t,tm,k)
   ph,pc=literal_dp(pp);th,tc=literal_dp(tt);graphs+=2
   need((ph,pc)==(evaluate(f['hpP'],k),evaluate(f['cycP'],k)),'parent DP mismatch')
   need((th,tc)==(evaluate(f['hpT'],k),evaluate(f['cycT'],k)),'child DP mismatch')
   need(th==2*ph+1,'doubling mismatch')
   sample={'k':k,'parent_H':ph,'parent_C':pc,'child_H':th,'child_C':tc}
   if f['name'] in ['B','C']:
    rr=delete(tt,0);rh,rc=literal_dp(rr);graphs+=1
    need((rh,rc)==(th,3),'reduced child not C3 companion')
    sample['replacement_H']=rh;sample['replacement_C']=rc
   samples.append(sample)
  out.append({**f,'samples':samples})
 print(json.dumps({'status':'PASS','literal_graphs':graphs,'module_sizes_checked':list(range(5)),'families':out},indent=2,sort_keys=True))
if __name__=='__main__':main()
