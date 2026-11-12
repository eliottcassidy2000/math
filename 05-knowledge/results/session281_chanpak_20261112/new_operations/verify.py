#!/usr/bin/env python3
"""Restored exact bounded central-trade search. Guards survive python -O."""
import collections,itertools,json,math

def require(q,msg):
 if not q:raise ValueError(msg)

def endpoints(rows):
 n=len(rows);a=[[0]*n for _ in range(1<<n)]
 for j in range(n):a[1<<j][j]=1
 for s in range(1,1<<n):
  for j in range(n):
   if s>>j&1:
    r=s^(1<<j)
    if r:a[s][j]=sum(a[r][i] for i in range(n) if rows[i]>>j&1)
 return a

def converse(rows):
 full=(1<<len(rows))-1
 return [full^(1<<i)^v for i,v in enumerate(rows)]

def representatives(n):
 pairs=list(itertools.combinations(range(n),2));reps={};class_sizes=collections.Counter()
 for bits in range(1<<len(pairs)):
  rows=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
  groups=[[i for i in range(n) if rows[i].bit_count()==d] for d in range(n)]
  key=min(sum(((rows[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(pairs)) for g in itertools.product(*(itertools.permutations(b) for b in groups)) for p in [sum(g,())])
  if key not in reps:reps[key]=rows
  class_sizes[key]+=1
 require(sum(class_sizes.values())==1<<len(pairs),'labelled universe')
 return reps,class_sizes

def cut_profiles(rows):
 n=len(rows);full=(1<<n)-1;E=endpoints(rows);B=endpoints(converse(rows))
 H=[sum(v) for v in E];H[0]=1;h=H[-1]
 Q=sum(H[s]*H[full^s] for s in range(1<<n))
 left=[sum(B[s][j]*H[full^s] for s in range(1<<n)) for j in range(n)]
 right=[sum(E[s][j]*H[full^s] for s in range(1<<n)) for j in range(n)]
 dl=[sum((-1)**s.bit_count()*B[s][j]*H[full^s] for s in range(1<<n)) for j in range(n)]
 dr=[sum((-1)**s.bit_count()*E[s][j]*H[full^s] for s in range(1<<n)) for j in range(n)]
 for out in range(1<<n):
  inc=full^out
  U=h+sum(left[j] for j in range(n) if out>>j&1);V=h+sum(right[j] for j in range(n) if inc>>j&1)
  d=h+sum(dl[j] for j in range(n) if out>>j&1);e=h+sum(dr[j] for j in range(n) if inc>>j&1)
  count=0
  for s in range(1<<n):
   p=1 if s==0 else sum(E[s][i] for i in range(n) if inc>>i&1)
   r=full^s;q=1 if r==0 else sum(B[r][j] for j in range(n) if out>>j&1)
   count+=p*q
  lifted=[rows[i]|((1<<n) if inc>>i&1 else 0) for i in range(n)]+[out]
  yield (h,Q,U,V,d,e),count,lifted

def direct_profile(rows,root):
 E=endpoints(rows);B=endpoints(converse(rows));full=(1<<len(rows))-1
 H=[sum(x) for x in E];H[0]=1;deletion=full^(1<<root)
 h=H[deletion];Q=sum(H[s]*H[deletion^s] for s in range(1<<len(rows)) if s&~deletion==0)
 U=sum(B[s][root]*H[full^s] for s in range(1<<len(rows)))
 V=sum(E[s][root]*H[full^s] for s in range(1<<len(rows)))
 d=sum((-1)**(s.bit_count()+1)*B[s][root]*H[full^s] for s in range(1<<len(rows)))
 e=sum((-1)**(s.bit_count()+1)*E[s][root]*H[full^s] for s in range(1<<len(rows)))
 return (h,Q,U,V,d,e),H[-1]

def fuse(A,r,B,s):
 aa=[i for i in range(len(A)) if i!=r];bb=[i for i in range(len(B)) if i!=s]
 rows=[0]*(1+len(aa)+len(bb));ma={r:0,**{v:i+1 for i,v in enumerate(aa)}};mb={s:0,**{v:len(aa)+i+1 for i,v in enumerate(bb)}}
 for T,m in [(A,ma),(B,mb)]:
  for u in range(len(T)):
   for v in range(len(T)):
    if T[u]>>v&1:rows[m[u]]|=1<<m[v]
 for u in aa:
  for v in bb:rows[ma[u]]|=1<<mb[v]
 return rows

def main():
 reps,counts=representatives(6);require(len(reps)==56,'order-six class count')
 fibers=collections.defaultdict(dict);tested=0;coprime=0
 for rows in reps.values():
  for key,H,lifted in cut_profiles(rows):
   tested+=1
   if math.gcd(key[0],15)!=1:continue
   coprime+=1;fibers[key].setdefault(H,lifted)
 output=[]
 for key,bank in sorted(fibers.items()):
  if len(bank)>1:
   hs=sorted(bank);require(all(b-a!=2 for a,b in itertools.combinations(hs,2)),'no coprime +2 trade')
   witnesses=[]
   for h in hs:
    rows=bank[h];require(direct_profile(rows,6)==(key,h),'whole-tournament endpoint profile')
    witnesses.append({'H':h,'rows':rows,'root':6})
   output.append({'carrier_without_H':key,'witnesses':witnesses})
 require([(tuple(x['carrier_without_H']),[v['H'] for v in x['witnesses']]) for x in output]==[((17,224,89,152,11,6),[65,69]),((17,224,152,89,6,11),[65,69]),((29,320,122,227,20,9),[91,95]),((29,320,227,122,9,20),[91,95])],'collision groups')
 contexts=[([0],0),([2,0],0),([2,0],1),([2,4,1],0)];controls=0
 for group in output:
  a,b=group['witnesses'];gap=b['H']-a['H'];require(gap==4,'gap')
  for L,r in contexts:
   for R,s in contexts:
    out=[]
    for w in (a,b):
     G=fuse(fuse(L,r,w['rows'],6),0,R,s);out.append(sum(endpoints(G)[-1]))
    require(out[1]-out[0]==gap*direct_profile(L,r)[0][0]*direct_profile(R,s)[0][0],'sandwich propagation');controls+=1
 require(sorted(58*i+60*j for i in range(6) for j in range(2))==[0,58,60,116,118,174,176,232,234,290,292,350],'mixed bank offsets')
 print(json.dumps({'status':'PASS','scope':'Every seven-vertex pointed tournament with deletion count coprime to15; isomorphic deletion representatives and every cut, no global omission claim','labelled_deletion_universe':32768,'deletion_isomorphism_classes':len(reps),'pointed_presentations':tested,'coprime_presentations':coprime,'distinct_coprime_signatures':len(fibers),'nonconstant_fibers':output,'literal_context_controls':controls},sort_keys=True,indent=2))
if __name__=='__main__':main()
