#!/usr/bin/env python3
"""Independent audit: forbidden triangles, literal permutations, motif orbit expansion."""
import argparse,array,hashlib,itertools,json,math
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent
count=Counter()
def need(x,msg):
 count[msg]+=1
 if not x:raise ValueError(msg)
def vertices(x):
 while x:
  b=x&-x;x-=b;yield b.bit_length()-1
def tournament(n,S):return [sum(1<<((v+d)%n) for d in S) for v in range(n)]
def induced(A,V):return [sum(((A[v]>>w)&1)<<j for j,w in enumerate(V)) for v in V]
def transitive_data(A):
 n=len(A);cycpair=[[] for _ in A]
 for i,j,k in itertools.combinations(range(n),3):
  if ((A[i]>>j)&1)==((A[j]>>k)&1)==((A[k]>>i)&1):
   for v,rest in [(i,(j,k)),(j,(i,k)),(k,(i,j))]:cycpair[v].append(sum(1<<u for u in rest))
 ok=bytearray(1<<n);ok[0]=1;F=Counter({0:1});roots=[Counter({0:1}) for _ in A];rootcontained=[Counter() for _ in A]
 for S in range(1,1<<n):
  bit=S&-S;v=bit.bit_length()-1;rest=S^bit
  if not ok[rest]:continue
  if any(rest&e==e for e in cycpair[v]):continue
  ok[S]=1;k=S.bit_count();F[k]+=1
  for u,row in enumerate(A):
   if S&row==S:roots[u][k]+=1
   if S>>u&1:rootcontained[u][k]+=1
 def tup(c):return tuple(c[i] for i in range(max(c)+1))
 return tup(F),[tup(c) for c in roots],[tup(c) for c in rootcontained]
def multiply(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 while len(c)>1 and c[-1]==0:c.pop()
 return tuple(c)
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,y in enumerate(b):c[i]+=y
 while len(c)>1 and c[-1]==0:c.pop()
 return tuple(c)
def shifted(a,k=1):return (0,)*k+tuple(a)
def bp(r):return tuple(math.comb(r,i) for i in range(r+1))
def lift(A,r):
 n=len(A);out=[row|(1<<n) for row in A]+[sum(1<<(n+1+i) for i in range(r))]
 for i in range(r):out.append(((1<<n)-1)|sum(1<<(n+1+j) for j in range(i+1,r)))
 return out
def ham_counts(A):
 n=len(A);M=1<<n;dp=array.array('Q',[0])*(M*n);h=[0]*M
 for v in range(n):dp[(1<<v)*n+v]=1
 for S in range(1,M):
  offset=S*n;h[S]=sum(dp[offset+v] for v in vertices(S))
  for v in vertices(S):
   z=dp[offset+v]
   if not z:continue
   for w in vertices(A[v]&((M-1)^S)):dp[(S|(1<<w))*n+w]+=z
 return h
def raw(A,p):return sum(((A[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(itertools.combinations(range(len(p)),2)))
def motif_lookup(records,n):
 lookup={}
 for d in records:
  key=(tuple(d['scores']),d['canonical_upper_bits']);A=[0]*n
  for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
   u,v=(i,j) if key[1]>>k&1 else(j,i);A[u]|=1<<v
  need(tuple(sorted(x.bit_count() for x in A))==key[0],'motif representative score check')
  for p in itertools.permutations(range(n)):
   code=raw(A,p);need(code not in lookup or lookup[code]==key,'distinct motif orbits')
   lookup[code]=key
 need(len(lookup)==1<<(n*(n-1)//2),'complete labelled motif universe')
 return lookup
def motif_hist(A,n,lookup):return Counter(lookup[raw(A,V)] for V in itertools.combinations(range(len(A)),n))
def expected_hist(rows):return Counter({(tuple(x['scores']),x['canonical_upper_bits']):x['count'] for x in rows})
def main():
 p=argparse.ArgumentParser();p.add_argument('--producer',type=Path,default=HERE.parent/'session38-root');p.add_argument('--write',action='store_true');args=p.parse_args()
 producer=args.producer;r=json.loads((producer/'transitive_profile_receipt.json').read_text());seeds=[];extensions=[]
 lookup5=motif_lookup(r['cases'][0]['global_motifs5'],5);lookup4=motif_lookup(r['cases'][0]['neighborhoods'][0]['motifs4'],4)
 for index,case in enumerate(r['cases']):
  A=tournament(15,case['steps']);need(A==case['rows'],'literal seed rows');f,roots,contained=transitive_data(A)
  need(list(f)==case['F'],'seed all-subset transitive polynomial');need(sorted(map(list,roots))==case['root_outneighbor_F'],'all seed root out profiles')
  need(all(c==contained[0] for c in contained),'vertex transitive rooted containment profile')
  need(motif_hist(A,5,lookup5)==expected_hist(case['global_motifs5']),'global five motif census')
  hout=[];r5=[];B=[]
  for v,rec in enumerate(case['neighborhoods']):
   U=induced(A,list(vertices(A[v])));h=sum(all(U[p[i]]>>p[i+1]&1 for i in range(6)) for p in itertools.permutations(range(7)));hout.append(h)
   need(h==rec['H'],'literal neighborhood Hamilton path permutations')
   need(motif_hist(U,4,lookup4)==expected_hist(rec['motifs4']),'local four motif census');need(motif_hist(U,5,lookup5)==expected_hist(rec['motifs5']),'local five motif census')
   regular=sum(all(sum((U[u]>>w)&1 for w in V)==2 for u in V) for V in itertools.combinations(range(7),5));r5.append(regular);need(regular==index,'direct regular five neighborhood count')
   hist=Counter((U[u]&U[w]).bit_count() for u,w in itertools.combinations(range(7),2));need(sorted(map(list,hist.items()))==sorted(rec['common_out_histogram']),'local common-out histogram');B.append(sum(math.comb(k,2)*z for k,z in hist.items()))
  need(B==[10]*15,'local B scalar')
  six=0
  for V in itertools.combinations(range(15),6):
   S=sum(1<<v for v in V)
   sources=[v for v in V if (A[v]&S).bit_count()==5]
   if sources:
    rest=S^(1<<sources[0])
    six+=all((A[v]&rest).bit_count()==2 for v in vertices(rest))
  need(six==sum(r5)==case['source_plus_regular_five_count'],'literal global source plus regular five')
  h=ham_counts(A);F2=sum(h[S]*h[((1<<15)-1)^S] for S in range(1,(1<<15)-1));need(h[-1]==case['H'],'independent forward global Hamilton DP');need(F2==case['ordered_two_path_covers'],'independent ordered two path covers')
  seeds.append({'F':list(f),'root_F':list(roots[0]),'root_containing_F':list(contained[0]),'H':h[-1],'H_out':hout,'regular5_out':r5,'source_plus_regular5_count':sum(r5),'F2':F2})
  for rr in range(1,5):
   L=lift(A,rr);ff,gg,cc=transitive_data(L);P=bp(rr)
   formula=add(multiply(add(P,(0,1)),f),shifted((0,)+P[1:]));need(ff==formula,'cyclic global F identity')
   predicted=[multiply((1,1),g) for g in roots]+[P]+[multiply(bp(j),f) for j in range(rr)]
   need(sorted(gg)==sorted(predicted),'cyclic root-out identity')
   # A further independently checked, precisely typed contained-root identity.
   cc_predict=[multiply(add(P,(0,1)),c) for c in contained]+[shifted(add(f,(0,)+P[1:]))]+[multiply(shifted(bp(rr-1)),add(f,(0,1)))]*rr
   need(cc==cc_predict,'cyclic root-containing identity')
   extensions.append({'seed':index,'r':rr,'order':len(L),'all_subsets_checked':1<<len(L),'F':list(ff),'root_profile_hash':hashlib.sha256(json.dumps(sorted(gg),separators=(',',':')).encode()).hexdigest()})
 need(seeds[0]['F']==seeds[1]['F'] and seeds[0]['root_F']==seeds[1]['root_F'] and seeds[0]['root_containing_F']==seeds[1]['root_containing_F'],'complete transitive profile seed collision')
 need(seeds[0]['H']!=seeds[1]['H'] and seeds[0]['H_out']!=seeds[1]['H_out'],'Hamilton information differs')
 for rr in range(1,5):
  pair=[e for e in extensions if e['r']==rr];need(pair[0]['F']==pair[1]['F'] and pair[0]['root_profile_hash']==pair[1]['root_profile_hash'],'literal lifted profiles collide')
 out={'scope':'Independent exact finite audit plus mathematical review of all-size polynomial identities; no producer imports.','seed_audit':seeds,'literal_lift_audit':extensions,'checks':dict(sorted(count.items())),'producer_sha256':{name:hashlib.sha256((producer/name).read_bytes()).hexdigest() for name in ('check_transitive_profile.py','transitive_profile_receipt.json','TRANSITIVE_PROFILE_OBSTRUCTION.md')}}
 dest=HERE/'receipt.json'
 if args.write:dest.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(json.loads(dest.read_text())==out,'frozen audit reproduction')
 print(json.dumps({'status':'PASS','checks':sum(count.values()),'F':seeds[0]['F'],'root_F':seeds[0]['root_F'],'root_containing_F':seeds[0]['root_containing_F'],'H_out':[s['H_out'][0] for s in seeds],'regular5_out':[s['regular5_out'][0] for s in seeds],'literal_lift_orders':sorted({e['order'] for e in extensions})},indent=2))
if __name__=='__main__':main()
