#!/usr/bin/env python3
"""Exact hybrid-induction certificate. Standard library only; all gates survive -O."""
import hashlib,itertools,json
from pathlib import Path

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def rows_for_word(w,apices=2):
 rows=[0]*(len(w)+apices)
 if apices==2:rows[0]=2
 for j,c in enumerate(w,apices):
  for a in range(apices):
   if c>>a&1:rows[a]|=1<<j
   else:rows[j]|=1<<a
  for i in range(apices,j):rows[i]|=1<<j
 return rows

def induced(rows,ix):
 return [sum(((rows[u]>>v)&1)<<j for j,v in enumerate(ix)) for u in ix]

def dp(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)]
 if not n:return 1
 for v in range(n):d[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    old=mask^(1<<v)
    if old:d[mask][v]=sum(d[old][u] for u in range(n) if rows[u]>>v&1)
 return sum(d[-1])

def permutation_count(rows):
 return sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in itertools.permutations(range(len(rows))))

def strong(rows):
 for source in range(len(rows)):
  seen=1<<source
  while True:
   new=seen
   for u,r in enumerate(rows):
    if seen>>u&1:new|=r
   if new==seen:break
   seen=new
  if seen!=(1<<len(rows))-1:return False
 return True

def valid(rows):
 n=len(rows)
 return all(not(rows[i]>>i&1) for i in range(n)) and all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n))

# First independent algorithm: S191's 36-state run-assignment transfer.
TRANS=[[[] for _ in range(36)] for _ in range(4)];ACCEPT=[]
for o,x,y,z in itertools.product(range(2),range(3),range(3),range(2)):
 i=o*18+x*6+y*2+z
 if x!=2 and (y==1 or (y==0 and o==0)):ACCEPT.append(i)
 for c in range(4):
  p=c>>o&1;q=c>>(1-o)&1
  TRANS[c][i].append(o*18+(2 if p else 1)*6+y*2+z)
  if y or p:TRANS[c][i].append(o*18+x*6+(2 if q else 1)*2+z)
  if z or q:TRANS[c][i].append(o*18+x*6+y*2+1)
INIT=(1,)+(0,)*17+(1,)+(0,)*17

def step(s,c):
 out=[0]*36
 for i,h in enumerate(s):
  for j in TRANS[c][i]:out[j]+=h
 return tuple(out)

# Second algorithm: literally enumerate each possible Hamilton word with the
# chain increasing within the three slots, turn its required incidences into
# a Boolean monomial, and evaluate all assignments by subset zeta transform.
# It shares neither statuses nor the transfer recurrence with the first.
def literal_all_counts(m):
 a=[0]*(4**m)
 for order in ((0,1),(1,0)):
  for slots in itertools.product(range(3),repeat=m):
   path=[]
   for slot in range(3):
    path.extend(j+2 for j,s in enumerate(slots) if s==slot)
    if slot<2:path.append(order[slot])
   yes=no=0;legal=True
   for u,v in zip(path,path[1:]):
    if u<2 and v<2:
     if (u,v)!=(0,1):legal=False;break
    elif u>=2 and v>=2:
     if u>v:legal=False;break
    elif u<2:yes|=1<<(2*(v-2)+u)
    else:no|=1<<(2*(u-2)+v)
   if not legal or yes&no:continue
   sub=no
   while True:
    a[yes|sub]+= -1 if sub.bit_count()%2 else 1
    if not sub:break
    sub=(sub-1)&no
 for b in range(2*m):
  bit=1<<b
  for mask in range(1<<(2*m)):
   if mask&bit:a[mask]+=a[mask^bit]
 return a

def join(*parts):
 rows=[];offset=0;total=sum(map(len,parts))
 for part in parts:
  later=((1<<total)-1)^((1<<(offset+len(part)))-1)
  rows.extend((r<<offset)|later for r in part);offset+=len(part)
 return rows

def rooted_profile(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    old=mask^(1<<v)
    if old:d[mask][v]=sum(d[old][u] for u in range(n) if rows[u]>>v&1)
 h=[sum(a) for a in d];h[0]=1;full=(1<<n)-1
 return [sum(d[s][v]*h[full^s] for s in range(1,1<<n)) for v in range(n)]

def literal_profile(rows):
 n=len(rows);f=[0]*n
 for p in itertools.permutations(range(n)):
  prefix=[True]*n;suffix=[True]*n
  for j in range(1,n):prefix[j]=prefix[j-1] and bool(rows[p[j-1]]>>p[j]&1)
  for j in range(n-2,-1,-1):suffix[j]=suffix[j+1] and bool(rows[p[j]]>>p[j+1]&1)
  for cut in range(1,n+1):
   if prefix[cut-1] and (cut==n or suffix[cut]):f[p[cut-1]]+=1
 return f

def wrapper(core,root,k):
 n=len(core);w=n+k;r=list(core)+[0]*(k+1)
 for z in range(n,n+k):
  r[z]=((1<<n)-1)|sum(1<<j for j in range(z+1,n+k))
  r[w]|=1<<z
 for i in range(n):
  if i==root:r[i]|=1<<w
  else:r[w]|=1<<i
 return r

def prime(n):return n>=2 and all(n%d for d in range(2,int(n**.5)+1))
def atom(n):
 return (prime(n) and n!=7) or (n%7==0 and prime(n//7) and n//7>=5) or n in (63,343)
def factors(n):
 out=[];p=3
 while p*p<=n:
  while n%p==0:out.append(p);n//=p
  p+=2
 if n>1:out.append(n)
 return out

def atom_factorization(n):
 a=0
 while n%7==0:n//=7;a+=1
 f=factors(n);out=[]
 if a>=2:
  if a%2:out.append(343);a-=3
  out.extend([49]*(a//2));a=0
 if a==1:
  p=next((p for p in f if p>=5),None)
  if p is not None:out.append(7*p);f.remove(p)
  else:need(f.count(3)>=2,'excluded 7 or 21');out.append(63);f.remove(3);f.remove(3)
 return out+f

def main():
 need(5**7>41**3 and 5**10>127**3,'Busch cutoffs')
 arrays=[literal_all_counts(m) for m in range(9)]
 counts=[0]*9;parent=[];target=[]
 def dfs(w,s,code):
  h=sum(s[i] for i in ACCEPT);m=len(w);counts[m]+=1
  need(h==arrays[m][code],'independent complete-spectrum mismatch')
  if h in (41,127):
   r=rows_for_word(w)
   if strong(r):
    need(h==dp(r),'target subset-DP count')
    if h==41 and m<=5:parent.append(w)
    if h==127:target.append(w)
  if m<8:
   for c in range(4):dfs(w+(c,),step(s,c),code|(c<<(2*m)))
 dfs((),INIT,0)
 need(counts==[4**m for m in range(9)],'complete universe')
 expected=['3200','3220','3320','10020','10201','10310','10321','11310','11321','12011','12031','13101','13231','20012','20030','20132','21332','30332','31011','31031','32331']
 need(sorted(parent,key=lambda w:(len(w),w))==[tuple(map(int,w)) for w in expected],'all strong H41 presentations')
 need(len(target)==40,'all strong H127 presentations')
 need(all(len(w)>=4 for w in parent),'no small H41 core')
 onechecks=0
 for m in range(7):
  for w in itertools.product(range(2),repeat=m):
   need(dp(rows_for_word(w,1))!=41,'one-apex H41 exclusion');onechecks+=1
 subsetchecks=0;subset_records=[]
 for w in target:
  r=rows_for_word(w);hist={}
  for m in (6,7):
   for ix in itertools.combinations(range(len(r)),m):
    h=dp(induced(r,ix));need(h!=41,'H41 induced parent found');subsetchecks+=1
    hist[h]=hist.get(h,0)+1
  subset_records.append({'word':''.join(map(str,w)),'order':len(r),'subset_histogram':hist})
 need(subsetchecks==3456,'subset universe')
 # Physical reset: reverse a single old arc before three specified births.
 reset=[
  [14,48,26,50,33,5],
  [14,48,58,50,33,1],
  [14,48,122,50,33,1,59],
  [14,48,250,50,33,1,187,59],
  [14,48,506,306,33,1,443,315,51],
 ]
 expected_h=[41,23,37,63,127]
 for r,h in zip(reset,expected_h):
  need(valid(r) and strong(r),'reset validity and strength')
  need(dp(r)==permutation_count(r)==h,'literal reset verification')
 need(reset[0][2]^(1<<5)==reset[1][2] and reset[0][5]^(1<<2)==reset[1][5] and all(reset[0][i]==reset[1][i] for i in (0,1,3,4)),'one old-arc reversal')
 for old,new in zip(reset[1:],reset[2:]):need(induced(new,range(len(old)))==old,'literal induced birth')
 seed_specs=[(3,'10',1),(5,'100',1),(23,'11000',1),(49,'10212',2),(63,'1001100',1),(67,'12301',2),(343,'1012000',2)]
 seeds=[]
 for h,w,a in seed_specs:
  r=rows_for_word(tuple(map(int,w)),a);need(dp(r)==h and strong(r),'physical seed')
  seeds.append({'H':h,'word':w,'apices':a,'rows':r})
 five=rows_for_word((1,0,0),1);product=join(five,five,five)
 need(dp(product)==125,'three 5 factors')
 # Infinite mixed fan: one old flip, one birth; lower branch factors smaller.
 lift=[[16,69,121,115,98,67,1],[18,68,121,115,98,67,1],[146,68,249,243,226,195,129,2]]
 profiles=[]
 for r,h,F in zip(lift,(41,55,127),(80,105,240)):
  f=rooted_profile(r);need(f==literal_profile(r),'independent rooted covers')
  need(dp(r)==h and f[4]==F and strong(r),'lift core profile')
  profiles.append({'H':h,'F4':F,'all_F':f,'rows':r})
 need(lift[0][0]^2==lift[1][0] and lift[0][1]^1==lift[1][1] and lift[0][2:]==lift[1][2:],'lift one-arc flip')
 need(induced(lift[2],range(7))==lift[1] and lift[2][7]==2,'lift sole-outneighbor birth')
 lift_checks=[]
 for k in range(1,6):
  wrapped=[]
  for r,h,F in zip(lift,(41,55,127),(80,105,240)):
   t=wrapper(r,4,k);need(valid(t) and strong(t) and dp(t)==h+F*2**k,'actual all-clock sample')
   wrapped.append(t)
  need(sum((a^b).bit_count() for a,b in zip(wrapped[0],wrapped[1]))==2,'physical wrapped flip')
  oldix=list(range(7))+list(range(8,8+k))+[8+k]
  need(induced(wrapped[2],oldix)==wrapped[1],'physical wrapped birth')
  N=41+80*2**k;q=25+48*2**k
  need(3*N+2==5*q and q<N and q not in (7,21),'composite lower route')
  need(3*N+4==127+240*2**k,'upper route')
  lift_checks.append({'k':k,'parent':N,'intermediate':55+105*2**k,'upper':3*N+4,'lower':3*N+2,'smaller_factor':q})
 # The same literal flip and sole-old-outneighbor birth is not repeatable.
 again=list(lift[2]);again[0]^=2;again[1]^=1
 nextcore=[a|(0 if i==1 else 1<<8) for i,a in enumerate(again)]+[2]
 need(dp(again)==119 and dp(nextcore)==269 and 269!=3*127+4,'nonreturn control')
 # Positive witness in the forbidden endpoint class and a simple H41 witness.
 controls=[]
 for h,w in [(41,'3200'),(127,'122132')]:
  r=rows_for_word(tuple(map(int,w)));need(dp(r)==h and strong(r),'positive witness')
  controls.append({'H':h,'word':w,'rows':r})
 bound=100001;atomcount=0
 for n in range(1,bound,2):
  if n in (7,21):continue
  f=atom_factorization(n);p=1
  for x in f:need(atom(x),'atom factor');p*=x
  need(p==n,'atom factorization product')
  if atom(n):
   atomcount+=1
   need(f==[n],'irreducible atom factorization')
   if n not in (3,5,23,49,63,67,343):
    need(n%3!=0,'nonseed atom ternary residue')
    q=(n-(2 if n%3==2 else 4))//3
    need(q>=3 and q<n and q%2==1 and q not in (7,21),'prime-lane descent')
 result={'status':'PASS','complete_word_universe':sum(counts),'counts_by_length':counts,'independent_algorithms':['36-state transfer','literal slots + Boolean polynomial zeta'],'strong_H41_words':expected,'strong_H127_word_count':40,'strong_H127_words':[''.join(map(str,w)) for w in target],'one_apex_presentations_through_length6':onechecks,'induced_six_seven_subset_checks':subsetchecks,'target_subset_records':subset_records,'reset_path':[{'H':h,'rows':r} for r,h in zip(reset,expected_h)],'infinite_lift_profiles':profiles,'infinite_lift_samples':lift_checks,'seed_witnesses':seeds,'positive_controls':controls,'product125_rows':product,'arithmetic_odd_bound_exclusive':bound,'arithmetic_atom_count':atomcount}
 raw=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode();Path(__file__).with_name('receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','words':sum(counts),'H41_core_words':len(parent),'H127_core_words':len(target),'induced_subset_checks':subsetchecks,'reset_counts':expected_h,'seed_witnesses':len(seeds),'lift_samples':len(lift_checks),'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
