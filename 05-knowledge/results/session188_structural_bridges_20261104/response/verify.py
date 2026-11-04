"""Exact quartic cover tomography. Standard library; no removable assertions."""
import hashlib,itertools,json
from collections import defaultdict

def need(ok,message):
 if not ok:raise RuntimeError(message)

def tournament(n,code):
 a=[0]*n
 for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if code>>b&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def opposite(a):return [((1<<len(a))-1)^r^(1<<v) for v,r in enumerate(a)]
def induced(a,vs):return [sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs]
def endpoints(a):
 n=len(a);D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if not D[s][v]:continue
   rem=a[v]&~s
   while rem:
    bit=rem&-rem;rem-=bit;D[s|bit][bit.bit_length()-1]+=D[s][v]
 return D

def hp(a):return sum(endpoints(a)[-1])
def covers(a):
 n=len(a);full=(1<<n)-1;E=endpoints(a);S=endpoints(opposite(a));C=[[0]*n for _ in a]
 for s in range(1,full):
  for u in range(n):
   if E[s][u]:
    for v in range(n):C[u][v]+=E[s][u]*S[full^s][v]
 return C

def extend(a,x,c):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for u in range(n):
  if x>>u&1:b[u]|=1<<n
  else:b[n]|=1<<u
  if c>>u&1:b[u]|=1<<(n+1)
  else:b[n+1]|=1<<u
 return b

def polynomial(a):
 n=len(a);result=defaultdict(int)
 for word in itertools.permutations(range(n+2)):
  terms={0:1}
  for u,v in zip(word,word[1:]):
   if u<n and v<n:
    if not a[u]>>v&1:break
   elif u>=n and v>=n:
    if u!=n:break
   else:
    bit=1<<((v-n)*n+u if u<n else (u-n)*n+v)
    nxt=defaultdict(int)
    for s,c in terms.items():
     if u<n:nxt[s|bit]+=c
     else:nxt[s]+=c;nxt[s|bit]-=c
    terms=nxt
  else:
   for s,c in terms.items():result[s]+=c
 return {s:c for s,c in result.items() if c}

def quartic_normal_form(a):
 n=len(a);out={};pairs=[(1<<i)|(1<<j) for i,j in itertools.combinations(range(n),2)]
 # Full bidegree (2,2) basis, including named zeros.
 for A in pairs:
  for B in pairs:out[A|(B<<n)]=0
 # One shared physical name: independently computed endpoint convolution.
 for i in range(n):
  keep=[v for v in range(n) if v!=i];C=covers(induced(a,keep));index={v:k for k,v in enumerate(keep)}
  for j,k in itertools.permutations(keep,2):
   mon=(1<<i)|(1<<j)|(1<<(n+i))|(1<<(n+k))
   out[mon]=C[index[j]][index[k]]+C[index[k]][index[j]]
 # Disjoint supports: direct enumeration of old three-path covers.
 for word in itertools.permutations(range(n)):
  bad={b for b in range(1,n) if not a[word[b-1]]>>word[b]&1}
  if len(bad)>2:continue
  for b,c in itertools.combinations(range(1,n),2):
   if not bad.issubset((b,c)):continue
   A=(1<<word[b-1])|(1<<word[b]);B=(1<<word[c-1])|(1<<word[c])
   if A&B:continue
   out[A|(B<<n)]+=1;out[B|(A<<n)]+=1
 return out

def finite_difference(a,b,A,B):
 n=len(a);need(n==len(b),'same named ground set');mon=A|(B<<n);acc=0;records=[];s=mon
 while True:
  x=s&((1<<n)-1);c=s>>n;fa=hp(extend(a,x,c));fb=hp(extend(b,x,c));sign=(-1)**(mon.bit_count()-s.bit_count());acc+=sign*(fa-fb)
  records.append([x,c,fa,fb,sign])
  if s==0:break
  s=(s-1)&mon
 return acc,records

def split_root(q):
 n=len(q);a=list(q)+[q[0]];a[0]|=1<<n
 for v in range(1,n):
  if q[v]&1:a[v]|=1<<n
 b=list(a);b[0]^=1<<n;b[n]^=1
 return a,b

def root_over_transitive(m,bits):
 q=[0]*(m+1)
 for i in range(1,m+1):
  for j in range(i+1,m+1):q[i]|=1<<j
  if bits>>(i-1)&1:q[0]|=1<<i
  else:q[i]|=1
 return q

def strong(a):
 full=(1<<len(a))-1
 for v in range(len(a)):
  r=1<<v
  while True:
   t=r
   for u in range(len(a)):
    if r>>u&1:t|=a[u]
   if t==r:break
   r=t
  if r!=full:return False
 return True

def main():
 normal=[];poly_cache={};coordinate_checks=0;literal_checks=0
 # Exhaustive named quartic normal-form controls through old order five.
 for n in range(1,6):
  count=1<<(n*(n-1)//2)
  for code in range(count):
   a=tournament(n,code);p=polynomial(a);q=quartic_normal_form(a)
   need(all(p.get(s,0)==v for s,v in q.items()),'complete quartic cover normal form')
   need(all(v>=0 and v%2==0 for v in q.values()),'quartic positivity and parity')
   need(all(s.bit_count()<=4 for s in p),'full-response degree')
   coordinate_checks+=len(q);poly_cache[n,code]=p
   # Check a noninterpolation attachment assignment by literal subset DP.
   if n<=4 or code%64==0:
    A=sum(1<<v for v in range(n) if v%2==0);B=((1<<n)-1)^A;mask=A|(B<<n)
    need(sum(v for s,v in p.items() if s&mask==s)==hp(extend(a,A,B)),'literal full response');literal_checks+=1
   digest=hashlib.sha256(json.dumps(sorted(q.items()),separators=(',',':')).encode()).hexdigest()
   normal.append({'n':n,'code':code,'rows':a,'quartic_coordinates':len(q),'quartic_sha256':digest})
  for code in range(count):
   p=poly_cache[n,code];r=poly_cache[n,count-1-code]
   need(all(p.get(s,0)==r.get(s,0) for s in set(p)|set(r) if s.bit_count()==4),'whole converse cancels entire quartic')
 # Twin identity in every rooted quotient through order four, and 20 fixed larger controls.
 twin=[]
 rooted=[tournament(n,c) for n in range(1,5) for c in range(1<<(n*(n-1)//2))]
 rooted += [tournament(5,(7919*k+17)%1024) for k in range(12)]
 rooted += [tournament(6,(15427*k+31)%32768) for k in range(8)]
 for q in rooted:
  a,b=split_root(q);n=len(a);pa=quartic_normal_form(a);pb=quartic_normal_form(b);C=covers(induced(q,list(range(1,len(q)))))
  checks=[]
  for i,j in itertools.permutations(range(1,len(q)),2):
   mon=(1<<(n-1))|(1<<i)|(1<<n)|(1<<(n+j));got=pa[mon]-pb[mon];want=C[j-1][i-1]-C[i-1][j-1]
   need(got==want,'exterior edge derivative');checks.append([i,j,got])
  twin.append({'quotient_rows':q,'T':a,'T_prime':b,'edge_derivative_coefficients':checks})
 # Exact literal exponential certificates, and all root incidences at m=3.
 family=[]
 cases=[(3,bits) for bits in range(8)]+[(m,1) for m in range(4,8)]
 for m,bits in cases:
  q=root_over_transitive(m,bits);a,b=split_root(q);n=len(a);A=(1<<(n-1))|(1<<1);B=1|(1<<m)
  value,terms=finite_difference(a,b,A,B);need(value==2**(m-2),'unbounded family finite control')
  if bits==1:need(strong(a) and strong(b),'strong family')
  family.append({'m':m,'root_incidence_mask':bits,'quotient_rows':q,'T':a,'T_prime':b,'a_mask':A,'c_mask':B,'coefficient':value,'strong':strong(a),'sixteen_literal_queries':terms})
 # Retain every coefficient of a concrete strong order-five hostile.
 q=root_over_transitive(3,1);a,b=split_root(q);p=polynomial(a);r=polynomial(b);op=polynomial(opposite(a))
 difference={s:p.get(s,0)-r.get(s,0) for s in set(p)|set(r) if p.get(s,0)!=r.get(s,0)}
 whole={s:p.get(s,0)-op.get(s,0) for s in set(p)|set(op) if p.get(s,0)!=op.get(s,0)}
 need(max(s.bit_count() for s in difference)==4,'partial reversal quartic')
 need(max(s.bit_count() for s in whole)<=3,'whole converse cubic boundary')
 full=1<<len(a)
 # Every named query is checked independently on this one graph pair.
 for A in range(full):
  for B in range(full):
   mask=A|(B<<len(a));need(sum(v for s,v in difference.items() if s&mask==s)==hp(extend(a,A,B))-hp(extend(b,A,B)),'hostile full truth table')
 return {'status':'PASS','scope':'Complete quartic normal-form census through order five; not a reconstruction census or hidden-core exclusion.','normal_form_tournaments':len(normal),'quartic_coordinate_checks':coordinate_checks,'literal_response_checks':literal_checks,'normal_form_records':normal,'twin_identity_cases':len(twin),'twin_identity_records':twin,'transitive_exterior_family':family,'hostile':{'T':a,'T_prime':b,'full_query_checks':full*full,'omission_convention':'All unlisted coefficients are zero in the entire Boolean basis.','partial_reversal_difference':[[s&31,s>>5,v] for s,v in sorted(difference.items())],'whole_converse_difference':[[s&31,s>>5,v] for s,v in sorted(whole.items())]}}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
