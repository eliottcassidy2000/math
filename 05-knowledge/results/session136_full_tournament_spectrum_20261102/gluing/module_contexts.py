from math import factorial,comb
from itertools import permutations
from functools import lru_cache

def graph(n,mask):
    a=[0]*n
    k=0
    for i in range(n):
      for j in range(i+1,n):
        u,v=(i,j) if mask>>k&1 else (j,i)
        a[u]|=1<<v;k+=1
    return tuple(a)

def h_all(a):
    n=len(a); size=1<<n; dp=[[0]*n for _ in range(size)]; h=[0]*size;h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,size):
      for v in range(n):
        d=dp[s][v]
        if d:
          h[s]+=d
          t=a[v]&~s
          while t:
            b=t&-t;t-=b;dp[s|b][b.bit_length()-1]+=d
    return h

def profile(a):
    n=len(a);h=h_all(a)
    @lru_cache(None)
    def go(s):
      if not s:return (1,)+(0,)*n
      b=s&-s;t=s^b;ans=[0]*(n+1);u=t
      while True:
        block=u|b;p=go(s^block)
        for k in range(n):ans[k+1]+=h[block]*p[k]
        if not u:break
        u=(u-1)&t
      return tuple(ans)
    return go((1<<n)-1)

def direct_profile(a):
    n=len(a);f=[0]*(n+1)
    for p in permutations(range(n)):
      bad=sum(not(a[p[i]]>>p[i+1]&1) for i in range(n-1))
      for k in range(bad+1,n+1):f[k]+=comb(n-1-bad,k-1-bad)
    return tuple(v//factorial(k) for k,v in enumerate(f))

import json
from itertools import product

def require(ok,msg):
 if not ok:raise RuntimeError(msg)

def hp(a):return h_all(a)[-1]
def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def substitute(q,fs):
 off=[0]
 for f in fs:off.append(off[-1]+len(f))
 a=[0]*off[-1]
 for i,f in enumerate(fs):
  for u,row in enumerate(f):
   a[off[i]+u]=row<<off[i]
   for j,g in enumerate(fs):
    if q[i]>>j&1:a[off[i]+u]|=((1<<len(g))-1)<<off[j]
 return tuple(a)

def coefficients(q):
 m=len(q)-1;k=m//2+1;full=(1<<m)-1
 @lru_cache(None)
 def go(mask,last,c):
  ans=[0]*(k+1)
  if mask==full and c:ans[c]+=1
  if c<k and last!=0 and (last<0 or q[last]&1):
   v=go(mask,0,c+1)
   for i in range(k+1):ans[i]+=v[i]
  for j in range(1,m+1):
   b=1<<(j-1)
   if not mask&b and (last<0 or q[last]>>j&1):
    v=go(mask|b,j,c)
    for i in range(k+1):ans[i]+=v[i]
  return tuple(ans)
 return tuple(v*factorial(i) for i,v in enumerate(go(0,-1,0)))

def active(q):
 reach=list(q)
 for k in range(len(q)):
  for i in range(len(q)):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return bool(reach[0]&1)

def cyclic_value(ps):
 fs=[tuple(factorial(k)*v for k,v in enumerate(p)) for p in ps]
 def f(i,k):return fs[i][k] if k<len(fs[i]) else 0
 ans=0
 for r in range(1,min(map(len,ps))):
  ans+=3*f(0,r)*f(1,r)*f(2,r)
  for i in range(3):
   js=[j for j in range(3) if j!=i]
   ans+=f(i,r+1)*f(js[0],r)*f(js[1],r)
   ans+=f(i,r)*f(js[0],r+1)*f(js[1],r+1)
 return ans

def join_profile(p,q):
 ans=[0]*(len(p)+len(q)-1)
 for a in range(1,len(p)):
  for b in range(1,len(q)):
   for t in range(min(a,b)+1):ans[a+b-t]+=p[a]*q[b]*comb(a,t)*comb(b,t)*factorial(t)
 return tuple(ans)
C3=(2,4,1);JOIN=(2,0)

def grammar_closure(B,seed_order=1):
 N=1
 while 2**N+1<=B:N+=1
 beta=B+1
 def state(a):
  p=profile(a)
  return ('s',p) if len(a)<N else ('l',min(beta,p[1]))
 single=(0,);ss=state(single);states={ss:'1'};sg={ss:single} if 1<N else {};rounds=0
 for n in range(2,seed_order+1):
  for mask in range(1<<(n*(n-1)//2)):
   a=graph(n,mask);st=state(a)
   if st not in states:
    states[st]='S'+str(n)+'_'+str(mask)
    if st[0]=='s':sg[st]=a
 while True:
  old=tuple(states);small=[s for s in old if s[0]=='s'];new={};ng={}
  def remember(s,e,a=None):
   if s not in states and s not in new:
    new[s]=e
    if a is not None:ng[s]=a
  def H(s):return s[1][1] if s[0]=='s' else s[1]
  for x,y in product(old,repeat=2):
   if x[0]=='s' and y[0]=='s' and len(x[1])+len(y[1])-2<N:
    p=join_profile(x[1],y[1]);a=substitute(JOIN,[sg[x],sg[y]])
    require(p==profile(a),'join polynomial')
    remember(('s',p),'J('+states[x]+','+states[y]+')',a)
   else:remember(('l',min(beta,H(x)*H(y))),'J('+states[x]+','+states[y]+')')
  for xyz in product(small,repeat=3):
   ps=[s[1] for s in xyz];n=sum(len(p)-1 for p in ps);value=cyclic_value(ps)
   e='C('+','.join(states[s] for s in xyz)+')'
   if n<N:
    a=substitute(C3,[sg[s] for s in xyz]);p=profile(a)
    require(p[1]==value,'cyclic small');remember(('s',p),e,a)
   else:remember(('l',min(beta,value)),e)
  if any(s[0]=='l' for s in old):
   firstlarge=next(s for s in old if s[0]=='l')
   remember(('l',beta),'C('+states[firstlarge]+',1,1)')
  if not new:break
  states.update(new);sg.update(ng);rounds+=1
 values=sorted({H(s) for s in states}-{beta})
 return {'bound':B,'seed_order':seed_order,'N':N,'rounds':rounds,'state_count':len(states),
  'small_states_by_order':{n:sum(s[0]=='s' and len(s[1])-1==n for s in states) for n in range(1,N)},
  'attained':values,'missing_odd':[h for h in range(1,B+1,2) if h not in values],
  'witnesses':{h:next(states[s] for s in states if H(s)==h) for h in values}}

def parse_witness(expr):
 def parse(i):
  if expr[i]=='1':return (0,),i+1
  if expr[i]=='S':
   j=i+1
   while j<len(expr) and (expr[j].isdigit() or expr[j]=='_'):j+=1
   n,mask=map(int,expr[i+1:j].split('_'));return graph(n,mask),j
  op=expr[i];require(op in 'JC' and expr[i+1]=='(','witness syntax')
  args=[];j=i+2
  while True:
   a,j=parse(j);args.append(a)
   if expr[j]==')':j+=1;break
   require(expr[j]==',','witness comma');j+=1
  require(len(args)==(2 if op=='J' else 3),'witness arity')
  return substitute(JOIN if op=='J' else C3,args),j
 a,j=parse(0);require(j==len(expr),'witness trailing');return a

def verify_witnesses(closure):
 orders=[]
 for h,expr in closure['witnesses'].items():
  a=parse_witness(expr);require(hp(a)==h,'direct attained certificate')
  orders.append(len(a))
 return {'verified':len(orders),'max_order':max(orders)}

def width_controls():
 reps=[transitive(4),substitute(JOIN,[C3,(0,)]),substitute(JOIN,[(0,),C3]),substitute(C3,[transitive(2),(0,),(0,)])]
 def iso(a,b):
  for p in permutations(range(len(a))):
   if all(bool(a[i]>>j&1)==bool(b[p[i]]>>p[j]&1) for i in range(len(a)) for j in range(len(a))):return True
  return False
 counts=[0]*4
 for mask in range(64):
  hits=[i for i,b in enumerate(reps) if iso(graph(4,mask),b)]
  require(len(hits)==1,'four-vertex classification');counts[hits[0]]+=1
 q=graph(5,10);fs=[(0,),transitive(2),(0,),(0,),(0,)];a=substitute(q,fs)
 literal=sum(all(a[p[i]]>>p[i+1]&1 for i in range(5)) for p in permutations(range(6)))
 require(hp(a)==literal==19,'width-five H19 witness')
 return {'four_vertex_recipe_counts':counts,'H19_recipe':{'quotient_n':5,'quotient_mask':10,'quotient_H':hp(q),'blown_vertex':1,'adjacency':a,'H':literal}}

def main():
 pc=0;pbank={}
 for n in range(1,6):
  for mask in range(1<<(n*(n-1)//2)):
   a=graph(n,mask);p=profile(a)
   require(p==direct_profile(a),'independent profile')
   require((p[2] if n>=2 else 0)>=2**(n-1)-1,'pc2 floor');pc+=1
   if n<=4:pbank[p]=a
 contexts=direct=threshold=0
 for qn in range(1,6):
  for mask in range(1<<(qn*(qn-1)//2)):
   q=graph(qn,mask);co=coefficients(q);ac=active(q)
   require(co[1]==hp(q),'coefficient1')
   require((len(co)>2 and co[2]>0)==ac,'a2 cyclic iff')
   if not ac:require(not any(co[2:]),'passive higher')
   contexts+=1
   for p,a in pbank.items():
    n=len(a);value=sum(v*(p[k] if k<len(p) else 0) for k,v in enumerate(co))
    actual=hp(substitute(q,[a]+[(0,)]*(qn-1)))
    require(value==actual,'transfer')
    if ac:require(value>=2**n+1,'module floor')
    direct+=1
    for B in (1,3,7,15,31,127):
     N=1
     while 2**N+1<=B:N+=1
     if n>=N:
      predicted=B+1 if ac else min(B+1,co[1]*min(B+1,p[1]))
      require(predicted==min(B+1,value),'large threshold');threshold+=1
 probes=[]
 for s in range(1,6):
  q=substitute(C3,[(0,),transitive(s),transitive(s)]);co=coefficients(q)
  require(co[-1]==factorial(s+1)*factorial(s)**2,'triangular top')
  probes.append({'s':s,'coefficients':co[1:]})
 padding=[]
 for n in range(1,9):
  v=hp(substitute(C3,[transitive(n),(0,),(0,)]));require(v==2**n+1,'padding hostile');padding.append(v)
 hostile=[]
 for mask in (40,76):
  a=graph(5,mask);p=profile(a)
  hostile.append({'mask':mask,'profile':p[1:],'triangle_response':hp(substitute(C3,[a,(0,),(0,)]))})
 result={'status':'PASS','profile_vs_permutation_checks':pc,'all_labelled_contexts_orders_1_to_5':contexts,
  'distinct_factor_profiles_through_4':len(pbank),'direct_transfer_checks':direct,
  'large_state_threshold_checks':threshold,'triangular_probes':probes,
  'H1_padding_triangle_responses':padding,'inherited_same_H_hostile_THM3134':hostile,
  'join_C3_grammar':grammar_closure(127),
  'five_seed_join_C3_grammar':grammar_closure(127,5),
  'substitution_width_controls':width_controls(),
  'outside_first_grammar_H13':{'n':5,'mask':41,'adjacency':graph(5,41),'H':hp(graph(5,41))}}
 for key in ('join_C3_grammar','five_seed_join_C3_grammar'):
  result[key]['direct_witness_audit']=verify_witnesses(result[key])
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
