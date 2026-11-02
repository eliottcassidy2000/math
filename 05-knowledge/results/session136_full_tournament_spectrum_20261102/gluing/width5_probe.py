#!/usr/bin/env python3
"""Optional exact width<=5 grammar at H<=127; imports companion only.
Strong-order pruning uses inherited Moon--Busch f(11)=225>127.
"""
from module_contexts import *

def prime(a):
 n=len(a);full=(1<<n)-1
 for m in range(1,full):
  if m.bit_count()<2:continue
  if all((a[v]&m) in (0,m) for v in range(n) if not(m>>v&1)):return False
 return True

def iso(a,b):
 return any(all(bool(a[i]>>j&1)==bool(b[p[i]]>>p[j]&1) for i in range(len(a)) for j in range(len(a))) for p in permutations(range(len(a))))

def prime5_reps():
 reps=[];multip=[]
 for mask in range(1024):
  a=graph(5,mask)
  if not prime(a):continue
  hit=next((i for i,(_,b) in enumerate(reps) if iso(a,b)),None)
  if hit is None:reps.append((mask,a));multip.append(1)
  else:multip[hit]+=1
 require(len(reps)==3,'prime5 classification')
 return reps,multip

def word_counter(q):
 @lru_cache(None)
 def end(c,last):
  if not c[last]:return 0
  if sum(c)==1:return 1
  d=list(c);d[last]-=1;d=tuple(d)
  return sum(end(d,v) for v in range(len(q)) if q[v]>>last&1)
 def W(c):return sum(end(c,v) for v in range(len(q)))
 return W

def valuation(W,ps):
 fs=[tuple(factorial(k)*p[k] for k in range(len(p))) for p in ps]
 ans=0
 for c in product(*(range(1,len(p)) for p in ps)):
  v=W(c)
  for i,k in enumerate(c):v*=fs[i][k]
  ans+=v
 return ans

def tuples_budget(small,arity,budget):
 def go(out,left,k):
  if not k:
   yield tuple(out);return
  for s in small:
   n=len(s[1])-1
   if n<=left-(k-1):yield from go(out+[s],left-n,k-1)
 yield from go([],budget,arity)

def core_normalize(a,N,capH):
 if capH==128:return substitute(C3,[transitive(7),(0,),(0,)])
 n=len(a);reach=list(a)
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 keep=[i for i in range(n) if reach[i]>>i&1]
 core=tuple(sum(1<<j for j,v in enumerate(keep) if a[u]>>v&1) for u in keep)
 if not core:return transitive(N)
 return substitute(JOIN,[transitive(max(0,N-len(core))),core]) if len(core)<N else core

def main():
 reps,multip=prime5_reps();words=[word_counter(q) for _,q in reps]
 N=7;B=127;beta=128;one=('s',(0,1));states={one};graphs={one:(0,)};rounds=0
 checks=0;byq=[0]*3
 while True:
  old=tuple(sorted(states,key=str));small=[s for s in old if s[0]=='s'];new=set();ng={}
  def H(s):return s[1][1] if s[0]=='s' else s[1]
  def add(s,a):
   if s not in states and s not in new:
    if callable(a):a=a()
    if s[0]=='l':a=core_normalize(a,N,s[1])
    new.add(s);ng[s]=a
  for x,y in product(old,repeat=2):
   if x[0]=='s' and y[0]=='s' and len(x[1])+len(y[1])-2<N:
    p=join_profile(x[1],y[1]);a=substitute(JOIN,[graphs[x],graphs[y]])
    require(p==profile(a),'join profile');add(('s',p),a)
   else:add(('l',min(beta,H(x)*H(y))),lambda:substitute(JOIN,[graphs[x],graphs[y]]))
  for xyz in tuples_budget(small,3,10):
   ps=[s[1] for s in xyz];n=sum(len(p)-1 for p in ps);v=cyclic_value(ps)
   if n<N:
    a=substitute(C3,[graphs[s] for s in xyz]);p=profile(a);require(p[1]==v,'C3');add(('s',p),a)
   else:add(('l',min(beta,v)),lambda:substitute(C3,[graphs[s] for s in xyz]))
  for i,(mask,q) in enumerate(reps):
   for ss in tuples_budget(small,5,10):
    ps=[s[1] for s in ss];n=sum(len(p)-1 for p in ps);v=valuation(words[i],ps);byq[i]+=1
    if n<N:
     a=substitute(q,[graphs[s] for s in ss]);p=profile(a);require(p[1]==v,'Q5 small');add(('s',p),a);checks+=1
    else:add(('l',min(beta,v)),lambda:substitute(q,[graphs[s] for s in ss]))
  if any(s[0]=='l' for s in old):add(('l',beta),substitute(C3,[transitive(7),(0,),(0,)]))
  if not new:break
  states.update(new);graphs.update(ng);rounds+=1
 vals=sorted({(s[1][1] if s[0]=='s' else s[1]) for s in states}-{beta})
 witnesses={}
 for v in vals:
  a=min((graphs[s] for s in states if (s[1][1] if s[0]=='s' else s[1])==v),key=lambda a:(len(a),a))
  require(hp(a)==v,'direct width5 certificate');witnesses[v]=a
 outside=graph(6,146)
 require(hp(outside)==31 and direct_profile(outside)[1]==31,'width6 H31')
 print(json.dumps({'positive_witnesses':witnesses,'max_witness_order':max(map(len,witnesses.values())),'H31_width6_witness':outside,'status':'PASS','B':B,'N':N,'strong_order_cutoff':10,'inherited_floor_at_11':225,
  'prime5_masks':[m for m,q in reps],'prime5_multiplicities':multip,
  'rounds':rounds,'states':len(states),'small_by_order':{n:sum(s[0]=='s' and len(s[1])-1==n for s in states) for n in range(1,N)},
  'evaluations_by_prime_quotient':byq,'direct_small_checks':checks,'attained':vals,'missing_odds':[h for h in range(1,128,2) if h not in vals]},indent=2,sort_keys=True))
if __name__=='__main__':main()
