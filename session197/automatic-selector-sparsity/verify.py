import json,math,random,collections,itertools

def need(b,s):
 if not b:raise RuntimeError(s)
def f(n):
 v=1
 while n:
  if n%2:n//=2
  else:n//=2;v+=n+n%2
 return v
def odd(n):return n//(n&-n)
def v2(n):return (n&-n).bit_length()-1
def suffix_q(t,S):return sum(1<<j for j in range(t) if not (S>>(t-1-j)&1))
def coeff(P,t,S,h,e,R):
 q=suffix_q(t,S);K=h*q+R*(1<<t)
 return (h+K)*(P+1),h*(P+1),h+K,h*(f(P)+f(S)+1)+R*(S+1)-e*((S+1)%2)
def value(P,t,S,a,b,h,e,R):
 x=((P+1)<<(a+b))-(1<<b);n=(x<<t)+S
 J=h*f(n)+R*(n+1)-e*((n+1)%2)
 A,B,C,D=coeff(P,t,S,h,e,R)
 need(J==(A<<(a+b))-(B<<a)-(C<<b)+D,'exact arbitrary prefix/suffix rectangle')
 need((A<<(a+b))+D-J==(B<<a)+(C<<b),'positive recovered response')
 need((A<<(a+b)).bit_length()-n.bit_length()==A.bit_length()-P.bit_length()-t,'regular auxiliary length offset')
 return J

def run(delta,w,state=0):
 for b in w:state=delta[state][int(b)]
 return state
def path(delta,start,targets,allowed=None):
 Q=collections.deque([(start,'')]);seen={start}
 while Q:
  v,w=Q.popleft()
  if v in targets:return w
  for b in range(2):
   u=delta[v][b]
   if u not in seen and (allowed is None or u in allowed):seen.add(u);Q.append((u,w+str(b)))
 raise RuntimeError('path does not exist')
def reach(delta,start):
 seen={start};Q=[start]
 for v in Q:
  for u in delta[v]:
   if u not in seen:seen.add(u);Q.append(u)
 return seen
def period(delta,state,bit):
 seen={};seq=[];v=state
 while v not in seen:seen[v]=len(seq);seq.append(v);v=delta[v][bit]
 mu=seen[v];p=len(seq)-mu;A=mu if mu else p
 return A,p,run(delta,str(bit)*A,state)
def main():
 formulas=suffixes=bands=0
 for P in range(1,24):
  for t in range(5):
   for S in range(1<<t):
    q=suffix_q(t,S)
    for x in range(0,40,2):
     need(f((x<<t)+S)==f(x)+q*x+f(S)-1,'even-prefix suffix identity');suffixes+=1
    for h,e,R in ((1,1,1),(3,1,3),(5,1,7),(13,2,18)):
     for a in range(1,9):
      for b in range(1,9):value(P,t,S,a,b,h,e,R);formulas+=1
 for B in range(1,33):
  for C in range(1,33):
   for delta0 in (-9,-4,-1,0,1,4,9):
    H=odd((B<<max(delta0,0))+(C<<max(-delta0,0)));ell=H.bit_length()-1
    lo=-max(0,v2(B)-v2(C)+ell);hi=max(0,v2(C)-v2(B)+ell)
    need(lo<=delta0<=hi,'chosen difference inside finite band')
    E=set()
    for d in range(-48,49):
     ok=odd((B<<max(d,0))+(C<<max(-d,0)))<=H
     if ok:need(lo<=d<=hi,'odd-part support bound');E.add(d)
    need(delta0 in E,'nonempty diagonal band');bands+=1
 rng=random.Random(1971545);rectangles=leaks=empty=rectpoints=leakblocks=0
 for nstates in range(2,8):
  for trial in range(120):
   delta=[[rng.randrange(nstates),rng.randrange(nstates)] for _ in range(nstates)]
   accepting={i for i in range(nstates) if rng.randrange(2)}
   if trial%2==0:delta[-1]=[nstates-1,nstates-1];accepting.discard(nstates-1)
   reachable=reach(delta,0);live={i for i in reachable if reach(delta,i)&accepting}
   if not live:empty+=1;continue
   components=[];todo=set(live)
   while todo:
    v=min(todo);component={u for u in todo if u in reach(delta,v) and v in reach(delta,u)}
    todo-=component;components.append(component)
   complete=[S for S in components if all(u in S for v in S for u in delta[v])]
   if complete:
    component=complete[0];q=min(component);p='1'+path(delta,0,{q});A,t1,q1=period(delta,q,1);B,t0,q2=period(delta,q1,0);T=math.lcm(t1,t0);s=path(delta,q2,accepting,component)
    P=int(p,2);t=len(s);S=int(s,2) if s else 0;AA,BB,CC,DD=coeff(P,t,S,1,1,1)
    d0=A-B;H=odd((BB<<max(d0,0))+(CC<<max(-d0,0)));ell=H.bit_length()-1
    lo=-max(0,v2(BB)-v2(CC)+ell);hi=max(0,v2(CC)-v2(BB)+ell)
    E={d for d in range(lo,hi+1) if odd((BB<<max(d,0))+(CC<<max(-d,0)))<=H}
    for i in range(6):
     for j in range(6):
      a=A+T*i;b=B+T*j;w=p+'1'*a+'0'*b+s
      need(run(delta,w[1:]) in accepting,'complete-component rectangle accepted')
      value(P,t,S,a,b,1,1,1)
      need((odd((BB<<a)+(CC<<b))<=H)==(a-b in E),'regular odd-part filter equals diagonal band');rectpoints+=1
    step=1+max(abs(d-d0) for d in E)//T
    for i in range(5):
     for j in range(5):need((d0+T*step*(i-j) in E)==(i==j),'explicit distinguishing prefixes')
    rectangles+=1
   else:
    s0=len(live)
    for v in live:
     stay=0
     for bits in itertools.product('01',repeat=s0):stay+=run(delta,bits,v) in live
     need(stay<1<<s0,'every live state has an escaping block');leakblocks+=1
    counts=[0]*nstates;counts[0]=1
    for length in range(31):
     b,r=divmod(length,s0);bound=((1<<s0)-1)**b*(1<<r)
     need(sum(counts[v] for v in accepting)<=bound,'strict entropy block bound')
     nxt=[0]*nstates
     for v,c in enumerate(counts):
      for u in delta[v]:nxt[u]+=c
     counts=nxt
    leaks+=1
 need(rectangles>0 and leaks>0,'both automaton regimes tested')
 primes=[p for p in range(2,401) if all(p%d for d in range(2,math.isqrt(p)+1))];primechecks=0
 for n in range(1,201):
  central=math.comb(2*n,n);need(central*(2*n+1)>=4**n,'central binomial lower bound');product=1;pi=0
  for p in primes:
   if p>2*n:break
   pi+=1;power=p;vp=0;count=0
   while power<=2*n:vp+=(2*n)//power-2*(n//power);count+=1;power*=p
   need(vp<=count and p**vp<=2*n,'prime valuation upper bound');product*=p**vp
  need(product==central and central<=(2*n)**pi,'elementary prime count bound');primechecks+=1
 print(json.dumps({'status':'PASS','rectangle_formula_checks':formulas,'even_suffix_checks':suffixes,'odd_part_band_controls':bands,'complete_component_fixtures':rectangles,'rectangle_automaton_points':rectpoints,'leaky_automaton_fixtures':leaks,'escaping_block_controls':leakblocks,'empty_language_fixtures':empty,'elementary_prime_bound_checks':primechecks,'scope':'Polynomial sparsity of regular selected subgraphs for fixed prefixes and finite explicit tags; no unrestricted spectrum omission.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
