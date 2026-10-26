#!/usr/bin/env python3
"""Native replay of derived guarded run macros. No downloaded code; stdlib only."""
import itertools, json, math
MACHINES={
 'sandwich':'1RB1LC_1LA1RD_1RA0LC_1LB0RD',
 'bouncer_bell':'1RB0LA_0RC1LA_1RD0RA_0LB1RB',
 'irregular':'1RB0RC_1LC1RA_1RA1LD_0LC0LA',
}
def require(b,msg):
 if not b: raise RuntimeError(msg)
def parse(code):
 return {(q,a):(int(row[3*a]),1 if row[3*a+1]=='R' else -1,row[3*a+2]) for q,row in zip('ABCDEFGHIJKLMNOPQRSTUVWXYZ',code.split('_')) for a in range(len(row)//3)}
TABLES={name:parse(code) for name,code in MACHINES.items()}
def native(name,ones,h,q,steps):
 ones=set(ones)
 for _ in range(steps):
  b,d,q1=TABLES[name][q,int(h in ones)]
  if b: ones.add(h)
  else: ones.discard(h)
  h+=d;q=q1
 return frozenset(ones),h,q

def macro(name,ones,h,q):
 old=frozenset(ones);ones=set(ones)
 if name=='sandwich':
  require(q in 'AB','sandwich state')
  d=1 if q=='A' else -1
  if h not in ones:
   ones.add(h);out=(h+d,'B' if q=='A' else 'A');steps=1;branch=q+'0';k=0
  else:
   # A1 sweeps left; B1 sweeps right.
   d=-d;k=0
   while h+d*(k+1) in ones:k+=1
   for i in range(1,k+1):ones.remove(h+d*i)
   ones.add(h+d*(k+1));out=(h+d*k,q);steps=k+2;branch=q+'1'
 elif name=='bouncer_bell':
  require(q=='A','bouncer_bell state')
  if h in ones:
   k=0
   while h-k in ones:ones.remove(h-k);k+=1
   out=(h-k,'A');steps=k;branch='erase'
  else:
   ones.add(h);k=0;j=h+1
   while j not in ones and j+1 not in ones and j+2 in ones:
    ones.add(j+1);j+=3;k+=1
   if j in ones:out=(j-1,'A');steps=3*k+2;branch='pump1'
   elif j+1 in ones:
    ones.remove(j+1);out=(j+2,'A');steps=3*k+3;branch='pump01'
   else:
    require(j+2 not in ones,'pump residual guard')
    ones.add(j+1);out=(j,'A');steps=3*k+5;branch='pump000'
 elif name=='irregular':
  require(q in 'AC','irregular state')
  if q=='A':
   k=0
   if h in ones:
    ones.remove(h);out=(h+1,'C');steps=1;branch='A1'
   elif h+1 in ones:
    ones.add(h);out=(h+2,'A');steps=2;branch='A01'
   else:
    ones.add(h);ones.add(h+1);out=(h,'C');steps=2;branch='A00'
  else:
   k=0;j=h
   while j in ones and j-1 not in ones:j-=2;k+=1
   if j not in ones:
    ones.add(j);out=(j+1,'A');steps=2*k+1;branch='sweep0'
   else:
    require(j-1 in ones,'sweep residual guard')
    ones.remove(j-1);out=(j-2,'A');steps=2*k+2;branch='sweep11'
 else: raise ValueError(name)
 return frozenset(ones),out[0],out[1],steps,branch,k

def check_macros():
 counts={name:{} for name in MACHINES};cases=0;max_steps=0
 for n in range(11):
  for mask in range(1<<n):
   tape=frozenset(i for i in range(n) if mask>>i&1)
   for h in range(-1,n+1):
    for name,states in [('sandwich','AB'),('bouncer_bell','A'),('irregular','AC')]:
     for q in states:
      t1,h1,q1,s,b,k=macro(name,tape,h,q)
      require(native(name,tape,h,q,s)==(t1,h1,q1),(name,n,mask,h,q,b,k))
      expected={
       'sandwich':{'A0':(1,1),'B0':(1,1),'A1':(k+1,1-k),'B1':(k+1,1-k)},
       'bouncer_bell':{'erase':(k,-k),'pump1':(k+1,k+1),'pump01':(k+2,k),'pump000':(k+2,k+2)},
       'irregular':{'A1':(1,-1),'A01':(1,1),'A00':(2,2),'sweep0':(1,1),'sweep11':(1,-1)},
      }[name][b]
      require((len(tape^t1),len(t1)-len(tape))==expected,'edit and mass formulas')
      counts[name][b]=counts[name].get(b,0)+1;cases+=1;max_steps=max(max_steps,s)
 for k in range(201):
  # Long sweeps whose output edit distance stays 1, and all three pump exits.
  for tail in ['1','01','000']:
   word='0'+'001'*k+tail;tape=frozenset(i for i,x in enumerate(word) if x=='1')
   t1,h1,q1,s,b,kk=macro('bouncer_bell',tape,0,'A')
   require(kk==k,'pump length');require(native('bouncer_bell',tape,0,'A',s)==(t1,h1,q1),'pump replay');cases+=1
  for prefix in ['0','11']:
   word=prefix+'01'*k;tape=frozenset(i for i,x in enumerate(word) if x=='1');h=len(word)-1
   t1,h1,q1,s,b,kk=macro('irregular',tape,h,'C')
   require(kk==k,'sweep length');require(native('irregular',tape,h,'C',s)==(t1,h1,q1),'sweep replay')
   require(len(tape^t1)==1,'sweep one edit');cases+=1
 return {'exhaustive_span_max':10,'exhaustive_states':{'sandwich':'AB','bouncer_bell':'A','irregular':'AC'},'all_head_positions':'-1..n','branch_counts':counts,'long_run_k_range':[0,200],'cases':cases,'max_steps_in_exhaustive_part':max_steps}

def graph(bits):
 m=len(bits);n=3*m
 def edge(i,j):
  if i//3!=j//3:return i//3<j//3
  a,b=i%3,j%3
  if bits[i//3] and {a,b}=={0,2}:return a==2
  return a<b
 return n,edge

def iso_dist(bits1,bits2):
 n,a=graph(bits1);_,b=graph(bits2);best=n*n
 for p in itertools.permutations(range(n)):
  cost=sum(a(i,j)!=b(p[i],p[j]) for i in range(n) for j in range(i+1,n))
  best=min(best,cost)
  if best==0:break
 return best

def check_embedding():
 pairs=0;path_counts=[]
 for m in range(3):
  words=list(itertools.product((0,1),repeat=m))
  for u in words:
   n,e=graph(u)
   c3=sum(e(a,b)==e(b,c)==e(c,a) for a,b,c in itertools.combinations(range(n),3))
   hp=sum(all(e(p[i],p[i+1]) for i in range(n-1)) for p in itertools.permutations(range(n)))
   require(c3==sum(u),'triangle count');require(hp==3**sum(u),'Hamiltonian count')
   path_counts.append({'word':''.join(map(str,u)),'c3':c3,'H':hp})
   for v in words:
    require(iso_dist(u,v)==sum(a!=b for a,b in zip(u,v)),'Hamming isometry');pairs+=1
 return {'max_word_length':2,'word_pairs':pairs,'all_permutations_checked_up_to':6,'invariants':path_counts}

def orbit_samples(limit=10000):
 out={}
 for name in MACHINES:
  ones=frozenset();h=0;q='A';steps=0;calls=0
  while steps<limit:
   t1,h1,q1,s,b,k=macro(name,ones,h,q)
   require(native(name,ones,h,q,s)==(t1,h1,q1),'blank orbit replay')
   ones,h,q=t1,h1,q1;steps+=s;calls+=1
  out[name]={'native_steps':steps,'macro_calls':calls,'head':h,'state':q,'ones':len(ones),'support_min':min(ones),'support_max':max(ones)}
 return out

def decode_scc(n,e):
 reach=[(1<<i)|sum(1<<j for j in range(n) if j!=i and e(i,j)) for i in range(n)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 remaining=set(range(n));comps=[]
 while remaining:
  i=min(remaining);c={j for j in remaining if reach[i]>>j&1 and reach[j]>>i&1}
  comps.append(c);remaining-=c
 return comps

def check_decoding():
 total=0
 for m in range(9):
  for word in itertools.product((0,1),repeat=m):
   n,e0=graph(word)
   e=lambda i,j:e0(n-1-i,n-1-j)
   # Reachability components, then explicitly order their condensation.
   comps=decode_scc(n,e);saved=list(comps)
   comps.sort(key=lambda c:-sum(e(min(c),min(d)) for d in saved if c!=d))
   out=[];zeros=0
   for c in comps+[set()]:
    if len(c)==1:zeros+=1;continue
    require(zeros%3==0,'singleton run not divisible by three')
    out.extend([0]*(zeros//3));zeros=0
    if c:require(len(c)==3,'non-C3 strong block');out.append(1)
   require(tuple(out)==word,'SCC decoder');total+=1
 def twins(word,a,b):
  n,e=graph(word)
  return all(e(a,x)==e(b,x) for x in range(n) if x not in (a,b))
 require(twins((0,0),2,3),'zero blocks merge twins')
 require(not twins((1,0),2,3),'one write destroys cross-boundary twins')
 return {'all_words_through_length':8,'words':total,'labels_reversed_before_decoding':True,'literal_cross_boundary_twins_00_split_under_10':True}

def check_hostiles():
 a=frozenset([1]);b=frozenset([0])
 ma=macro('bouncer_bell',a,0,'A');mb=macro('bouncer_bell',b,0,'A')
 require(len(a)==len(b)==1,'equal scalar source')
 require(ma[:3]!=mb[:3],'source position hostile')
 require(iso_dist((0,1),(1,0))==2,'distinct tournament words')
 require(native('sandwich',set(),0,'A',1)!=native('sandwich',set(),0,'B',1),'missing state hostile')
 # All of 000 is needed for the terminal 4-step B rule.
 t=frozenset([3]);wrong_native=native('bouncer_bell',t,0,'A',5)
 wrong_target=(frozenset([0,2,3]),1,'A')
 require(wrong_native!=wrong_target,'unguarded 000 rule rejected on 001')
 return {'same_H_and_c3':{'words':['01','10'],'H':3,'c3':1,'distance':2,'bouncer_bell_heads_after_macro':[ma[1],mb[1]],'bouncer_bell_clocks':[ma[3],mb[3]]},'missing_machine_state_rejected':True,'drop_terminal_000_guard_rejected':True,'unbounded_clock_distance_one':'C on 0(01)^k: 2k+1 steps, one changed bit'}

if __name__=='__main__':
 print(json.dumps({'machines':MACHINES,'macros':check_macros(),'embedding':check_embedding(),'blank_orbit':orbit_samples(),'hostiles':check_hostiles(),'decoding':check_decoding()},sort_keys=True,indent=2))
