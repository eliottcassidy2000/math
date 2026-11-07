import json, pathlib, subprocess, random
from check_suffix import f,rev
D=pathlib.Path(__file__).resolve().parent
states=[tuple(map(int,l.split())) for l in (D/'CORE_STATES.txt').read_text().splitlines()]
cores=[(8,1,3,6),(2,4,9,19,7),(2,8,3,21,7)]
def value(c,p,r,b):
 h,q,A,G,T,S,E=states[c]
 return (q*p+A)*((2<<b)-rev(r,b))+h*(2*r+1-f(r))-p*G-h*f(p)-T-(p&1)*S-E*(r>>(b-1) if b else 0)
def graph(c,p,r,b):
 if c<16:core,cut=cores[0],c
 elif c<48:core,cut=cores[1],c-16
 else:core,cut=cores[2],c-48
 a=p.bit_length();k=len(core);n=a+k+b;out=[0]*(n+2)
 selected=[]
 for i in range(a):
  if p>>(a-1-i)&1:selected.append(i)
 for i in range(k):
  if cut>>i&1:selected.append(a+i)
 for i in range(b):
  if r>>(b-1-i)&1:selected.append(a+k+i)
 for i in range(n):
  for j in range(i+1,n):
   if a<=i<j<a+k:
    if core[i-a]>>(j-a)&1:out[i]|=1<<j
    else:out[j]|=1<<i
   else:out[i]|=1<<j
  out[i]|=1<<n
  if i in selected:out[n+1]|=1<<i
  else:out[i]|=1<<(n+1)
 out[n]|=1<<(n+1)
 return out
def hp(out):
 n=len(out);full=(1<<n)-1;dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if not dp[mask][i]:continue
   nxt=out[i]&~mask&full
   while nxt:
    bit=nxt&-nxt;nxt-=bit;j=bit.bit_length()-1
    dp[mask|bit][j]+=dp[mask][i]
 return sum(dp[-1])
def strong(out):
 n=len(out)
 for start in range(n):
  seen=1<<start;todo=seen
  while todo:
   bit=todo&-todo;todo-=bit;i=bit.bit_length()-1
   new=out[i]&~seen;seen|=new;todo|=new
  if seen!=(1<<n)-1:return False
 return True
# Exhaust every legal parameter which could give a target <=X using alpha bound.
X=3000; support=set();enumerated=0
for c,(h,q,A,G,T,S,E) in enumerate(states):
 for p in range((X-h-A)//q+1):
  if p==0 and A==0:continue
  a0=q*p+A;b=0
  while (1<<b)*a0<=X-h:
   for r in range(min((1<<b)-1,(X-h-(1<<b)*a0)//h)+1):
    val=value(c,p,r,b);enumerated+=1
    if val<=X:support.add(val)
   b+=1
raw=subprocess.check_output([str(D/'inverse'),str(D/'CORE_STATES.txt')]+list(map(str,range(1,X+1))),text=True)
rows=[json.loads(l) for l in raw.splitlines()]
assert {d['N'] for d in rows if d['status']=='FOUND'}==support
literal=[]
for d in [z for z in rows if z['status']=='FOUND'][::20][:60]:
 out=graph(d['core'],d['P'],d['R'],d['b'])
 if len(out)>15:continue
 assert hp(out)==d['N'];assert strong(out)
 literal.append({'N':d['N'],'outgoing_masks':out})
for d in rows:
 if d['status']=='FOUND':assert value(d['core'],d['P'],d['R'],d['b'])==d['N']
large=[json.loads(l) for l in (D/'INVERSE_CASES.jsonl').read_text().splitlines()]
largegraphs=[]
for d in large:
 if d['status']=='FOUND':
  assert value(d['core'],d['P'],d['R'],d['b'])==d['N']
  out=graph(d['core'],d['P'],d['R'],d['b']);assert strong(out)
  largegraphs.append(dict(d,outgoing_masks=out))
(D/'LITERAL_GRAPHS.json').write_text(json.dumps(literal,indent=2)+'\n')
(D/'CASE_GRAPHS.json').write_text(json.dumps(largegraphs,indent=2)+'\n')
(D/'COMBINED_VERIFY.json').write_text(json.dumps({'status':'PASS','all_targets_checked':[1,X],'independent_exhaustive_parameter_choices':enumerated,'represented_targets':len(support),'small_literal_subset_dp_graphs':len(literal),'large_formula_and_strongness_checks':len(largegraphs),'omissions47and63':'ABSENT_FROM_GRAMMAR after complete all-length search','caveat':'Large graph counts use the inherited proved construction formula, not subset DP.'},indent=2)+'\n')
print((D/'COMBINED_VERIFY.json').read_text())
