"""Standalone exact controls for the physical append-only binary-child theorem."""
from collections import deque
from pathlib import Path
import hashlib,json
T=(1,0,2,1)
def require(ok,message):
 if not ok:raise RuntimeError(message)
def step(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D
 raise ValueError(c)
def state(w):
 s=T
 for c in w:s=step(s,c)
 return s
def F(s):return s[0]-s[1]+s[2]
def h(s):return s[2]-s[0]+s[1]-2*s[3]
def pos(s):
 D,E,U,d=s;return D-E-d,U-D+E-d,E,d-E
def pos_step(y,c):
 a,f,e,p=y
 if c=='O':return a+e,a+2*f+e+2*p,0,e+p
 if c=='Z':return 2*a+f+e+p,f,e+p,0
 return f+p,0,a+e+p,e

def bank(cap):
 seen={T:''};queue=deque([T])
 while queue:
  s=queue.popleft()
  for c in'OZJ':
   q=step(s,c)
   require(q[2]>s[2],'strict source cover growth')
   if q[2]<=cap and q not in seen:seen[q]=seen[s]+c;queue.append(q)
 return seen

def predicted(cap):
 minus=set();plus=set();end=set()
 m=0
 while state('Z'*m)[2]<=cap:
  minus.add(state('Z'*m));m+=1
 r=0
 while state('O'*r)[2]<=cap:end.add(state('O'*r));r+=1
 if state('Z')[2]<=cap:end.add(state('Z'))
 j=1
 while state('Z'*j+'J')[2]<=cap:end.add(state('Z'*j+'J'));j+=1
 bases=['O','ZO'];j=1
 while state('Z'*j+'JO')[2]<=cap:bases.append('Z'*j+'JO');j+=1
 for w in bases:
  k=0
  while state(w+'Z'*k)[2]<=cap:plus.add(state(w+'Z'*k));k+=1
 return minus,plus,end

def graph(w):
 rows=[2,0];root=0
 for c in w:
  n=len(rows);rows.append(0)
  for v in range(n):
   if c!='O'and v==root:rows[n]|=1<<v
   else:rows[v]|=1<<n
  if c=='J':root=n
 return rows

def counts(rows):
 n=len(rows);full=(1<<n)-1
 def calc(starts,cycles):
  dp=[{}for _ in range(1<<n)]
  for v in starts:dp[1<<v][v]=1
  for mask in range(1,full+1):
   for last,count in dp[mask].items():
    todo=rows[last]&~mask&full
    while todo:
     bit=todo&-todo;todo-=bit;v=bit.bit_length()-1
     d=dp[mask|bit];d[v]=d.get(v,0)+count
  return sum(c for v,c in dp[full].items()if not cycles or rows[v]&1)
 return calc(range(n),False),calc([0],True)

def main():
 cap=5000;B=bank(cap)
 for s in B:
  y=pos(s);require(min(y)>=0,'physical nonnegativity')
  require(F(s)>=s[2]+1,'source cover scalar bound')
  for c in'OZJ':
   q=step(s,c)
   require(pos(q)==pos_step(y,c),'independent coordinate recurrence')
   require(2*F(q)>=3*F(s),'one-letter growth lower bound')
  require(F(step(s,'Z'))==2*F(s)+h(s),'Z binary offset')
  require(F(step(s,'O'))==F(step(s,'J'))==2*F(s)-(s[0]-2*s[1]),'O/J binary offset')
 observed=({s for s in B if h(s)==-1},{s for s in B if h(s)==1},{s for s in B if s[0]-2*s[1]==1})
 expected=predicted(cap)
 for a,b,label in zip(observed,expected,['h_minus','h_plus','endpoint']):
  require(a==b,'full physical slice classification '+label)
 target25={s:w for s,w in B.items()if F(s)==25}
 require(set(target25)=={(5,0,20,1),(15,5,15,5)},'all-height parent25 profiles')
 child_values={F(step(s,c))for s in target25 for c in'OZJ'}
 require(child_values=={45,63},'all parent25 one-letter outputs')
 for n in range(5,100000,2):require(9*n>4*(2*n+1),'two-letter exclusion')
 limit=10000001
 def add_dyadic(out,base,shift):
  k=0
  while base*(1<<k)+shift<=limit:out.add(base*(1<<k)+shift);k+=1
 PLUS=set();MINUS=set()
 for base in[6,10]:add_dyadic(PLUS,base,-1)
 j=1
 while 6*((1<<j)+1)-1<=limit:add_dyadic(PLUS,6*((1<<j)+1),-1);j+=1
 r=1
 while (1<<r)+1<=limit:MINUS.add((1<<r)+1);r+=1
 j=1
 while 3*((1<<j)+1)<=limit:MINUS.add(3*((1<<j)+1));j+=1
 require(PLUS&MINUS=={5,9,17},'binary parent intersection')
 literal=[]
 for n,wm,wp in[(5,'Z','O'),(9,'ZZ','ZO'),(17,'ZZZ','ZJO')]:
  for sign,w in[(-1,wm),(1,wp)]:
   require(F(state(w))==n and h(state(w))==sign,'signed parent source')
   cw=w+'Z';require(F(state(cw))==2*n+sign,'actual binary child')
   rows=graph(cw+'J');H,C=counts(rows)
   require((H,C)==(2*n+sign,1),'literal binary child graph')
   literal.append(dict(parent=n,sign=sign,word=cw+'J',H=H,C=C,vertices=len(rows),rows=rows))
 for n,w in[(25,'OJOO'),(49,'JZJJZ'),(51,'OOZZ')]:
  require(F(state(w))==n,'literal target source')
  rows=graph(w+'J');H,C=counts(rows);require((H,C)==(n,1),'literal hostile target graph')
  literal.append(dict(target=n,word=w+'J',H=H,C=C,vertices=len(rows),rows=rows))
 formulas=0
 for j in range(1,41):
  for k in range(41):
   for word,value in [('O'+'Z'*k,6*(1<<k)-1),('ZO'+'Z'*k,10*(1<<k)-1),('Z'*j+'JO'+'Z'*k,6*((1<<j)+1)*(1<<k)-1)]:
    s=state(word);require(F(s)==value and h(s)==1,'unbounded-family formula sample');formulas+=1
 receipt=dict(status='PASS',physical_cover_cap=cap,physical_profiles=len(B),
  slice_profile_counts=dict(zip(['h_minus','h_plus','endpoint'],map(len,observed))),
  parent25_profiles=[dict(profile=list(s),word=w)for s,w in sorted(target25.items())],
  parent25_one_letter_outputs=sorted(child_values),
  intersection_check_bound=limit,intersection=sorted(PLUS&MINUS),
  formula_checks=formulas,literal_graphs=literal,
  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  scope='Exact finite controls support the all-order physical append-only classification. No eventual coverage assertion.')
 print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
