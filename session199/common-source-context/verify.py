#!/usr/bin/env python3
"""Exact physical common-context verifier. Standard library only."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def need(p,msg):
 if not p:raise RuntimeError(msg)
def move(v,c):
 D,E,U,d=v
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D
 raise ValueError(c)
def apply(v,w):
 for c in w:v=move(v,c)
 return v
T=(1,0,2,1)
def sig(v):
 D,E,U,d=v;return D+2*d,E+d,U-d
def decode_sig(s):
 A,B,C=s
 if B%2:d=B if A>=3*B else A-B-C
 else:d=B//2 if 2*A<=3*B+2*C else A-B-C
 return A-2*d,B-d,C+d,d

def n_raw(v,w):
 D,E,U,d=v;F,G,V,e=w
 return (D-E+U)*F-(D-E+d)*G+(D+3*U-d)*V-(U+E)*e

def digits(v):
 A,B,C=sig(v);return A-2*B,6*(A-B)+8*C,4*(A+3*C)
def encoded(v,x):
 q0,q1,q2=digits(v);return q0+x*q1+x*x*q2

def decoded(n,x):
 q0=n%x;q1=(n//x)%x;q2=n//(x*x)
 need(q2%4==0,'leading digit multiple of four');need(q1%2==0,'middle digit even')
 L=q2//4;S=q1//2
 C=3*L+3*q0-2*S;A=-8*L-9*q0+6*S;B=-4*L-5*q0+3*S
 return decode_sig((A,B,C))

def bank(cap):
 out={T:''};queue=[T]
 for s in queue:
  for c in 'OZJ':
   t=move(s,c)
   if t[2]<=cap and t not in out:out[t]=out[s]+c;queue.append(t)
 return out

def graph(word):
 rows=[2,0];root=0
 for c in word:
  n=len(rows);new=0
  for v in range(n):
   if c!='O' and v==root:new|=1<<v
   else:rows[v]|=1<<n
  rows.append(new)
  if c=='J':root=n
 return rows,root

def induced(rows,vertices):
 return [sum(((rows[v]>>u)&1)<<j for j,u in enumerate(vertices))for v in vertices]
def converse(rows):
 mask=(1<<len(rows))-1
 return [mask^r^(1<<v)for v,r in enumerate(rows)]
def graft_N(vword,wword):
 left,c=graph(vword);n=len(left)
 # Birth a beating only c.
 for v in range(n):
  if v!=c:left[v]|=1<<n
 left.append(1<<c);entries={c,n};n+=1
 right,r=graph(wword);right=converse(right)
 rows=left+[s<<n for s in right]
 for v in range(n):
  for u in range(len(right)):
   if u==r and v in entries:rows[n+u]|=1<<v
   else:rows[v]|=1<<(n+u)
 return rows

def paths_all(rows,root=None):
 n=len(rows);N=1<<n;dp=[None]*N;total=[0]*N;total[0]=1
 for v in range(n):
  if root is None or v==root:
   dp[1<<v]={v:1};total[1<<v]=1
 for mask in range(1,N):
  cur=dp[mask]
  if cur is None:continue
  total[mask]=sum(cur.values())
  for v,num in cur.items():
   avail=rows[v]&(N-1-mask)
   while avail:
    bit=avail&-avail;avail-=bit;u=bit.bit_length()-1;m=mask|bit
    if dp[m]is None:dp[m]={}
    dp[m][u]=dp[m].get(u,0)+num
 return total,dp[-1]or{}

def direct_profile(rows,root):
 n=len(rows);full=(1<<n)-1;tot,ends=paths_all(rows);starts,_=paths_all(rows,root)
 R=sum(starts[m]*tot[full^m]for m in range(1,full+1)if m>>root&1)
 return (tot[-1],ends.get(root,0),R,tot[full^(1<<root)]),starts[-1]
def direct_HC(rows):
 total,_=paths_all(rows);_,ends=paths_all(rows,0)
 return total[-1],sum(h for v,h in ends.items()if rows[v]&1)

def run():
 result={}
 # Symbolic identities reduced to exact linearity: a basis independently
 # checks the original four-coordinate form, rather than reusing signature form.
 for k in range(1,11):
  x=2**k;w=apply(T,'O'*k+'JJ'+'O'*k)
  need(w==(1+4*x,0,2*x+4*x*x,1+2*x),'partner response')
  need(sig(w)==(3+8*x,1+2*x,4*x*x-1),'partner signature')
  for j in range(4):
   v=tuple(int(i==j)for i in range(4))
   need(n_raw(v,w)==encoded(v,x),'linear identity on basis')
 result['basis_cases']=40
 # Finite bounded bank. Monotone U makes this exact; no order cutoff.
 cap=5000;b=bank(cap);k=(14*cap-1).bit_length();x=1<<k
 values={};three=set()
 probes=[T,move(T,'J'),apply(T,'JJ')]
 for v,word in b.items():
  q0,q1,q2=digits(v);D,E,U,d=v
  need(1<=q0<=cap,'low digit bound');need(0<=q1<=14*cap-2,'middle digit bound')
  need(0<q2<=16*cap-4,'leading bound');need(q0<x and q1<x,'no carries')
  need(decode_sig(sig(v))==v,'physical signature inverse')
  n=encoded(v,x);need(decoded(n,x)==v,'single-context inverse')
  need(n not in values,'single-context injectivity');values[n]=v
  r0,r1,r2=[n_raw(v,w)for w in probes]
  A=3*r0+2*r1-2*r2;B=(-r0+6*r1-3*r2)//2;C=(-3*r0+r2)//2
  need((A,B,C)==sig(v),'three-small-context inverse');three.add((r0,r1,r2))
 result['finite_exact']={'U_cap':cap,'profiles':len(b),'x':x,'k':k,'single_context_values':len(values),'three_context_values':len(three),'max_encoded':max(values)}
 # Literal graphs, including every word through length three, and enough
 # partner contexts to test each O/J boundary.
 small={''}
 for _ in range(3):small|={w+c for w in tuple(small)for c in 'OZJ'if len(w)<3}
 profiles=0;graphs=0
 for word in sorted(small,key=lambda w:(len(w),w)):
  rows,r=graph(word);v,s=direct_profile(rows,r)
  need(s==1 and v==apply(T,word),'literal source profile');profiles+=1
  if len(word)<=2:
   for k0 in range(1,3):
    ww='O'*k0+'JJ'+'O'*k0
    h,c=direct_HC(graft_N(word,ww))
    need(c==2 and h==encoded(v,2**k0),'literal C2 common context');graphs+=1
 # Every pair of the three specified small contexts has a physical collision.
 irredundancy=[]
 for inds,lword,rword,expected in [((0,1),'OJ','ZJ',(41,73)),((0,2),'OOJJOOZ','ZOJOZOZ',(1041,3375)),((1,2),'OJOZ','ZJOO',(277,501))]:
  lv=apply(T,lword);rv=apply(T,rword)
  need(lv!=rv,'irredundancy distinct profiles')
  need(tuple(n_raw(lv,probes[i])for i in inds)==expected,'left irredundancy values')
  need(tuple(n_raw(rv,probes[i])for i in inds)==expected,'right irredundancy values')
  for word in (lword,rword):
   for i,n in zip(inds,expected):
    need(direct_HC(graft_N(word,'J'*i))==(n,2),'literal irredundancy graph')
  omitted=next(i for i in range(3)if i not in inds)
  need(n_raw(lv,probes[omitted])!=n_raw(rv,probes[omitted]),'omitted context separates')
  irredundancy.append({'probes':inds,'left_word':lword,'left':lv,'right_word':rword,'right':rv,'common_values':expected})
 result['three_probe_irredundancy']=irredundancy
 # Full no-carry context at the smallest nontrivial bound.
 fullcap=2;k0=(14*fullcap-1).bit_length();xx=2**k0
 rows=graft_N('','O'*k0+'JJ'+'O'*k0);h,c=direct_HC(rows)
 need((h,c)==(encoded(T,xx),2),'full bounded selector graph')
 result['literal_checks']={'source_words':profiles,'small_C2_graphs':graphs,'full_selector_order':len(rows),'full_selector_H':h,'full_selector_C':c,'full_selector_bound':fullcap,'full_selector_x':xx}
 # Hostile: a fixed undersized probe really collapses different sources.
 a=apply(T,'OZJ');b0=apply(T,'ZOZ');need(a!=b0,'hostile distinct')
 hostile_w=apply(T,'O'+'JJ'+'O')
 need(n_raw(a,hostile_w)==n_raw(b0,hostile_w)==931,'hostile collision')
 for word in ('OZJ','ZOZ'):
  need(direct_HC(graft_N(word,'OJJO'))==(931,2),'literal hostile collision')
 result['undersized_probe_hostile']={'x':2,'left':a,'left_word':'OZJ','right':b0,'right_word':'ZOZ','common_H':931}
 # Formal cone points remain invalid despite entirely valid digit syntax.
 fake=(1,0,3,1);need(fake not in bank(3),'fake source rejected')
 need(digits(fake)==(1,28,36),'formal digit hostile')
 result['nonphysical_digit_hostile']={'profile':fake,'digits':digits(fake),'reason':'D=1 forces a transitive source and U a power of two'}
 result['dependency_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((HERE/'dependencies').glob('*.md'))}
 print(json.dumps(result,indent=2,sort_keys=True))
 return result
if __name__=='__main__':run()
