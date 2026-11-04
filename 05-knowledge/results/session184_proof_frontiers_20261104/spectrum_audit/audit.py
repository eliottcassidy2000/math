from pathlib import Path
from itertools import combinations
from math import prod
import json
BASE=Path('/workspace/shared/session184/spectrum')
def need(v,m):
 if not v:raise RuntimeError(m)
def blocks(w):
 return [(1,2,a,a,0,a,a) if k==1 else (3,12,3*a,3*a,a*(a-1),a,a) for k,a in w]
def sig(w):
 b=blocks(w);n=len(b)
 h=prod(x[0] for x in b);p=prod(x[1] for x in b)
 f=sum(b[i][2]*prod(b[j][1] if j<i else b[j][0] for j in range(n) if j!=i) for i in range(n))
 g=sum(b[i][3]*prod(b[j][1] if j>i else b[j][0] for j in range(n) if j!=i) for i in range(n))
 J=sum(b[i][4]*prod(b[j][0] for j in range(n) if j!=i) for i in range(n))
 J+=sum(b[i][5]*b[i+1][6]*prod(b[j][0] for j in range(n) if j not in (i,i+1)) for i in range(n-1))
 J+=sum(b[i][3]*b[j][2]*prod(b[l][1] if i<l<j else b[l][0] for l in range(n) if l not in (i,j)) for i in range(n) for j in range(i+1,n))
 return (h,p,f,g,J,prod(x[0] for x in b[:-1])*b[-1][5],b[0][6]*prod(x[0] for x in b[1:]))
def inv(F):
 out=set();k=0
 while F%3**k==0:
  u=F//3**k; L=u.bit_length()
  # Independent exhaustive choice of disjoint two-bit intervals.
  for width in (L,L+1):
   for starts in combinations(range(width-1),k):
    if any(b==a+1 for a,b in zip(starts,starts[1:])):continue
    w=[];pos=0
    while pos<width:
     length=2 if pos in starts else 1
     w.append((3 if length==2 else 1,(u>>pos)&((1<<length)-1)))
     pos+=length
    if w[-1][1] and sig(w)[2]==F:out.add(tuple(w))
  k+=1
 return out
def candidates(odd):return {w for F in (odd,2*odd) for w in inv(F) if w[0][1]}
counts={f:len(candidates(f)) for f in (21,63,189,567,1701)}
need(list(counts.values())==[5,10,18,29,50],('counts',counts))
queue=[('',((1,1),)*3)];nodes={};all_candidates=set();tests=0
while queue:
 a,w=queue.pop();z=sig(w);c=z[0]-z[4];odd=3*z[2]
 while odd%2==0:odd//=2
 cs=candidates(odd);all_candidates|=cs
 children={str(d):[] for d in (2,4)}
 for d in (2,4):
  for v in cs:
   tests+=1;s=sig(v)
   if s[0]-s[4]==3*c+d:
    children[str(d)].append(a+str(d));queue.append((a+str(d),v))
 nodes[a]={'word':w,'signature':z,'c':c,'children':children}
need(set(nodes)=={'','2','4','22','42','422','4222'},'language')
r=json.loads((BASE/'receipt.json').read_text())
for node in r['return_nodes']:
 own=nodes[node['address']]
 need(tuple(map(tuple,node['word']))==own['word'],'node word')
 need(tuple(node['signature'])==own['signature'],'node signature')
 need(node['c']==own['c'] and node['children']==own['children'],'node transitions')
need(tests==280 and len(all_candidates)==112,'coverage')
cert=json.loads((BASE/'candidates.json').read_text());need(len(cert)==112,'certificate length')
need({tuple(map(tuple,x['word'])) for x in cert}==all_candidates,'candidate set')
for x in cert:
 s=sig(tuple(map(tuple,x['word'])))
 need(tuple(x['signature'])==s and x['c']==s[0]-s[4],'candidate signature')
# Arbitrary zero-port prefixes, including triangle prefixes, check normalization law.
normal=0
for w in all_candidates:
 z=sig(w);c=z[0]-z[4]
 for m in range(8):
  ss=sig(((1,0),)*m+w)
  need(ss[2]==2**m*z[2] and ss[0]-ss[4]==c,'singleton normalization')
  tt=sig(((3,0),)+((1,0),)*m+w)
  need(tt[0]-tt[4]==3*c,'triangle obstruction')
  normal+=2
# Leading coefficient on mixed singleton/triangle patterns, different mark counts.
valuation=0
for k in range(2,15):
 for a in range(4):
  for b in range(2):
   for separated in (False,True):
    w=[]
    for j in range(k):
     w.append((3,(a+j)%4))
     if separated and j<k-1:w.append((1,b))
    z=sig(w);A=sum(x[1]*y[1] for x,y in zip(w,w[1:]) if x[0]==y[0]==3)
    for t in (1,2,9):
     H=z[0]-z[4]+2**t*z[2]
     need(H%3**(k-2)==0 and (H//3**(k-2)+A)%3==0,'valuation coefficient')
     valuation+=1
print(json.dumps({'status':'PASS','candidate_counts':counts,'return_tests':tests,'unique_candidate_words':len(all_candidates),'normalization_checks':normal,'leading_coefficient_checks':valuation,'language':sorted(nodes)},sort_keys=True,indent=2))
