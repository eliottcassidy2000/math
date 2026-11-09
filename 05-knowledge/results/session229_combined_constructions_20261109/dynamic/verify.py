"""Portable standard-library independent routes, explicit checks (no asserts)."""
from pathlib import Path
import importlib.util
_spec = importlib.util.spec_from_file_location("dynamic_core_probe", Path(__file__).resolve().with_name("probe.py"))
_probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_probe)
states, compose = _probe.states, _probe.compose
import itertools as it,json

def need(x,msg):
 if not x: raise RuntimeError(msg)
def f(p):
 bits=list(map(int,bin(p)[2:]));return sum(2**(j-i-1)+(j==i+1) for i in range(len(bits)) for j in range(i+1,len(bits)) if bits[i] and bits[j])
def value(st,P,R,b):
 h,q,A,G,T,S,E=st
 rev=int(format(R,'0%db'%b)[::-1],2) if b else 0
 return (q*P+A)*(2**(b+1)-rev)+h*(2*R+1-f(R))-P*G-h*f(P)-T-(P%2)*S-E*(R>>(b-1) if b else 0)
def literal(B,mask,P,R,b):
 p=list(map(int,bin(P)[2:])) if P else [];r=list(map(int,format(R,'0%db'%b))) if b else []
 blocks=[[0]]*len(p)+[B]+[[0]]*b;masks=p+[mask]+r
 n=sum(map(len,blocks));out=[0]*(n+2);o=0
 for C,m in zip(blocks,masks):
  for i,c in enumerate(C):
   out[o+i]=(c<<o)|sum(1<<j for j in range(o+len(C),n+1))
   if m>>i&1:out[n+1]|=1<<(o+i)
   else:out[o+i]|=1<<(n+1)
  o+=len(C)
 out[n]|=1<<(n+1)
 return out

def hp(B):
 n=len(B);d=[{} for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for m in range(1,1<<n):
  for v,c in d[m].items():
   rem=B[v]&~m
   while rem:
    bit=rem&-rem;rem-=bit;v2=bit.bit_length()-1;dest=d[m|bit];dest[v2]=dest.get(v2,0)+c
 return sum(d[-1].values())
def permutation_state(B,m):
 n=len(B);h=q=A=G=T=S=E=0
 for p in it.permutations(range(n)):
  if all(B[u]>>v&1 for u,v in zip(p,p[1:])):
   h+=1;S+=m>>p[0]&1;E+=m>>p[-1]&1
  for k in range(n+1):
   x,y=p[:k],p[k:]
   if all(B[u]>>v&1 for z in (x,y) for u,v in zip(z,z[1:])):
    q+=1;a=bool(y and m>>y[0]&1);g=bool(x and m>>x[-1]&1);A+=a;G+=g;T+=a*g
 return h,q,A,G,T,S,E
B=[24,17,3,6,12];st=states(B)
need(st[3]==(15,120,42,42,10,6,6),'first core')
need(st[5]==(15,120,42,42,8,6,6),'second core')
for m in range(32):need(st[m]==permutation_state(B,m),'independent permutation state')
checks=0
for a in range(4):
 for P in (range(2**(a-1),2**a) if a else [0]):
  for b in range(4-a):
   for R in range(2**b):
    for mask in (3,5):need(hp(literal(B,mask,P,R,b))==value(st[mask],P,R,b),'literal closure');checks+=1
    need(value(st[5],P,R,b)-value(st[3],P,R,b)==2,'universal difference')
# Universal triple, checked by literal cover splits and whole graphs.
triples=[([48,49,3,7,12,28],9,20),([48,49,3,7,12,28],3,22),([40,37,17,6,11,28],3,24)]
def strong(C):
 for v in range(len(C)):
  reached=1<<v
  while True:
   grown=reached
   for u in range(len(C)):
    if reached>>u&1:grown|=C[u]
   if grown==reached:break
   reached=grown
  if reached!=(1<<len(C))-1:return False
 return True
need(strong(B),'B5 strong')
triple_checks=0
for C,m,T in triples:
 expected=(45,432,129,129,T,15,15)
 need(states(C)[m]==expected and permutation_state(C,m)==expected,'triple state independent')
 need(strong(C),'triple core strong')
 for P,R,b in [(0,0,0),(0,0,2),(0,1,1),(1,0,0),(1,1,1),(2,2,2),(3,0,1)]:
  graph=literal(C,m,P,R,b)
  need(strong(graph) and hp(graph)==value(expected,P,R,b),'literal triple closure')
  triple_checks+=1
# Complete small universe: all labeled tournaments and all cuts, n<=4.
minimal=[]
for n in range(1,5):
 bank={};cnt=0
 for z in range(1<<(n*(n-1)//2)):
  C=[0]*n
  for k,(i,j) in enumerate(it.combinations(range(n),2)):
   u,v=(i,j) if z>>k&1 else (j,i);C[u]|=1<<v
  for h,q,A,G,T,S,E in states(C):
   key=h,q,A,G,S,E
   need(all(abs(T-t)!=2 for t in bank.get(key,[])),'smaller universal switch')
   bank.setdefault(key,set()).add(T);cnt+=1
 minimal.append([n,cnt,len(bank)])
H=[40,1,3,6,15,30];ss=states(H);fib={}
for m,s in enumerate(ss):fib.setdefault(s[2],[]).append(m)
diffs=[]
for A,ms in sorted(fib.items()):
 if len(ms)==2:
  x,y=ms;delta=tuple(v-u for u,v in zip(ss[x],ss[y]));diffs.append([ms,list(delta)])
  for P in range(8):
   for b in range(4):
    for R in range(2**b):
     e=R>>(b-1) if b else 0
     need(value(ss[y],P,R,b)-value(ss[x],P,R,b)==-P*delta[3]-delta[4]-(P%2)*delta[5]-e*delta[6],'H19 delta')
need(hp(literal(H,19,0,1,1))==345 and hp(literal(H,42,0,1,1))==343,'H19 minimal twin')
C=states([16,1,3,7,14]);u=compose(C[2],C[28]);v=compose(C[4],C[0]);delta=tuple(b-a for a,b in zip(u,v))
need(delta==(0,0,0,-3021,-558,9,-36),'carry difference')
# Vertex-transitive larger cyclic cores are hostile to a general interval inference.
cyclic=[]
for n in (5,7,9):
 C=[sum(1<<((u+j)%n) for j in range(1,(n+1)//2)) for u in range(n)];s=states(C)
 fibers=[]
 for k in range(2,n//2+1):
  entries=[x for m,x in enumerate(s) if m.bit_count()==k]
  need(len({(h,q,A,G,S,E) for h,q,A,G,T,S,E in entries})==1,'transitive core coordinates')
  fibers.append([k,sorted({x[4] for x in entries})])
 cyclic.append([n,list(s[0][:2]),fibers])
print(json.dumps({'status':'PASS','literal_closures':checks,'literal_triple_closures':triple_checks,'minimality_universes':minimal,'H19_fibers':diffs,'carry_difference':delta,'cyclic_hostiles':cyclic},indent=2))
