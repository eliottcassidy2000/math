"""Independent exact fixed-block controls. No author code imports."""
import itertools,functools,json,math,random

def need(b,m):
 if not b:raise RuntimeError(m)
def path_counter(rows):
 @functools.cache
 def suffix(last,left):
  if not left:return 1
  todo=left&rows[last];v=0
  while todo:
   bit=todo&-todo;todo-=bit;v+=suffix(bit.bit_length()-1,left^bit)
  return v
 def H(S=None):
  if S is None:S=(1<<len(rows))-1
  if S==0:return 1
  todo=S;v=0
  while todo:
   bit=todo&-todo;todo-=bit;v+=suffix(bit.bit_length()-1,S^bit)
  return v
 return H

def sccs(rows):
 # Tarjan, independent of the author's transitive-closure SCC check.
 n=len(rows);idx=[-1]*n;low=[-1]*n;stack=[];on=set();out=[];timer=0
 def visit(v):
  nonlocal timer
  idx[v]=low[v]=timer;timer+=1;stack.append(v);on.add(v)
  todo=rows[v]
  while todo:
   bit=todo&-todo;todo-=bit;w=bit.bit_length()-1
   if idx[w]<0:visit(w);low[v]=min(low[v],low[w])
   elif w in on:low[v]=min(low[v],idx[w])
  if low[v]==idx[v]:
   C=0
   while True:
    w=stack.pop();on.remove(w);C|=1<<w
    if w==v:break
   out.append(C)
 for v in range(n):
  if idx[v]<0:visit(v)
 return out

def glue(block,m,pat):
 s=len(block);n=1+s*m;rows=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if i==0:winner=0 if pat>>(j-1)&1 else j
   else:
    a,b=(i-1)//s,(j-1)//s
    winner=i if a<b or (a==b and block[(i-1)%s]>>((j-1)%s)&1) else j
   loser=j if winner==i else i;rows[winner]|=1<<loser
 return tuple(rows)

def bin_covers(block):
 s=len(block);out=[]
 for perm in itertools.permutations(range(s)):
  for cut in range(s+1):
   A,B=perm[:cut],perm[cut:]
   if all(block[a]>>b&1 for seq in [A,B] for a,b in zip(seq,seq[1:])):out.append((A,B))
 return out

def literal_injection(block,m,pat,rows):
 if m<2:return 0
 s=len(block);paths=[p for p in itertools.permutations(range(s)) if all(block[a]>>b&1 for a,b in zip(p,p[1:]))]
 first=next(p for p in paths if pat>>p[0]&1)
 last=next(p for p in paths if not(pat>>(s*(m-1)+p[-1])&1))
 first=tuple(1+a for a in first);last=tuple(1+s*(m-1)+a for a in last)
 covers=bin_covers(block);found=set()
 for choices in itertools.product(covers,repeat=m-2):
  A=tuple(1+s*(i+1)+v for i,(a,b) in enumerate(choices) for v in a)
  B=tuple(1+s*(i+1)+v for i,(a,b) in enumerate(choices) for v in b)
  path=A+last+(0,)+first+B
  need(len(path)==len(rows) and len(set(path))==len(rows),'injection vertices')
  need(all(rows[a]>>b&1 for a,b in zip(path,path[1:])),'injection arcs')
  found.add(path)
 need(len(found)==len(covers)**(m-2),'injection count')
 return len(found)

summaries=[];total=0;injection_paths=0
for name,block,h0,beta0,limit in [('C3',(2,4,1),3,12,3),('D4',(8,1,3,6),5,32,2)]:
 s=len(block);H=path_counter(block);full=(1<<s)-1;h=H();beta=sum(H(A)*H(full^A) for A in range(1<<s))
 need((h,beta)==(h0,beta0),'core constants')
 need(beta==len(bin_covers(block)),'permutation beta check')
 need(len(sccs(block))==1 and h>=s,'strong block')
 need(beta>=2**s+2*(h-1)>2**s,'strict beta bound')
 universes=[(m,list(range(1<<(s*m))),'complete') for m in range(limit+1)]
 if name=='D4':
  rng=random.Random(20261106);pats=sorted({0,(1<<12)-1,1,1<<11,0x555,0xaaa}|{rng.randrange(1<<12) for _ in range(32)})
  universes.append((3,pats,'declared deterministic sample'))
 for m,pats,scope in universes:
  strong=0;low=None;hi=0;injections=0
  for pat in pats:
   rows=glue(block,m,pat);count=path_counter(rows);val=count();parts=sccs(rows);root=next(C for C in parts if C&1)
   indices={i for i in range(m) if root&(full<<(1+s*i))}
   ell=len(indices)
   need(all((root&(full<<(1+s*i)))==(full<<(1+s*i)) for i in indices),'whole root blocks')
   need(not indices or indices==set(range(min(indices),max(indices)+1)),'root interval')
   for C in parts:
    if C!=root:need(any(C==(full<<(1+s*i)) for i in range(m)),'nonroot SCC exactly block')
   need(val==h**(m-ell)*count(root),'root factorization')
   need(val==math.prod(count(C) for C in parts),'all SCC product')
   if ell>=2:need(count(root)>=beta**(ell-2),'two-bin floor')
   if len(parts)==1:
    strong+=1;low=val if low is None else min(low,val)
    if m>=2:injections+=literal_injection(block,m,pat,rows)
   hi=max(hi,val);total+=1
  injection_paths+=injections
  summaries.append({'block':name,'blocks':m,'scope':scope,'presentations':len(pats),'strong':strong,'minimum_strong_H':low,'maximum_H':hi,'literal_injected_paths':injections})
# The excluded singleton is exactly critical.
H=path_counter((0,));singleton_beta=sum(H(A)*H(1^A) for A in [0,1]);need(singleton_beta==2,'singleton critical')
# Independent literal compilation of prime-atom seeds.
def compile_word(w):
 rows=[2,0];r=0
 for c in w:
  n=len(rows);rows.append(0)
  for i in range(n):
   if c!='O' and i==r:rows[n]|=1<<i
   else:rows[i]|=1<<n
  if c=='J':r=n
 return tuple(rows)
seeds=[]
for n,w in [(49,'OJJOOJ'),(63,'OZZOOJ'),(343,'OOJOZJJOJ')]:
 rows=compile_word(w);H=path_counter(rows)();need(H==n,'literal prime seed')
 seeds.append({'word':w,'H':H,'vertices':len(rows)})
need(total==896,'presentation count')
print(json.dumps({'status':'PASS','scope':'Finite controls only; all-order proof audited separately','presentations':total,'literal_injection_paths':injection_paths,'universes':summaries,'singleton_beta':singleton_beta,'literal_prime_atom_seeds':seeds},indent=2,sort_keys=True))
