import importlib.util,itertools,random
spec=importlib.util.spec_from_file_location('bank',str(__import__('pathlib').Path(__file__).resolve().parent / 'bank.py'))
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
def direct(word,ab):
 m=len(word); total=0
 def edge(u,v):
  if u>=m and v>=m:return bool(ab) if (u-m,v-m)==(0,1) else not bool(ab)
  if u<m and v<m:return u<v
  if u<m:return bool(word[u]>>(v-m)&1)
  return not bool(word[v]>>(u-m)&1)
 for order in [(m,m+1),(m+1,m)]:
  for gaps in itertools.product(range(3),repeat=m):
   ls=[]
   for j in range(3):
    ls.extend(i for i in range(m) if gaps[i]==j)
    if j<2:ls.append(order[j])
   total+=all(edge(u,v) for u,v in zip(ls,ls[1:]))
 return total
n=0
for m in range(6):
 for word in itertools.product(range(4),repeat=m):
  for ab in (0,1):
   result=direct(word,ab)
   if result!=b.automaton(word,ab): raise AssertionError((word,ab,result))
   n+=1
print('direct assignment agreement',n)
print('counterexample type word (1,0,1), 0 beats 1: strong,H', b.strong(b.rows((1,0,1),1)),direct((1,0,1),1))
# Random unrestricted tournament with arbitrary chosen feedback deletion S:
# compare claimed bound with direct endpoint dynamic programming.
random.seed(156); cases=0
for k in range(1,5):
 for m in range(1,9):
  for z in range(100):
   n=k+m; A=[0]*n
   def add(i,j):A[i]|=1<<j
   for i in range(m):
    for j in range(i+1,m):add(i,j)
   for i in range(m,n):
    for j in range(i+1,n): add(i,j) if random.randrange(2) else add(j,i)
   typ=[]
   for i in range(m):
    mask=0
    for a in range(k):
     if random.randrange(2):add(i,m+a);mask|=1<<a
     else:add(m+a,i)
    typ.append(mask)
   if b.strong(A):
    t=max(typ.count(i) for i in set(typ)); h=b.dp(A)
    if h<2**((t+1)//2)+1: raise AssertionError((k,m,typ,h,t))
    cases+=1
print('random strong independent DP bound checks',cases)
