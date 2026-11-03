from itertools import combinations, permutations

def make(n,code):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if code>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def ear(a,mask):
 n=len(a); b=a[:]+[mask]
 for i in range(n):
  if not(mask>>i&1):b[i]|=1<<n
 return b

def literal(a):
 n=len(a)
 return sum(all((a[i]>>j)&1 for i,j in zip(p,p[1:])) for p in permutations(range(n)))

def is_strong(a):
 n=len(a)
 for v in range(n):
  seen={v}; stack=[v]
  while stack:
   for w in range(n):
    if a[stack[-1]]>>w&1 and w not in seen:seen.add(w);stack.append(w);break
   else:stack.pop()
  if len(seen)!=n:return False
 return True
base=ear(make(5,759),2)
if literal(base)!=31 or not is_strong(base):raise RuntimeError("base failure")
for pm,h,m,target in [(50,93,32,281),(24,95,126,285)]:
 parent=ear(base,pm); child=ear(ear(parent,m),m)
 H=literal(parent); Hc=literal(child)
 if (H,Hc)!=(h,target) or not is_strong(parent) or not is_strong(child):raise RuntimeError('witness failure')
 print('parent_mask',pm,'H',H,'two_masks',m,m,'H_child',Hc,'strong',True)
