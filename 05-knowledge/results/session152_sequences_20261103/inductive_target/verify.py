from pathlib import Path
exec((Path(__file__).resolve().parent.parent / 'cut_realizability' / 'check.py').read_text().split('for n in range(2,5):')[0])
def ear(a,m):
 n=len(a); return [r+[1-(m>>i&1)] for i,r in enumerate(a)]+[[(m>>i)&1 for i in range(n)]+[0]]
def permhp(a):
 n=len(a);return sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
def strong(a):
 n=len(a)
 for i in range(n):
  seen={i}; todo=[i]
  while todo:
   u=todo.pop()
   for j in range(n):
    if a[u][j] and j not in seen:seen.add(j);todo.append(j)
  if len(seen)!=n:return False
 return True
for code,masks in [(759,[2]),(759,[2,50]),(759,[2,24]),(759,[2,28]),(759,[2,28,2]),(759,[2,28,78]),(8,[10,19]),(8,[10,33])]:
 a=adjacency(5,code)
 for m in masks:a=ear(a,m)
 h=hp(a);p=permhp(a)
 if h!=p or not strong(a): raise AssertionError((code,masks,h,p))
 print(code,masks,'H',h,'DP=permutation; strong')
