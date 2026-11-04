"""Independent permutation audit: no subset DP used in perm_response."""
import itertools,json
from pathlib import Path
from search import response,strong
def birth(rows,ins):
 n=len(rows);return [r|((1<<n) if ins>>i&1 else 0) for i,r in enumerate(rows)]+[((1<<n)-1)^ins]

def check(condition):
 if not condition: raise RuntimeError("verification failed")

def perm_response(rows,w):
 n=len(rows);h=0;f=0
 for p in itertools.permutations(range(n)):
  bad=[i for i in range(n-1) if not (rows[p[i]]>>p[i+1]&1)]
  h+=not bad
  for k in range(1,n+1):
   if p[k-1]==w and all(i==k-1 for i in bad):f+=1
 return h,f

def attach(B,w,s):
 # B labels0..n-1, prefix n..n+s-1, root n+s.
 n=len(B);r=n+s;rows=B[:]+[0]*(s+1)
 for i in range(s):rows[n+i]=((1<<n)-1)|sum(1<<(n+j) for j in range(i+1,s))
 rows[w]|=1<<r;rows[r]=((1<<r)-1)^(1<<w)
 return rows

def main():
 B=[8,1,3,6];X=birth(B,13)
 variants={2:birth(X,15),4:birth(X,1)}
 receipt={'suffixes':[],'clock_tests':[]}
 for d,Y in variants.items():
  check(strong(Y))
  check(Y[:4] != B) # New-bit changes only; old induced relation preserved.
  check([a&15 for a in Y[:4]]==B)
  hp,fp=perm_response(Y,0);hd,fd=response(Y)
  check((hp,fp)==(hd,fd[0])==(15+d,21))
  receipt['suffixes'].append({'d':d,'rows':Y,'port':0,'H':hp,'E2':fp})
  for s in range(1,5):
   source=attach(B,0,s);target=attach(Y,0,s)
   hs,fs=response(source);ht,ft=response(target)
   check(strong(source) and strong(target))
   check(hs==5+7*2**s and ht==3*hs+d)
   check(fs[-1]==hs and ft[-1]==ht)
   if s<=2:check(perm_response(target,len(target)-1)==(ht,ht))
   receipt['clock_tests'].append({'d':d,'s':s,'source_H':hs,'target_H':ht,'source_root_E2':fs[-1],'target_root_E2':ft[-1]})
 Path(__file__).with_name('verification.json').write_text(json.dumps(receipt,indent=2))
 print('PASS: exact suffix counts, actual induced embeddings, both branches, four clocks, independent permutations at first two clocks')
if __name__=='__main__':main()
