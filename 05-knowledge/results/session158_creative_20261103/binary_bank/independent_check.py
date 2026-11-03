import itertools
from fractran import compile_word

def oracle(w):
    n=len(w); total=0
    for mask in range(1<<n):
        before=[i for i in range(n) if mask & (1<<i)]
        after=[i for i in range(n) if not mask & (1<<i)]
        total += (not before or w[before[-1]]=='0') and (not after or w[after[0]]=='1')
    return total

def direct(prog,n):
    while True:
        for a,b in prog:
            # no prior assumptions on fraction reduction
            q,r=divmod(n*a,b)
            if r==0:
                n=q;break
        else: return n

def formula(w):
    return 2 + sum(2**(i-j-1) for i in range(len(w)) for j in range(i) if w[j]=='1' and w[i]=='0') + sum(w[i:i+2]=='01' for i in range(len(w)-1))
for L in range(2,9):
    for x in itertools.product('01',repeat=L-2):
        w='1'+''.join(x)+'0';prog,n,final=compile_word(w)
        expected=oracle(w)
        if expected != formula(w): raise RuntimeError(('formula mismatch',w,expected,formula(w)))
        if direct(prog,n) != final*5**expected: raise RuntimeError(('program mismatch',w))
print('PASS: independently interpreted 127 programs against apex-path subset bijection and closed formula')
P='495/14 7/11 13/7 34/39 13/17 19/13 207/38 19/23 29/19 62/87 29/31 370/29 205/74 37/41'.split()
P=[tuple(map(int,s.split('/'))) for s in P]; n=350; count=0
while True:
 for a,b in P:
  v,r=divmod(n*a,b)
  if not r: n=v;count+=1;break
 else: break
if (n,count)!=(72265625,32): raise RuntimeError(('printed instance mismatch',n,count))
print('PASS: hard-coded printed 1010 fractions reach',n,'in',count,'steps')
