"""Verify period25 primitivity with six pairs, without enumerating word states."""
from pathlib import Path
import json

def require(ok,message):
    if not ok: raise RuntimeError(message)

def prime(p):
    return p>=2 and all(p%d for d in range(2,int(p**0.5)+1))

def evaluate(w,n,r,N):
    require(0<=w<1<<n,'word width')
    require(not(w & (w>>1)),'adjacent ones')
    require(not((w&1) and (w>>(n-1)&1)),'cyclic adjacent ones')
    return sum(pow(r,i,N) for i in range(n) if w>>i&1)%N

def reverse(w,n):
    return sum(((w>>i)&1)<<((-i)%n) for i in range(n))

n=25;N=167761;r=135722
require(11*101*151==N and all(prime(p) for p in (11,101,151)),'factorization')
require((r*r-r-1)%N==0 and pow(r,n,N)==1,'chosen root')
require(all(pow(r,d,N)!=1 for d in (1,5)),'exact full clock')
a,b=0,1
for _ in range(n-1):a,b=b,a+b
# a=F24,b=F25; alpha=(F24-1)+F25 phi.
require((a-1+b*r)%N==0,'period ideal')
require(abs((a-1)**2+(a-1)*b-b*b)==N,'ideal norm')
cert=json.loads(Path(__file__).with_name('sparse_certificate.json').read_text())
proper={11,101,151,11*101,11*151,101*151}
require({c['d'] for c in cert}==proper,'complete proper nontrivial divisor list')
for c in cert:
    d=c['d']
    for label in ('x','y'):
        w=c['word_'+label]
        require(evaluate(w,n,r,N)==c[label],'source reader')
        require(evaluate(reverse(w,n),n,r,N)==c['R'+label],'reflected reader')
    require((c['x']-c['y'])%d==0,'same candidate block')
    require((c['Rx']-c['Ry'])%d!=0,'reflection separates candidate block')
print(json.dumps({'period':n,'states':N,'source_word_pairs_checked':len(cert),'word_enumeration':False,'conclusion':'<add1,word_reflection> primitive; every nonconstant CRT observation has full future carrier'},sort_keys=True))
