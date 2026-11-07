import json
from math import comb
from fractions import Fraction as F
from pathlib import Path

def rev(n,a):
    return int(f'{n:0{a}b}'[::-1],2) if a else 0

def f(n):
    pos=[i for i in range(n.bit_length()) if n>>i&1]
    return sum(2**(j-i-1)+(j==i+1) for i in pos for j in pos if i<j)

def delta(n,a): return F(n*rev(n,a),1<<a)-f(n)

def runs(n):
    word=bin(n)[2:] if n else ''
    p=0; total=F(0); count=0; i=0
    while i<len(word):
        if word[i]=='0': p*=2;i+=1;continue
        j=i
        while j<len(word) and word[j]=='1':j+=1
        u=F(rev(p,p.bit_length()),1<<p.bit_length())
        length=j-i
        total+=u+(1-u)/(1<<length)
        p=(p<<length)+(1<<length)-1
        count+=1;i=j
    return total,count

word_tests=0
for a in range(14):
    for n in range(1<<a):
        d=delta(n,a); dr,r=runs(n)
        if d!=dr or not (0<=d<=F(3*r,4)) or (n and d<=0):raise Exception(('word',a,n,d,dr))
        word_tests+=1
for m in range(1,41):
    n=(4**m-1)//3
    if delta(n,2*m-1)!=F(2*m,3)-F(2,9)+F(2,9*4**m):raise Exception(('alternating',m))
states=[list(map(int,l.split())) for l in Path(__file__).with_name('CORE_STATES.txt').read_text().splitlines() if l.strip()]
tests=0
for h,q,A,G,T,S,E in states:
 for a in range(9):
  for b in range(9-a):
   for p in range(1<<a):
    for r in range(1<<b):
     t=rev(p,a);s=rev(r,b);v=F(t,1<<a);u=F(r,1<<b)
     alpha=(1<<b)*(q*p+A)+h*r
     gamma=(1<<a)*(q*s+G)+h*t
     Q=q*(1<<(a+b))
     first=(r>>(b-1)) if b else 0
     tau=p*G+h*f(p)+T+(p%2)*S+(q*p+A)*s+h*f(r)+E*first
     D=F(alpha*gamma,Q)-tau
     B=F(A*G,q)+F(h*A,q)*v+F(h*G,q)*u+F(h*h,q)*u*v-T-(p%2)*S-E*first
     if D!=h*(delta(p,a)+delta(r,b))+B:raise Exception('factor')
     if not -T-S-E<=B<=F(A*G+h*A+h*G+h*h,q):raise Exception('bound')
     output=(q*p+A)*((1<<(b+1))-s)+h*(2*r+1-f(r))-p*G-h*f(p)-T-(p%2)*S-E*first
     if output!=alpha*(2-F(gamma,Q))+h+D:raise Exception('output')
     tests+=1
run_count_tests=0
for n in range(15):
 counts=[0]*(n+1)
 for x in range(1<<n): counts[runs(x)[1]]+=1
 for j in range(n+1):
  expected=comb(n+1,2*j) if 2*j<=n+1 else 0
  if counts[j]!=expected:raise Exception(('run-count',n,j))
  run_count_tests+=1
print(json.dumps({'status':'PASS','binary_word_tests':word_tests,'alternating_words':40,'core_padding_tests':tests,'run_count_tests':run_count_tests,'coverage_claim':False},indent=2))
