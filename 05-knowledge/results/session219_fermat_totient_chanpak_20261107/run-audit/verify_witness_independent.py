import json,math
from pathlib import Path
P,Q=641,6700417

def need(v,msg):
    if not v: raise RuntimeError(msg)

def step(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)

def word(w,s=(1,0,2,1)):
    for c in w:s=step(s,c)
    return s

def F(s):return s[0]-s[1]+s[2]

def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
need(prime(P) and prime(Q),'prime factors')
need(P*Q==2**32+1,'F5 factor')
for p in (P,Q):
    need(pow(2,32,p)==p-1 and pow(2,64,p)==1,'order64')
# Exhaustive small/state-independent identity checks, including zero representative as length64.
for D in range(1,5):
 for E in range(4):
  for U in range(1,5):
   for d in range(1,4):
    s=(D,E,U,d)
    for k in (1,2,3,31,64,65):
     need(word('O'*k,s)==(D,0,2**k*U,d),'O closed form')
     need(word('Z'*k,s)==(D-E+(2**k-1)*U-(2**k-2)*d,d,d+2**k*(U-d),d),'Z closed form')
     for c in 'OZ':
      need(tuple(x%(P*Q) for x in word(c*k,s))==tuple(x%(P*Q) for x in word(c*(k+64),s)),'positive run period')
w='O'*31+'Z'*15+'O'*9+'Z'*5+'O'*9
s=word(w);h=F(s)
need(h==1182825546590390812645 and h==176530139331685*Q,'five-run witness')
need(h%P==275 and 2+len(w)+1==72,'five-run residue/order')
rest='Z'*15+'O'*9+'Z'*5+'O'*9
beta=F(word(rest,(1,0,0,1))); alpha=F(word(rest,(1,0,1,1)))-beta
need(alpha>0 and h==alpha*2**32+beta,'infinite five-run family')
prefix='OOZJOZJJOJJOZ';sp=word(prefix)
need(sp==(4909,859,6013,859),'root-change prefix')
hr=F(word(prefix+'O'*20))
need(hr==6305092397==941*Q and hr%P==380,'root-change witness')
need(2+len(prefix)+20+1==36,'root-change order')
for t in range(6):
    ht=alpha*2**(32+64*t)+beta
    rt=4909+6013*2**(20+64*t)
    need(ht%Q==0 and ht%P==275,'five-run family residue')
    need(rt%Q==0 and rt%P==380,'root-change family residue')
out={'status':'PASS','five_run':{'word_runs':[['O',31],['Z',15],['O',9],['Z',5],['O',9]],'state':s,'H':h,'order':72,'q_quotient':h//Q,'mod641':h%P,'family_formula':{'alpha':alpha,'beta':beta,'exponent':'32+64t','order':'72+64t'}},'root_change':{'prefix':prefix,'state':sp,'H0':hr,'order0':36,'mod641':hr%P},'scope':'exact arithmetic and inherited physical grammar; no independent enumeration of all paths in order36 or72 graphs'}
Path(__file__).with_name('witness_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
