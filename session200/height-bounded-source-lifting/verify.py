"""Explicit quadratic circuit replay, modular closure, and CRT hostile."""
from itertools import product
from math import isqrt, prod
from pathlib import Path
import hashlib,json

T=(1,0,2,1)
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
A={
'I':((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)),
'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))}
LET='IOZJ'

def check(q,msg):
 if not q: raise RuntimeError(msg)
def dot(v,w):return sum(x*y for x,y in zip(v,w))
def step(s,c):
 D,E,U,d=s
 if c=='I':return s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D
 raise ValueError(c)
def K(v,w):return dot(v,[dot(row,w) for row in Q])
def cap(n):return isqrt(4*(n-1))
def limit(B):
 L=0
 while 5*3**(L+1)<=3*B*2**(L+1):L+=1
 return L

def bank(B):
 out={T:''} if B>=2 else {};queue=list(out)
 for s in queue:
  for c in 'OZJ':
   t=step(s,c)
   check(t[2]>s[2],'cover grows')
   if t[2]<=B and t not in out:out[t]=out[s]+c;queue.append(t)
 return out

def binary(x,ell):return [(x>>i)&1 for i in range(ell)]
def circuit(n,w1,w2):
 B=cap(n);L=limit(B);ell=B.bit_length();data={}
 for side,word in enumerate((w1,w2)):
  check(len(word)<=L,'word limit')
  s=T
  for i,c in enumerate(word+'I'*(L-len(word))+'!'):
   for j,x in enumerate(s):
    check(0<=x<=B,'range witness')
    for kind,z in [('b',x),('c',B-x)]:
     for h,v in enumerate(binary(z,ell)):data[f'{side}.{i}.{j}.{kind}.{h}']=v
   if c!='!':
    for a in LET:data[f'{side}.{i}.e.{a}']=int(a==c)
    s=step(s,c)
 return data

def equations(n,d):
 B=cap(n);L=limit(B);ell=B.bit_length();out=[];ends=[]
 for z in d.values():out.append(z*(z-1))
 for side in range(2):
  states=[]
  for i in range(L+1):
   x=[]
   for j in range(4):
    X=sum((1<<h)*d[f'{side}.{i}.{j}.b.{h}'] for h in range(ell))
    Y=sum((1<<h)*d[f'{side}.{i}.{j}.c.{h}'] for h in range(ell))
    out.append(X+Y-B);x.append(X)
   states.append(x)
  out.extend(a-b for a,b in zip(states[0],T))
  for i in range(L):
   e=[d[f'{side}.{i}.e.{a}'] for a in LET]
   out.append(sum(e)-1)
   out.extend(e[j]*e[k] for j in range(4) for k in range(j+1,4))
   for j in range(4):out.append(states[i+1][j]-sum(e[k]*dot(A[a][j],states[i]) for k,a in enumerate(LET)))
  ends.append(states[-1])
 out.append(K(*ends)-n)
 return out

def prime_after(n):
 p=n+1
 while any(p%d==0 for d in range(2,isqrt(p)+1)):p+=1
 return p

def modular_states(n,p):
 B=cap(n);L=limit(B);R=(1<<B.bit_length())-1
 allowed={x%p for x in range(R+1) if any((x+y-B)%p==0 for y in range(R+1))}
 states={tuple(x%p for x in T)}
 if any(x not in allowed for x in T):return set(),allowed
 for _ in range(L):
  states={q for s in states for a in LET for q in [tuple(x%p for x in step(s,a))] if all(x in allowed for x in q)}
 return states,allowed

def crt(vals,mods):
 M=prod(mods)
 return sum(a*(M//p)*pow(M//p,-1,p) for a,p in zip(vals,mods))%M

def main():
 report={};N=401;witnesses=0;modchecks=0
 for n in range(3,N+1):
  B=cap(n);b=bank(B);L=limit(B);p=prime_after(4*n)
  check(all(len(w)<=L for w in b.values()),'all physical words bounded')
  values={K(v,w):(b[v],b[w]) for v in b for w in b}
  check(all(0<x<=B*B for x in values),'physical scalar cap')
  for s in b:
   check(all(0<=x<=B for x in s),'physical cone')
   for a in LET:check(step(s,a)==tuple(dot(r,s) for r in A[a]),'independent transition matrix')
  if n in values:
   d=circuit(n,*values[n]);eq=equations(n,d)
   check(not any(eq),'integer circuit')
   check(not any(x%p for x in eq),'prime circuit')
   check(len(d)==16*B.bit_length()*(L+1)+8*L,'variable count')
   witnesses+=1
  check((n in values)==any((v-n)%p==0 for v in values),'one-prime scalar equivalence')
  if n<=101:
   states,allowed=modular_states(n,p)
   check(allowed==set(range(B+1)),'finite field binary box exact')
   check(states==set(b),'complete field circuit trace inventory')
   check(any((K(v,w)-n)%p==0 for v in states for w in states)==(n in values),'complete finite field circuit acceptance')
   modchecks+=1
 report['integer_targets_checked_through']=N
 report['actual_integer_circuit_witnesses']=witnesses
 report['complete_modular_circuit_target_checks']=modchecks
 # Exhaustive formal box residual bound, independent of physical grammar.
 residuals=0
 for B in range(2,9):
  for s in product(range(B+1),repeat=4):
   for a in LET:
    t=step(s,a)
    for z in t:check(-2*B<=-z<=2*B and -2*B<=B-z<=2*B,'selected residual interval')
    residuals+=4
 report['formal_box_residual_checks']=residuals
 # Explicit terminal-only wrap at a prime that exceeds the transition bound.
 n=67;p=67;d=circuit(n,'OZO','OOO');eq=equations(n,d)
 check(not any(x%p for x in eq),'small prime false positive')
 check(all(x==0 for x in eq[:-1]) and eq[-1]==134,'only terminal equation fails integrally')
 b=bank(cap(n));check(all(K(v,w)!=n for v in b for w in b),'C1 gap67 complete cutoff')
 report['terminal_wrap_hostile']={'target':67,'prime':67,'cover_cap':16,'words':['OZO','OOO'],'actual_count':201,'terminal_integer_residual':134}
 # A complete bounded CRT hostile, assembled at polynomial-variable level.
 n=21;mods=[3,5,7,13];words=[('',''),('OZ',''),('OZ','ZJ'),('OZ','OO')]
 ds=[circuit(n,*w) for w in words]
 for p,d in zip(mods,ds):check(not any(x%p for x in equations(n,d)),'CRT local circuit')
 M=prod(mods);joined={key:crt([d[key] for d in ds],mods) for key in ds[0]}
 check(not any(x%M for x in equations(n,joined)),'all combined polynomial equations')
 check(joined['0.0.e.I']==910 and joined['0.0.e.O']==456,'nontrivial selectors')
 b=bank(cap(n));values=sorted({K(v,w) for v in b for w in b})
 check(values==[3,5,9,11,15,17,19,23,25,27,29,31,33,35,37,39,45,47,51],'complete gap21 bank')
 check(21 not in values and M>cap(21)**2,'CRT high product false positive')
 check(all(not all((h-21)%p==0 for p in mods) for h in values),'no common physical witness')
 report['crt_hostile']={'target':21,'cover_cap':8,'length_cap':limit(8),'primes':mods,'modulus':M,'words':words,'values':[3,11,35,47],'first_left_selectors':{a:joined[f'0.0.e.{a}'] for a in LET},'all_polynomial_equations_hold_mod_product':True,'integer_target_absent':True,'complete_pair_values':values}
 report['proof_sha256']=hashlib.sha256(Path(__file__).with_name('PROOF.md').read_bytes()).hexdigest()
 report['verifier_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 report['status']='PASS: all-order lifting proved in text; finite controls exact'
 print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':main()
