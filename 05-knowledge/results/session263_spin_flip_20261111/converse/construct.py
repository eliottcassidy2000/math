#!/usr/bin/env python3
"""Literal new converse operations and finite independent DP controls."""
import importlib.util,json,itertools
from pathlib import Path
P=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('ordinal',P.parents[1]/'session262_interface_supply_20261111/constructor/ordinal_interfaces.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def op(T):return tuple(sum(1<<j for j in range(len(T)) if T[j]>>i&1) for i in range(len(T)))
def flip(p):return {k:p[{'S':'E','E':'S','U':'V','V':'U'}.get(k,k)] for k in p}
def neutral(T,O):
 T,O=m.ordinal((0,),0,T,O)
 return m.ordinal(T,O,(0,),0)
def word(bits):
 A=((10,12,9,0),8) # C3 => marked sink
 return m.compile_word([(A if x==0 else (op(A[0]),A[1])) for x in bits])
SPECIES=(2,4,8,16,32,64)
def triangle(n):
 if n<1:raise ValueError('nonempty transitive module required')
 M=tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
 return tuple(x|(1<<n) for x in M)+(1<<(n+1),(1<<n)-1)
def allocation_source(counts,left):
 """Emit arbitrary six-species source and exact profile without full-graph DP."""
 if len(counts)!=6 or len(left)!=6 or any(not isinstance(x,int) for x in [*counts,*left]) or any(n<0 or a<0 or a>n for n,a in zip(counts,left)):
  raise ValueError('need six nonnegative counts and six bounded left counts')
 empty=lambda h,q:dict(H=h,Q=q,S=0,E=0,U=0,V=0,W=0)
 atom=[(triangle(r),0,empty(2**r+1,2*3**r+4*2**r-2)) for r in SPECIES]
 cap=((0,),0,empty(1,2));root=((0,),1,dict(H=1,Q=2,S=1,E=1,U=1,V=1,W=0))
 seq=[cap]+[atom[i] for i,a in enumerate(left) for _ in range(a)]+[root]+[atom[i] for i,(n,a) in enumerate(zip(counts,left)) for _ in range(n-a)]+[cap]
 T,O,p=seq[0]
 for B,OB,pb in seq[1:]:
  T,O=m.ordinal(T,O,B,OB);p=m.combine(p,pb)
 return T,O,p
def main():
 bank=[((0,),1),((2,0),1),((2,4,1),1),((10,12,9,0),8),((2,4,9,19,7),9)]
 controls=[]
 for T,O in bank:
  p=m.profile(T,O);assert m.profile(op(T),O)==flip(p)
  R,RO=m.ordinal(T,O,op(T),O);r=m.profile(R,RO)
  expected=dict(H=p['H']**2,Q=p['Q']**2,S=p['S']*p['H'],E=p['S']*p['H'],U=p['U']*p['Q']+p['H']*p['V'],V=p['U']*p['Q']+p['H']*p['V'],W=2*p['H']*p['W']+p['U']**2+p['E']**2)
  assert r==expected and r['H']-r['W']<=0
  N,NO=neutral(T,O);n=m.profile(N,NO)
  assert n==dict(H=p['H'],Q=4*p['Q'],S=0,E=0,U=2*p['U'],V=2*p['V'],W=p['W'])
  controls.append(dict(seed_rows=T,marks=O,profile=p,double=r,neutral=n))
 words=[]
 for k in range(1,4):
  seen=set()
  for bits in itertools.product([0,1],repeat=k):
   T,O,p=word(bits);literal=m.profile(T,O);assert p==literal
   U=3**k*sum((1+3*x)*8**(k-1-i) for i,x in enumerate(bits))
   V=3**k*sum((4-3*x)*8**i for i,x in enumerate(bits))
   assert (p['H'],p['U'],p['V'])==(3**k,U,V)
   assert U not in seen;seen.add(U)
   if k>=2:assert p['H']-p['W']<0
   words.append(dict(bits=bits,H=p['H'],U=U,V=V,D=p['H']-p['W']))
 # Single-mark scalar-identical allocation bank: TT1 and triangle chain blowup.
 F=triangle
 alloc=[]
 atoms=[F(2),F(4)]
 for bits in itertools.product([0,1],repeat=2):
  left=[(a,0) for a,b in zip(atoms,bits) if b==0];right=[(a,0) for a,b in zip(atoms,bits) if b==1]
  T,O,p=m.compile_word(left+[((0,),1)]+right)
  assert p==m.profile(T,O) and p['W']==0
  capped,mark,cp=allocation_source([1,1,0,0,0,0],[1-bits[0],1-bits[1],0,0,0,0])
  assert cp['U']==2*p['U'] and cp['V']==2*p['V']
  assert cp==m.profile(capped,mark) and cp['H']==cp['H']-cp['W'] and cp['S']==cp['E']==0
  alloc.append(dict(bits=bits,H=p['H'],U=p['U'],V=p['V'],capped_profile=cp))
 # Exact valuation certificate, no factorization oracle needed: trial primality.
 ms=[2,4,8,16,32,64];primes=[5,17,257,65537,641,274177]
 def vp(x,p):
  c=0
  while x%p==0:c+=1;x//=p
  return c
 assert all(all(p%d for d in range(2,int(p**.5)+1)) for p in primes)
 values=[(2**r+1,2*3**r+4*2**r-2) for r in ms]
 matrix=[[vp(q,p)-vp(h,p) for h,q in values] for p in primes]
 expected=[[-1 if i==j else (1 if i==1 and j==2 else 0) for j in range(6)] for i in range(6)]
 assert matrix==expected
 certificate=dict(module_lengths=ms,prime_rows=primes,H_Q=values,valuation_matrix=matrix,determinant=1)
 receipt=dict(valuation_certificate=certificate,status='FINITE-EXACT new operations only',converse_controls=controls,binary_word_controls=words,allocation_controls=alloc)
 (P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
