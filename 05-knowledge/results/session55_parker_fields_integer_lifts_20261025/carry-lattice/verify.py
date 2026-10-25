#!/usr/bin/env python3
"""Exact carry-syndrome controls; stdlib only, no assert-based checks."""
import itertools, json, hashlib, pathlib, argparse
P=pathlib.Path(__file__).resolve().parent
PIV=(4,0,2)
NON=(1,3,5,6,7,8)
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def check(ok,msg):
 if not ok: raise RuntimeError(msg)
def magic(c,a,b): return (a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a)
def syndrome(r):
 z=magic(r[4],r[0],r[2]);return tuple(r[i]-z[i] for i in NON)
def sums(r): return tuple(sum(r[i] for i in line) for line in LINES)
def run():
 # Recover every vector with the three free coordinates plus six residuals.
 n=0
 for r in itertools.product(range(-1,2),repeat=9):
  z=list(magic(r[4],r[0],r[2]));d=syndrome(r)
  for i,v in zip(NON,d):z[i]+=v
  check(tuple(z)==r,'integer splitting')
  check((len(set(sums(r)))==1)==(not any(d)),'line equivalence')
  n+=1
 # Universal integral matrix identity checked on standard basis.
 basis=[]
 for j in range(9):
  v=tuple(int(i==j) for i in range(9));basis.append(syndrome(v))
 for j,i in enumerate(NON): check(basis[i]==tuple(int(k==j) for k in range(6)),'surjective primitive syndrome')
 counts={}
 for m in (2,3,4,5,7,9):
  grids={tuple(t%m for t in magic(c,a,b)) for c,a,b in itertools.product(range(m),repeat=3)}
  check(len(grids)==m**3,'parameter count')
  for r in grids:
   d=syndrome(r);check(all(x%m==0 for x in d),'modular membership')
   q=tuple(x//m for x in d)
   # Perturb each selected root vector arbitrarily: error law is linear.
   h=tuple((i*i+3*i)%11 for i in range(9))
   v=tuple(r[i]+m*h[i] for i in range(9))
   check(syndrome(v)==tuple(m*(qj+sj) for qj,sj in zip(q,syndrome(h))),'carry identity')
  counts[str(m)]=len(grids)
 # Weak no-wrap counterexample supplied independently by direct arithmetic.
 r=(28,3,17,5,16,27,15,29,4);p=59
 roots=tuple(next(x for x in range(p) if x*x%p==v) for v in r)
 check(len(set(r))==9 and sums(r)==(48,)*8 and not any(syndrome(r)),'weak no-wrap')
 h=tuple((x*x-v)//p for x,v in zip(roots,r))
 actual=tuple(x*x for x in roots)
 check(syndrome(actual)==tuple(p*v for v in syndrome(h)),'root carry identity')
 check(any(syndrome(actual)),'modular root hostile')
 # Every finite congruence test can accept nine distinct actual squares.
 congruence=[]
 for M in (2,6,30,210,30030):
  roots2=tuple(1+M*i for i in range(1,10));v=tuple(x*x for x in roots2)
  check(len(set(v))==9 and all(d%M==0 for d in syndrome(v)) and any(syndrome(v)),'finite-sieve hostile')
  congruence.append({'modulus':M,'roots':roots2,'syndrome':syndrome(v)})
 return {'integer_vectors_checked':n,'modular_magic_counts':counts,'syndrome_matrix_columns':basis,'weak_no_wrap':{'p':p,'entries':r,'least_roots':roots,'root_carries':h,'root_square_line_sums':sums(actual),'root_square_syndrome':syndrome(actual)},'finite_sieve_hostiles':congruence}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
 out=json.dumps(run(),indent=2,sort_keys=True)+'\n'
 if args.write:(P/'certificate.json').write_text(out)
 else:check((P/'certificate.json').read_text()==out,'certificate drift')
 print('PASS: integral splitting, universal identities, modular counts, carries, hostiles')
if __name__=='__main__':main()
