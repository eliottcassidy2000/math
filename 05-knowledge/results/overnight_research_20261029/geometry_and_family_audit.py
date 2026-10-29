"""Exhaust quadratic refinements and verify the public family specialization."""
from pathlib import Path
from itertools import combinations,permutations
from collections import Counter
from fractions import Fraction as F
import json,importlib.util

ROOT=Path(__file__).resolve().parent
def require(x,label):
 if not x:raise RuntimeError(label)
def canon(x):return min(x,x^63)
def pairing(x,y):return (x&y).bit_count()%2
def main():
 even=[x for x in range(64) if x.bit_count()%2==0]
 vec=sorted({canon(x) for x in even});nonzero=vec[1:]
 require(len(vec)==16,'quotient dimension four')
 lag={tuple(sorted((x,y,canon(x^y)))) for x,y in combinations(nonzero,2) if not pairing(x,y)}
 require(len(lag)==15,'all isotropic planes')
 types=Counter();profiles=set()
 for A in sorted({canon(x) for x in range(64) if x.bit_count()%2==1}):
  def q(x):return (x.bit_count()//2+(x&A).bit_count())%2
  require(all(q(x)==q(x^63) for x in even),'quadratic quotient descent')
  require(all(q(canon(x^y))==q(x)^q(y)^pairing(x,y) for x in vec for y in vec),'polarization')
  profiles.add(tuple(q(x) for x in vec))
  kept=[s for s in lag if all(q(x)==0 for x in s)];types[len(kept)]+=1
  if A.bit_count()==3:
   require(len(kept)==6,'K33 matchings')
   for row in kept:
    pairs=[x if x.bit_count()==2 else x^63 for x in row]
    require(all((x&A).bit_count()==1 for x in pairs),'crossing-edge selector')
 require(types=={0:6,6:10} and len(profiles)==16,'all quadratic refinements')
 actions=set()
 for p in permutations(range(6)):
  def act(x):return canon(sum(1<<p[i] for i in range(6) if (x>>i)&1))
  action=tuple(act(x) for x in vec);actions.add(action)
  require({tuple(sorted(act(x) for x in row)) for row in lag}==lag,'S6 kernel action')
 require(len(actions)==720,'faithful S6')
 spec=importlib.util.spec_from_file_location('height_producer',ROOT/'fibration_height.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
 # Recompute invariants from the public minimal curve, not from its family coefficients.
 source=ROOT.parent/'session97_overnight_arithmetic_carriers_20261029/ranks/curve4239.json'
 data=json.loads(source.read_text());a1,a2,a3,a4,a6=map(F,data['ainvs'])
 b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6
 c4=b2*b2-24*b4;c6=-b2**3+36*b2*b4-216*b6
 t=F(-2868768,362209);u=F(4,131195359681)
 require(h.evaluate(h.c4,t)==u**4*c4,'named specialization c4')
 require(h.evaluate(h.c6,t)==u**6*c6,'named specialization c6')
 # Independent exact symbolic path over Q, separate from modular Euclid.
 import sympy as sp
 T=sp.symbols('T');f=sp.Poly.from_list(h.c4[::-1],T);g=sp.Poly.from_list(h.c6[::-1],T)
 require(sp.gcd(f,g).degree()==0,'exact rational polynomial gcd')
 require(sp.gcd(f**3,f**3-g**2).degree()==0,'rational j-map degree')
 out={'all_lagrangians':15,'quadratic_refinement_kept_counts':dict(types),'distinct_refinements':16,
      'faithful_S6_actions':720,'family_coprimality_over_Q':True,'curve4239_specialization_invariants_match':True,
      'sympy_version':sp.__version__,'scope':'The rational gcd is an independent arithmetic path; no generic rank or section independence claim.'}
 (ROOT/'geometry_and_family_audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
