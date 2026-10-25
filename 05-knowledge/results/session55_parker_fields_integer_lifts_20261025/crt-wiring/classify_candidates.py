#!/usr/bin/env python3
"""Finite CRT exhaustion conditional on separately proved prime-power classification."""
import argparse,json,pathlib
from verify import *
ODD_PRIME_PARKER=[3,5,7,11,13,17,19,23,31,43,47,67]
R=ODD_PRIME_PARKER[1:]+[25,49]
ODD_PARKER=sorted([3,9]+R+[3*x for x in R])

def build():
 local={q:minimal(enumerate_masks(q)) for q in sorted(set([2**a for a in range(1,7)]+[q for n in ODD_PARKER for q in factors(n)]))}
 rows=[]
 for odd in [1]+ODD_PARKER:
  for a in range(7):
   n=2**a*odd;fs=factors(n);states={ALL:[]}
   for q in fs:states=compose(states,local[q])
   w=None
   if 0 in states:
    ls=states[0];m=tuple(crt(v,fs) for v in zip(*ls))
    if not valid(n,m) or len(set(m))!=9:raise RuntimeError('invalid CRT witness')
    w={'entries':m,'square_roots':[roots(n)[x] for x in m],'local_arrays':ls}
   rows.append({'n':n,'dyadic_exponent':a,'odd_part':odd,'parker':0 not in states,'minimal_masks':[hex(x) for x in states],'witness':w})
 rows.sort(key=lambda r:r['n']);P=[r['n'] for r in rows if r['parker']]
 return {'scope':'All remaining candidates after prime-power and odd-part reduction; dependency theorem must be supplied separately','odd_prime_parker':ODD_PRIME_PARKER,'odd_parker':ODD_PARKER,'candidate_count':len(rows),'parker_count_including_1':len(P),'largest_parker':max(P),'parker_moduli':P,'local_antichains':{str(q):[{'mask':hex(x),'entries':m} for x,m in D.items()] for q,D in local.items()},'rows':rows}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 out=build();target=pathlib.Path(__file__).with_name('classification_candidates.json')
 if args.verify:
  if json.loads(target.read_text())!=json.loads(json.dumps(out)):raise RuntimeError('candidate certificate mismatch')
 else:target.write_text(json.dumps(out,indent=2)+'\n')
 print({k:v for k,v in out.items() if k not in ['rows','local_antichains']})
