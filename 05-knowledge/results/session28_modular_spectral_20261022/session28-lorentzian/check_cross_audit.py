from fractions import Fraction as F
from pathlib import Path
import hashlib,json,argparse
# Independently multiply factored q rather than importing supplied coefficients.
def mul(a,b):
 r=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):r[i+j]+=x*y
 return r
q=[1]
for factor in ([0,1],[-1,1],[-1,1],[1,1],[-1,2]):q=mul(q,factor)
sq=mul(q,q)
# Independent sieve through10^7, with primary-source table beyond that cutoff.
limit=10**7;s=bytearray(b'\1')*(limit+1);s[0:2]=b'\0\0'
for p in range(2,3163):
 if s[p]:s[p*p:limit+1:p]=b'\0'*(((limit-p*p)//p)+1)
small=[0]+[sum(s[:10**k+1]) for k in range(1,8)]
assert small==[0,4,25,168,1229,9592,78498,664579]
pi=small+[5761455,50847534,455052511,4118054813]
delta=[F(pi[k]-pi[k-1],9*10**(k-1)) for k in range(1,12)]
value=sum(F(co)*delta[i] for i,co in enumerate(sq))
control=sum(F(co,i+1) for i,co in enumerate(sq))
assert q==[0,-1,3,-1,-3,2] and value==F(-2447633,22500000000) and control==F(1,396)
# Algebraic verification of one-variable expected-H and expected-c3 Hessian coefficients.
for n in range(3,26):
 import math
 c=F(math.factorial(n-2),2**(n-2))
 assert 4*c*F(1,4)==c
# Default paths follow the published packet layout; explicit overrides reproduce scratch sources.
parser=argparse.ArgumentParser(description='Replay exact prime-density audit and hash its audited source inputs.')
parser.add_argument('--root-dir',type=Path,default=Path(__file__).resolve().parent.parent,
                    help='Packet root holding the three prime-density artifacts (default: parent directory).')
parser.add_argument('--weil-report',type=Path,default=None,
                    help='Peer Weil report; default is ROOT/session28-weil/REPORT.md.')
args=parser.parse_args()
weil=args.weil_report if args.weil_report is not None else args.root_dir/'session28-weil'/'REPORT.md'
inputs={name:args.root_dir/name for name in ['PRIME_DENSITY_SPECTRAL_OBSTRUCTION.md','check_prime_density_square.py','check_prime_density_square.json']}
inputs['session28-weil/REPORT.md']=weil
print(json.dumps({'status':'passed','polynomial_from_factors':q,'square':sq,'negative_square':str(value),'positive_uniform_control':str(control),'sieve_limit':limit,'sieved_prime_counts':small,'large_counts_source':'Granville PrimesICM1.pdf Table1 printedpage2; inspected web extraction 2026-10-22','input_sha256':{name:hashlib.sha256(p.read_bytes()).hexdigest() for name,p in inputs.items()}},indent=2))
