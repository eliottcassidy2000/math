import itertools,json,math
from pathlib import Path
O=Path(__file__).resolve().parent

def need(c,m):
 if not c:raise ValueError(m)
def mul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
def ev(a,u,v):
 return sum(c*u**(len(a)-1-i)*v**i for i,c in enumerate(a))
A=[[k*x for x in mul(q,q)] for k,q in [(9,[18,50,37]),(36,[1,13,16]),(36,[9,53,64]),(4,[3,-25,-48])]]
R=[[160,-160,-4492,-10700,-7335],[640,5408,18068,27508,15705],[320,5728,27052,48908,30375],[1440,5856,6314,1035],[1440,4704,3626,1035],[1440,7008,9002,1035],[1440,-672,-37484,-93612,-66015],[2880,53088,249356,445548,273375],[5760,51744,174388,258324,141345],[45,498,2023,3498,2160],[45,534,2161,3624,2160],[45,570,2299,3750,2160]]
constants=[18,-9,9,2,-1,-1,2,1,-1,-64,128,-64]
primes=[13,13,13,7,7,7,11,11,11,13,17,13]
cert=json.loads((O/'ap-certificate.json').read_text());expected_cases=[(t,m) for t in itertools.combinations(range(4),3) for m in t]
need(len(cert['cases'])==len(expected_cases)==12,'case count')
checks=projective=0
for n,((triple,mid),f,k,p,c) in enumerate(zip(expected_cases,R,constants,primes,cert['cases'])):
 need(p>1 and all(p%d for d in range(2,math.isqrt(p)+1)),'nonprime')
 need(tuple(c['triple_zero_based'])==triple and c['middle_zero_based']==mid,'wrong index')
 balance=[2*A[mid][j]-sum(A[i][j] for i in triple if i!=mid) for j in range(5)]
 factored=mul([2,3],f) if n in [3,4,5] else f
 need(balance==[k*x for x in factored],'coefficient identity')
 residual=c['residual_certificates'][0]
 need(residual['degree']==len(f)-1 and residual['prime']==p and residual['exponent']==1 and c['constant']==k,'certificate factor metadata')
 need(c['linear_factors']==(['2*u + 3*v'] if n in [3,4,5] else []),'linear list')
 vals=[ev(f,t,1)%p for t in range(p)];inf=ev(f,1,0)%p
 need(vals==residual['affine_residues'] and inf==residual['infinity_residue'],'residue mismatch')
 need(all(vals) and inf,'projective zero');projective+=p+1
 for u in range(p):
  for v in range(p):
   if (u,v)==(0,0):continue
   need(ev(f,u,v)%p!=0,'vector zero');checks+=1
need([ev(a,-3,2) for a in A]==[900,900,12996,900],'boundary')
print(json.dumps({'status':'PASS','coefficient_factor_identities':12,'projective_checks':projective,'nonzero_vector_checks':checks,'boundary':[900,900,12996,900]},indent=2))
