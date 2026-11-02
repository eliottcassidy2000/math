"""Exact finite sufficient certificates; failure means NOT_CERTIFIED only."""
from profile_tools import covers

def validate_tournament(a):
 n=len(a)
 if any(not isinstance(x,int) for x in a):raise ValueError('integer tournament masks required')
 if n<1 or any(x<0 or x>>n for x in a) or any(a[i]>>i&1 for i in range(n)) or any(((a[i]>>j)&1)+((a[j]>>i)&1)!=1 for i in range(n) for j in range(i+1,n)):raise ValueError('nonempty tournament masks required')

def cyclic_kernel(d):
 if not d or any(not isinstance(x,int) for x in d):raise ValueError('nonempty integer module difference required')
 z=[0]+list(d)+[0,0];n=len(d)+1;K=[[0]*(n+1) for _ in range(n+1)]
 for r in range(1,n+1):
  K[r][r]=z[r-1]+3*z[r]+z[r+1]
  if r<n:K[r][r+1]=K[r+1][r]=z[r]+z[r+1]
 return K

def quadratic(K,v):
 w=[0]+list(v);w+= [0]*max(0,len(K)-len(w))
 return sum(K[i][j]*w[i]*w[j] for i in range(1,len(K)) for j in range(1,len(K)))

def positive_certificate(d):
 K=cyclic_kernel(d)
 if any(x<0 for row in K for x in row):return {'status':'NOT_CERTIFIED','reason':'negative matrix entry'}
 if K[1][1]<=0:return {'status':'WEAK_NONNEGATIVITY_ONLY','K11':K[1][1]}
 return {'status':'SUFFICIENT_CERTIFICATE','scope':'every pair of nonempty side tournaments','lower_bound_coefficient':K[1][1]}

def negative_certificate(d,v):
 K=cyclic_kernel(d);N=len(v)
 if any(not isinstance(x,int) for x in v):raise ValueError('integer base profile required')
 if not N or any(x<=0 for x in v):return {'status':'NOT_CERTIFIED','reason':'base profile lacks positive contiguous support'}
 if any(v[i]*v[i]<v[i-1]*v[i+1] for i in range(1,N-1)):return {'status':'NOT_CERTIFIED','reason':'base profile not log-concave'}
 q=quadratic(K,v)
 if q>=0:return {'status':'NOT_CERTIFIED','reason':'base quadratic is not negative','base_quadratic':q}
 for c in range(1,N+1):
  if any(K[i][j]<0 for i in range(1,min(c,len(K))) for j in range(1,len(K))):continue
  if any(K[i][j]>0 for i in range(c+1,len(K)) for j in range(1,len(K))):continue
  return {'status':'SUFFICIENT_CERTIFICATE','scope':'both side color polynomials are P_base times nonzero monomial-nonnegative polynomials','pivot':c,'base_order':N,'base_quadratic':q}
 return {'status':'NOT_CERTIFIED','reason':'no allowed one-crossing pivot','base_quadratic':q}

def certify_actual(old,new,base=None):
 validate_tournament(old);validate_tournament(new)
 if len(old)!=len(new):raise ValueError('equal module orders required for this fixed-order API')
 f=covers(old);g=covers(new);d=[b-a for a,b in zip(f,g)];out={'actual_profile_difference':d,'positive':positive_certificate(d)}
 if base is not None:
  validate_tournament(base);v=covers(base);out['base_profile']=v;out['negative']=negative_certificate(d,v)
 return out
