from math import isqrt
import json

def fa(n):
 if n==0:return 0
 return (3*n-isqrt(5*n*n)-1)//2

def floor_phi(n):return (n+isqrt(5*n*n))//2

def floor_beta(n):return (3*n+isqrt(5*n*n))//2

def S(n):
 ans=0;steps=0
 while n:
  k=fa(n);ans+=k*(n+1)-3*k*(k+1)//2;n=k;steps+=1
 return ans,steps

sigma={'C':'CA','A':'B','B':'CA'};w='C'
while len(w)<20000:w=''.join(sigma[x] for x in w)
counts=dict.fromkeys('CAB',0);tri=dict.fromkeys('CAB',0);checks=0
for n,c in enumerate(w[:20000],1):
 counts[c]+=1
 expected={'C':fa(n+2),'A':fa(n+1)};expected['B']=n-sum(expected.values())
 if counts!=expected:raise ValueError('prefix counts')
 for x in tri:tri[x]+=counts[x]
 ctri=S(n+2)[0];atri=S(n+1)[0];btri=n*(n+1)//2-ctri-atri
 if tri!={'C':ctri,'A':atri,'B':btri}:raise ValueError('triangle counts')
 if ctri-atri!=counts['C']:raise ValueError('boundary difference')
 checks+=3
for k in range(1,4000):
 for c,pos in [('C',floor_beta(k)-1),('A',floor_beta(k)),('B',2*floor_phi(k)+k)]:
  if w[pos-1]!=c:raise ValueError('Beatty position')
  checks+=1
for n in range(1000):
 if S(n)[0]!=sum(fa(i) for i in range(1,n+1)):raise ValueError('floor sum')
 checks+=1
# Matrix quotient and elementary factor projection, with no numeric approximation.
M=((1,0,1),(1,0,1),(0,1,0))
def mv(v):return tuple(sum(row[i]*v[i] for i in range(3)) for row in M)
if mv((1,0,-1))!=(0,0,0):raise ValueError('killed mode')
for x in range(-3,4):
 for y in range(-3,4):
  if mv((x,x,y))!=(x+y,x+y,x):raise ValueError('image quotient')
  checks+=1
for v in [(1,0),(0,1),(1,1)]:
 x=v;orbit=[]
 for _ in range(3):orbit.append(x);x=((x[0]+x[1])%2,x[0])
 if x!=v or len(set(orbit))!=3:raise ValueError('F2 quotient action')
 checks+=1
projected={(a,b) for a in [0,1] for b in [0,1] if not a*0 and not 0*b}
reconnected={(a,b) for a in [0,1] for b in [0,1] if not a*b}
if projected-reconnected!={(1,1)}:raise ValueError('zero elimination')
large=[]
for n in [10**20,10**100,10**500]:
 value,steps=S(n);k=fa(n)
 if value!=k*(n+1)-3*k*(k+1)//2+S(k)[0]:raise ValueError('large recursion')
 # Exact radical bracketing equivalent to k < n alpha < k+1.
 if not (3*n-2*k)**2>5*n*n or not (3*n-2*(k+1))**2<5*n*n:raise ValueError('large floor')
 large.append(dict(n=str(n),S=str(value),steps=steps))
print(json.dumps(dict(status='PASS',checks=checks,canonical_prefix=w[:36],large_controls=large,zero_projected_pairs=sorted(projected),zero_reconnected_pairs=sorted(reconnected)),indent=2))
