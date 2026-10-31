import itertools,json
checks=0
def ck(v):
 global checks
 checks+=1
 if not v:raise RuntimeError('identity failed')
for Q in (3,5,9,25,49,121):
 h=(Q-1)//2
 for n in range(0,201):
  for r in range(Q):
   d=h+r;N=Q*n+d
   ck(2*N+1==Q*(2*n+1)+2*r)
   ck((N-h)//Q==n and N-Q*n==d)
   if Q in(9,25,49,121):ck(N%2==(n+r)%2)
for a in(32,41,42,53,112):
 for length in range(7):
  for word in itertools.product((25,49),repeat=length):
   n=a;prod=1
   for Q in word:n=Q*n+(Q-1)//2;prod*=Q
   ck(n==((2*a+1)*prod-1)//2)
ck(49*41+24==2033);ck(49*42+24==2082);ck(25*42+12==1062);ck(49*32+24==1592)
print(json.dumps({'status':'PASS','checks':checks,'radices':[3,5,9,25,49,121],'mixed_depth':6},sort_keys=True))
