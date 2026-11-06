import itertools,json
S={0:{1}}; A={1:{1},2:set()}; minima={1:1}
for n in range(1,6):
 s=set(); a=set(); pairs=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  out=[0]*n
  for k,(i,j) in enumerate(pairs):
   x,y=(i,j) if bits>>k&1 else (j,i)
   out[x]|=1<<y
  h=sum(all(out[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in itertools.permutations(range(n)))
  s.add(h)
  reach=[o|1<<i for i,o in enumerate(out)]
  for k in range(n):
   for i in range(n):
    if reach[i]>>k&1:reach[i]|=reach[k]
  if all(x==(1<<n)-1 for x in reach):a.add(h)
 S[n]=s; A[n]=a
 predicted={x*y for r in range(1,n+1) for x in A[r] for y in S[n-r]}
 assert predicted==s,(n,predicted,s)
 for h in s:minima.setdefault(h,n)
 print(n,'spectrum',sorted(s),'strong',sorted(a))
assert minima=={1:1,3:3,5:4,9:5,11:5,13:5,15:5}
print('PASS: exhaustive independent permutation counts and SCC spectrum recurrence through n=5; minima',minima)
