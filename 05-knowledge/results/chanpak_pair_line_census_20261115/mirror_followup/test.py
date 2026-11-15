import collections, math
def sig(w,n):
 b=[(w>>i)&1 for i in reversed(range(n))]
 rev=sum(x<<i for i,x in enumerate(b))
 S=w-rev;E=b[-1]-b[0]
 I=0; prefix=0
 for j,z in enumerate(b):
  if not z and j:
   I+=prefix//2+b[j-1]
  prefix=2*prefix+z
 return S,E,I
def run(lim):
 for n in range(3,lim+1):
  shells=collections.defaultdict(dict)
  for w in range(1<<n):
   s,e,i=sig(w,n);shells[s,e][i]=w
  best=(0,0,0,0)
  for (s,e),ij in shells.items():
   for i in ij:
    if i-1 in ij:continue
    j=i
    while j+1 in ij:j+=1
    best=max(best,(j-i+1,s,e,i))
  k,s,e,i=best
  print(n,k,s,e,i,[format(shells[s,e][j],f'0{n}b') for j in range(i,i+k)],flush=True)
if __name__=='__main__':
 import sys
 run(int(sys.argv[1]))
