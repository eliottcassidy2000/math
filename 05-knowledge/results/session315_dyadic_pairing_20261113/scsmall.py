from bridge import *
def sc(m,odd):
    v=list(combinations(range(m),2)); N=2*m+odd
    for mask in range(1<<(2*len(v)+(m if odd else 0))):
        A=[0]*N
        def put(i,j,bit=True):
            if not bit:i,j=j,i
            A[i]|=1<<j
        for i in range(m):
            put(i,m+i)
            if odd:
                b=(mask>>(2*len(v)+i))&1
                put(i,2*m,b); put(2*m,m+i,b)
        for k,(i,j) in enumerate(v):
            a=(mask>>(2*k))&1; g=(mask>>(2*k+1))&1
            put(i,j,a); put(j+m,i+m,a)
            put(i,j+m,g); put(j,i+m,g)
        yield mask,A
if __name__=='__main__':
  H={}
  for N in range(1,9):
    if N==1: continue
    counts=Counter()
    for mask,t in sc(N//2,N%2):
        h=sum(table(t)[-1]); counts[h]+=1
        if h not in H:H[h]=(N,mask,t)
    print(N,[ (k,counts[k]) for k in (3,5,7,9,11,13,15,17,19,21,23,25,27)])
  print({h:H.get(h) for h in [11,17,19,23,25,27]})
  support=set()
  for n in range(1,9):
      support |= {sum(table(t)[-1]) for _,t in sc(n//2,n%2)}
  attained={h*k for h in support for k in (1,9,25) if h*k<75}
  holes=set(range(1,75,2))-attained
  if holes!={7,19,21,23,35,39,43,47,53,55,59,63}:
      raise ValueError(holes)
  if sum(table([48,1,3,7,14,30])[-1])!=23:
      raise ValueError('ordinary 23')
  print('global SC holes below 75:', sorted(holes))
