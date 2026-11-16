from functools import lru_cache
def data(a):
 n=len(a)
 allmask=(1<<n)-1
 @lru_cache(None)
 def f(S,j):
  if S==1<<j:return 1
  R=S^(1<<j)
  return sum(f(R,i) for i in range(n) if R>>i&1 and a[i]>>j&1)
 @lru_cache(None)
 def H(S):
  return sum(f(S,j) for j in range(n) if S>>j&1) if S else 1
 def profile(x):
  p=[0,0]
  def dfs(S,j):
   p[S.bit_count()&1]+=H(allmask^S)
   for i in range(n):
    if a[j]>>i&1 and not S>>i&1:dfs(S|1<<i,i)
  dfs(1<<x,x)
  return H(allmask),H(allmask^(1<<x)),*p
 return H(allmask),tuple(profile(x) for x in range(n))

def need(ok):
    if not ok:
        raise AssertionError

def apex(w):
    n=len(w)
    return tuple(sum(1<<j for j in range(i+1,n)) + ((w[i]=='0')<<n)
                 for i in range(n))+(sum((c=='1')<<i for i,c in enumerate(w)),)

def extend(a,x):
    n=len(a)
    return tuple(r|((i!=x)<<n) for i,r in enumerate(a))+(1<<x,)

def bridge(u,v,k):
    U,V=apex(u),apex(v)
    m,l=len(U),len(V)
    n=m+k+l
    A=[row|sum(1<<j for j in range(m,n)) for row in U]
    A+=[sum(1<<j for j in range(i+1,n)) for i in range(m,m+k)]
    A+=[row<<(m+k) for row in V]
    A[m-1]^=1<<(n-1)
    A[n-1]|=1<<(m-1)
    return tuple(A)

for k in (1,2,3):
    bank=[
      ("10100"+"1"*k+"00","1"+"0"*(k+1)+"100000",17+63*2**k,33+63*2**(k+1)),
      ("10"+"1"*k+"000","1"*(k+2)+"000",-3+21*2**k,-5+7*2**(k+2)),
      ("1000"+"1"*k+"0110","1000"+"1"*(k+2)+"00",5+81*2**k,7+27*2**(k+2)),
    ]
    for a,b,f,g in bank:
        need(data(apex(a))[0]==f and data(apex(b))[0]==g)
    for u,v,t,w,F,G in [('11010','1000','','110101100',171,405),
                         ('100010','1000','100','1110100',243,525)]:
        need(data(bridge(u,v,k))[0]==F+G*2**k)
        need(data(bridge(t,w,k+1))[0]==2*F-1+G*2**(k+1))

for a,x,S in [
  ((12,21,24,18,1),0,(13,3,6,9)),
  ((24,45,57,48,34,1),0,(23,5,11,16)),
  ((8,65,3,6,79,31,45),3,(41,9,20,29)),
  ((48,77,121,113,98,66,1),0,(59,11,29,40)),
  ((24,77,121,112,98,67,1),0,(75,11,37,48)),
  ((12,77,120,112,99,67,1),0,(95,11,47,58)),
  ((62,124,120,112,32,64,17),4,(27,17,13,22))]:
    need(data(a)[1][x]==S)
    for j in range(2):
        h,g,r,s=S
        S=h+2*r,g,r+s-g,r+s
        a=extend(a,x)
        need(data(a)[1][x]==S)

if __name__=='__main__':
    from itertools import combinations
    from pathlib import Path
    import sys
    if len(sys.argv)>1:
        source=Path(sys.argv[1])
    else:
        source=Path(__file__).resolve().parent.parent/'actual_hamilton_multiplicities_20261021/session24-multiplicities/n7_representatives.txt'
    reps=[tuple(map(int,line.split())) for line in source.read_text().splitlines() if line]
    need(len(reps)==456 and all(len(a)==7 for a in reps))
    masks=set()
    for a in reps:
        for k in range(1,8):
            for S in combinations(range(7),k):
                masks.add(tuple(sum((a[i]>>j&1)<<q for q,j in enumerate(S)) for i in S))
    count=0
    for a in masks:
        h,profiles=data(a)
        if h==27:
            count+=1
            need(not any(2*r==26 and s-r==g for _,g,r,s in profiles))
    print(f"PASS: literal lifts and ports; H27 exclusion on {count} inherited/deleted mask controls")
    import csv
    with source.with_name("n8_classes.csv").open() as f:
        reps8=list(csv.DictReader(f))
    need(len(reps8)==6880)
    count53=0
    extra=[a for a in masks if data(a)[0]==53]
    for row in reps8:
        if int(row["H"])!=53:
            continue
        code=int(row["code"])
        extra.append(tuple((code>>(8*i))&255 for i in range(8)))
    need(all(len(a)==7 for a in extra[:len(extra)-sum(int(row["H"])==53 for row in reps8)]))
    for a in extra:
        h,profiles=data(a)
        need(h==53)
        count53+=1
        need(not any(2*r==52 for _,g,r,s in profiles))
    need(count53>0)
    print(f"PASS: no H53 first-hit port in {count53} inherited atlas controls")
    for K,a,x in [(2,(0,),0),(6,(2,4,1),0),(18,(2,4,25,19,3),1),(54,(4,33,34,7,15,25),2)]:
        h,g,r,s=data(a)[1][x]
        need(2*r==h-1 and h+r+s==K)
    d=[5*2**(1+2*j)*K-1 for j,K in enumerate((54,18,6,2))]
    need(d==[539,719,959,1279] and all(3*y==4*x+1 for x,y in zip(d,d[1:])))

    b=extend((62,124,120,112,32,64,17),4)
    c=list(b)
    c[0]^=2
    c[1]^=1
    c=tuple(r|((not (181>>i&1))<<8) for i,r in enumerate(c))+(181,)
    need(data(b)[0]==53 and data(c)[0]==105)
    need(c==(188,509,248,496,32,192,401,16,181))
