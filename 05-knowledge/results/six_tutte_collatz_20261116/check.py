from itertools import permutations
def require(ok):
    if not ok:
        raise AssertionError
def starts(a):
    n=len(a)
    require(all(not (a[i]>>i&1) for i in range(n)))
    require(all((a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(n) for j in range(i)))
    counts=[0]*n
    for p in permutations(range(n)):
        if all(a[i]>>j&1 for i,j in zip(p,p[1:])):
            counts[p[0]]+=1
    return tuple(counts)
for h,a in [(11,(6,28,24,17,1)),(13,(12,21,24,18,1)),(27,(30,52,56,50,32,1)),(41,(14,28,56,48,33,3))]:
    require(sum(starts(a))==h)
a=(50,252,225,501,420,128,177,257,103)
require(starts(a)==(33,81,53,224,50,18,51,18,39))
require(all((a[i]>>j&1)==(a[8-j]>>(8-i)&1) for i in range(9) for j in range(9)))
require(any(all(a[i]>>j&1 for i,j in zip(p,p[1:]+p[:1])) for p in permutations(range(9))))
p=tuple(sum(1<<((i+r)%7) for r in [1,2,4]) for i in range(7))
require(sum(starts(p))==189)
def v2(n): return (n&-n).bit_length()-1
def verify_colouring(q,N):
    colors={}
    def put(x,y,color): colors[tuple(sorted((x,y)))]=color
    def up(n,sign):
        return (0 if q!=3 or n==1 else ((n-1)//2).bit_length()%2)^(sign<0)
    for n in range(1,N+1,2):
        for sign in [-1,1]:
            m=q*n+sign
            if m<=N: put(n,m,up(n,sign))
    for m in range(1,N//2+1):
        c=2
        if v2(m)%2:
            want=[]
            for x in (m,2*m):
                for sign in [-1,1]:
                    if (x-sign)%q==0:
                        n=(x-sign)//q
                        if n>0 and n%2 and (n,x)!=(1,2):
                            want.append(up(n,sign))
            require(len(set(want))<=1)
            c=1-want[0] if want else 0
        put(m,2*m,c)
    used={}
    for (x,y),c in colors.items():
        for v in (x,y):
            s=used.setdefault(v,set())
            require(c not in s)
            s.add(c)
for q in (3,5,7,9,11,15):
    verify_colouring(q,60000)
require([pow(12,6,49),pow(3,6,49)]==[22,43])
print("Witnesses, anti-automorphism, six clocks and six multiplier colourings pass.")

a=(148,125,120,241,96,193,1,86)
require(starts(a)==(7,68,29,33,13,13,7,19))
require(all((a[i]>>j&1)==(a[7-j]>>(7-i)&1) for i in range(8) for j in range(8)))
require(any(all(a[i]>>j&1 for i,j in zip(p,p[1:]+p[:1])) for p in permutations(range(8))))
print("The even-order H189 witness passes.")
