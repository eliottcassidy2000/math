#!/usr/bin/env python3
"""
fib_tri_square_prime_heron — mac-mini-2026-06-16-S1

Independent master verification of the Fibonacci / triangular / square / prime
constellation raised by the human, with the tournament-parity PROJECT HOOK that
the staircase tile-count m = C(n-1,2) = T_{n-2} is a TRIANGULAR number.

Threads verified:
  T0  m = C(n-1,2) = T_{n-2}; which n make m Fibonacci (n in {3,4,8,12})
  T1  Heronian m=1: exactly 5 integer triangles with Area=Perimeter; all inradius 2;
       general Area = M*Perimeter  <=>  inradius r = 2M.
  T2  octal odd squares: (2k+1)^2 = 8*T_k + 1  =>  octal = (octal T_k) followed by '1'
  T3  Fibonacci INTERSECT triangular = {1,3,21,55};  (2n-1)(4n-1) = T_{4n-2}
  T4  2x6 grid marbles up to Klein-four (Z2xZ2): k=0..3 -> 1,3,21,55; k=4 -> 135
  T5  Fibonacci rank of apparition alpha(p) | p - (5/p)  (Euler's criterion on 5)
  T6  Heron square = (x+y+z)*x*y*z (product of 4 naturals, glue s=x+y+z);
       Lagrange 4 squares; Gauss 3 triangular (EYPHKA)
No external deps; pure stdlib.
"""
from math import isqrt, gcd
from itertools import combinations

def line(t): print("\n" + "="*70 + "\n" + t + "\n" + "="*70)

# ---------------------------------------------------------------- T0
line("T0  m = C(n-1,2) = T_{n-2}; when is the tile-count Fibonacci?")
def C(n,k):
    if k<0 or k>n: return 0
    num=1
    for i in range(k): num=num*(n-i)//(i+1)
    return num
def T(k): return k*(k+1)//2
fibset=set(); a,b=1,1
while a< 10**6:
    fibset.add(a); a,b=b,a+b
fibset.add(0)  # F(0)=0
rows=[]
for n in range(3,16):
    m=C(n-1,2); tri=T(n-2)
    assert m==tri, (n,m,tri)
    rows.append((n,m,m in fibset))
print(" n  m=C(n-1,2)=T_{n-2}  Fibonacci?")
for n,m,f in rows:
    print(f"{n:2d}  {m:6d}             {'YES <--' if f else ''}")
fib_n=[n for n,m,f in rows if f]
print("tile-count m is BOTH triangular AND Fibonacci at n =", fib_n,
      "(m =", [T(n-2) for n in fib_n], ")")

# ---------------------------------------------------------------- T1
line("T1  Heronian m=1: Area = Perimeter  <=>  inradius 2; exactly 5")
def heron_area_sq16(a,b,c):
    # 16*Area^2 = (a+b+c)(-a+b+c)(a-b+c)(a+b-c)
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)
def integer_area(a,b,c):
    A2_16=heron_area_sq16(a,b,c)
    if A2_16<=0: return None
    if A2_16 % 16: return None
    A2=A2_16//16
    r=isqrt(A2)
    return r if r*r==A2 else None
def perfect_triangles(M, cap=400):
    out=[]
    for a in range(1,cap):
        for b in range(a,cap):
            for c in range(b,a+b):
                if c>=cap: break
                A=integer_area(a,b,c)
                if A is None: continue
                if A==M*(a+b+c):
                    out.append((a,b,c,A,a+b+c))
    return out
five=perfect_triangles(1)
print("m=1 (Area=Perimeter):")
for a,b,c,A,P in five:
    s=P/2; r=A/s
    print(f"  ({a:2d},{b:2d},{c:2d})  area={A}  perim={P}  inradius r=A/s={r:g}")
print("count =", len(five), "(expected 5)")
print("all inradius==2 :", all(abs(A/(P/2)-2)<1e-9 for *_,A,P in [(0,)+t for t in five]))
print("counts of integer triangles with r=2M (M=1..5):",
      {M: len(perfect_triangles(M)) for M in range(1,6)})

# ---------------------------------------------------------------- T2
line("T2  octal odd squares:  (2k+1)^2 = 8*T_k + 1")
def oct_str(x): return format(x,'o')
ok=True
for k in range(0,2001):
    if (2*k+1)**2 != 8*T(k)+1: ok=False; break
print("(2k+1)^2 == 8*T_k+1 for k=0..2000:", ok)
print(" k   (2k+1)^2  octal     T_k  octal(T_k)  octal(sq)==octal(T_k)+'1'")
for k in range(0,12):
    sq=(2*k+1)**2
    chk = oct_str(sq)==oct_str(T(k))+'1' if k>0 else oct_str(sq)=='1'
    print(f"{k:2d}  {sq:7d}  {oct_str(sq):>7}   {T(k):4d}  {oct_str(T(k)):>6}      {chk}")
print("=> every odd square is 1 mod 8; strip the octal '1' -> triangular number.")

# ---------------------------------------------------------------- T3
line("T3  Fibonacci INTERSECT triangular = {1,3,21,55};  (2n-1)(4n-1)=T_{4n-2}")
def is_tri(x):
    if x<0: return False
    d=8*x+1; r=isqrt(d); return r*r==d
fibs=[]; a,b=1,1
seen=set()
for _ in range(20000):
    fibs.append(a); a,b=b,a+b
    if a>10**4000: break
tri_fibs=sorted({f for f in fibs if is_tri(f)})
print("Fibonacci numbers that are triangular (searched F up to ~10^4000):", tri_fibs)
print(" n  (2n-1)(4n-1)  == T_{4n-2}?  Fibonacci?")
for n in range(0,7):
    v=(2*n-1)*(4*n-1); t=T(4*n-2)
    print(f"{n:2d}  {v:6d}        {v==t}          {'FIB' if v in fibset else ''}")

# ---------------------------------------------------------------- T4
line("T4  2x6 grid, k marbles up to Klein-four {e, hflip, vflip, 180}")
cells=[(r,cc) for r in range(2) for cc in range(6)]
idx={c:i for i,c in enumerate(cells)}
def hflip(c): r,cc=c; return (r,5-cc)      # reverse columns
def vflip(c): r,cc=c; return (1-r,cc)      # swap rows
def rot(c):   r,cc=c; return (1-r,5-cc)    # 180
syms=[lambda c:c, hflip, vflip, rot]
perm=[[idx[s(c)] for c in cells] for s in syms]
def burnside(k):
    tot=0
    for p in perm:
        # count k-subsets fixed by permutation p = choose orbits fully in/out
        # build cycles of p
        seen=[False]*12; cyc=[]
        for i in range(12):
            if not seen[i]:
                j=i; l=0
                while not seen[j]:
                    seen[j]=True; j=p[j]; l+=1
                cyc.append(l)
        # subsets fixed: pick a subset of cycles whose total length = k
        from functools import lru_cache
        # dp over cycle lengths
        dp=[0]*(13); dp[0]=1
        for l in cyc:
            ndp=dp[:]
            for s in range(12,-1,-1):
                if dp[s] and s+l<=12: ndp[s+l]+=dp[s]
            dp=ndp
        tot+=dp[k] if k<=12 else 0
    return tot//4
# brute orbit check for small k
def canon(subset):
    return min(tuple(sorted(idx[s(cells[i])] for i in subset)) for s in syms)
def brute(k):
    return len({canon(comb) for comb in combinations(range(12),k)})
seq_b=[burnside(k) for k in range(13)]
seq_brute=[brute(k) for k in range(7)]
print("Burnside k=0..12 :", seq_b)
print("brute   k=0..6   :", seq_brute, "(must match)")
print("k=0..3 ->", seq_b[:4], "= Fib∩triangular {1,3,21,55};  k=4 ->", seq_b[4],
      "(135, neither). break confirmed.")

# ---------------------------------------------------------------- T5
line("T5  Fibonacci rank of apparition alpha(p) | p - (5/p)  (Euler criterion)")
def legendre(a,p):
    a%=p
    if a==0: return 0
    ls=pow(a,(p-1)//2,p)
    return -1 if ls==p-1 else 1
def fib_mod(k,p):
    a,b=0,1
    for _ in range(k): a,b=b,(a+b)%p
    return a
def alpha(p):
    k=1
    while True:
        if fib_mod(k,p)%p==0: return k
        k+=1
def primes_upto(N):
    sieve=[True]*(N+1); sieve[0:2]=[False,False]
    for i in range(2,isqrt(N)+1):
        if sieve[i]:
            for j in range(i*i,N+1,i): sieve[j]=False
    return [i for i in range(2,N+1) if sieve[i]]
print("alpha(2)=%d alpha(3)=%d alpha(5)=%d alpha(7)=%d alpha(11)=%d  (expect 3,4,5,8,10)"
      % (alpha(2),alpha(3),alpha(5),alpha(7),alpha(11)))
print("  p  alpha  (5/p)  p-(5/p)  divides? mod5")
law=True
for p in primes_upto(60):
    if p==5: continue
    al=alpha(p); l=legendre(5,p); tgt=p-l
    div=(tgt%al==0); law&=div
    print(f"{p:3d}  {al:4d}  {l:+d}    {tgt:4d}     {div}    {p%5}")
print("LAW alpha(p) | p-(5/p) holds for all p<60:", law,
      "   (5 QR <=> p ≡ ±1 mod 5 <=> alpha|p-1)")

# ---------------------------------------------------------------- T6
line("T6  Heron square = (x+y+z)*x*y*z;  Lagrange 4-sq; Gauss 3-tri (EYPHKA)")
print("Ravi x=s-a,y=s-b,z=s-c -> s=x+y+z, Area^2=(x+y+z)*x*y*z for the 5 m=1 triangles:")
for a,b,c,A,P in five:
    s=(a+b+c)//2; x,y,z=s-a,s-b,s-c
    print(f"  ({a},{b},{c}) s={s} (x,y,z)=({x},{y},{z}) (x+y+z)xyz={(x+y+z)*x*y*z} =? A^2={A*A}  {(x+y+z)*x*y*z==A*A}")
def reps_4sq(N):
    c=0
    r=isqrt(N)
    for w in range(r+1):
        for x in range(w,isqrt(N-w*w)+1):
            rem=N-w*w-x*x
            for y in range(x,isqrt(rem)+1):
                z2=rem-y*y
                zz=isqrt(z2)
                if zz>=y and zz*zz==z2: c+=1
    return c
def reps_3tri(N):
    c=0; tris=[]; k=0
    while T(k)<=N: tris.append(T(k)); k+=1
    for i,ti in enumerate(tris):
        for j in range(i,len(tris)):
            rem=N-ti-tris[j]
            if rem<0: break
            if is_tri(rem) and rem>=tris[j]: c+=1
    return c
miss4=[N for N in range(0,200) if reps_4sq(N)==0]
miss3=[N for N in range(0,200) if reps_3tri(N)==0]
print("naturals 0..199 NOT a sum of 4 squares (Lagrange):", miss4, "(empty => theorem holds)")
print("naturals 0..199 NOT a sum of 3 triangular (Gauss): ", miss3, "(empty => EYPHKA holds)")
print("\nThe recursive duality: square = PRODUCT of 4 naturals (Heron, glue=sum),")
print("natural = SUM of 4 squares (Lagrange); '4'=2^2 (parity squared), '3'=triangle sides.")
print("\nDONE.")
