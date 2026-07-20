from math import comb, gcd, lcm
def pisano(m):
    a,b=0,1
    for i in range(1,m*m*6+10):
        a,b=b,(a+b)%m
        if a==0 and b==1: return i
print("(1) PISANO PERIOD of 10 -- the user's claim: last digit repeats every 60")
print(f"    pi(10) = {pisano(10)}   -> CLAIM IS CORRECT")
print(f"    structure: pi(2)={pisano(2)}, pi(5)={pisano(5)}, lcm = {lcm(pisano(2),pisano(5))}")
print()
print("(2) '1001 = three sixties' -- does this check out arithmetically?")
print(f"    3 x 60 = 180.  1001 != 180.  1001/60 = {1001/60:.4f}.  NOT three sixties.")
print(f"    1001 factors: 7 x 11 x 13 = {7*11*13}")
print(f"    Pisano period of 1001 = lcm(pi(7),pi(11),pi(13)) = lcm({pisano(7)},{pisano(11)},{pisano(13)}) = {lcm(pisano(7),pisano(11),pisano(13))}")
print(f"    (direct: pi(1001) = {pisano(1001)})   -- so 1001's Pisano period is 560, not 180")
print()
print("(3) BUT: 1001 IS A PASCAL ENTRY, AND IN ROW 14.")
for k in range(0,15):
    if comb(14,k)==1001: print(f"    C(14,{k}) = {comb(14,k)}")
print(f"    row 14 of Pascal: {[comb(14,k) for k in range(15)]}")
print("    -> 1001 = C(14,4) = C(14,10).  14 is THIS REPO'S APEX (LRC(14), lam=1/14).")
print()
print("(4) FIBONACCI FROM SHIFTED PASCAL (shallow diagonals) -- verify")
def fib(n):
    a,b=0,1
    for _ in range(n): a,b=b,a+b
    return a
ok=all(sum(comb(n-1-k,k) for k in range(0,(n)//2+1))==fib(n) for n in range(1,40))
print(f"    F(n) = sum_k C(n-1-k, k) for n=1..39: {ok}")
print(f"    e.g. n=10: {[comb(9-k,k) for k in range(5)]} sums to {sum(comb(9-k,k) for k in range(5))} = F(10) = {fib(10)}")
print()
print("(5) WHY PERIODICITY TRANSFERS: Lucas' theorem makes Pascal mod p self-similar,")
print("    so the shallow-diagonal sums inherit a period.  Check ord and periods:")
for p in [2,5,7,11,13]:
    print(f"    p={p:3d}: pi(p)={pisano(p):4d}   p-1={p-1:3d}   p+1={p+1:3d}   "
          f"pi | p-1? {(p-1)%pisano(p)==0}   pi | 2(p+1)? {(2*(p+1))%pisano(p)==0}")
print()
print("(6) Does 14 / 1001 / 60 cohere? sanity checks:")
print(f"    C(14,4)=1001=7*11*13 ; 14=2*7 ; ord_1001(10)={next(d for d in range(1,20) if pow(10,d,1001)==1)}  (10^3 = -1 mod 1001)")
print(f"    the repo's two extremal families are {{1..13}} and {{1..11,13,24}} -- contain 7,11,13")
