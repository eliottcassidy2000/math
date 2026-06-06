def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
def sigma_proper(n): return sum(d for d in range(1,n) if n%d==0)
print("=== PERFECT NUMBERS and the arc: 6 is the first perfect number ===")
perfects=[n for n in range(2,9000) if sigma_proper(n)==n]
print(f"  perfect numbers < 9000: {perfects}  (6 = 1+2+3, the first)")
print(f"  Euclid-Euler: even perfect = 2^(p-1)*(2^p-1), 2^p-1 a Mersenne prime:")
for p in [2,3,5,7]:
    M=2**p-1
    print(f"    p={p}: M_p=2^{p}-1={M} {'(Mersenne PRIME)' if isprime(M) else '(composite)'}, perfect = 2^{p-1}*{M} = {2**(p-1)*M}")

print("\n=== the 6 = 2*3 meeting point IS the first perfect number ===")
print(f"  6 = 2*3 = 2*M_2 (M_2=3, first Mersenne prime) = doubled prime (2*3) = tripled prime (3*2, S630) = FIRST PERFECT.")
print(f"  6 = the hexagonal / pi-3 / dZ=1/6 / lcm(2,3) -- and now the first perfect number. Everything meets at 6.")

print("\n=== the forbidden H-values are MERSENNE products: 7 = M_3 = Phi_3(2), 21 = M_2*M_3 ===")
print(f"  7 = 2^3-1 = M_3 (2nd Mersenne prime) = Phi_3(2) (S628, cube-root) = the forbidden H atom")
print(f"  21 = 3*7 = M_2 * M_3 = the other forbidden H. So forbidden H {{7,21}} = {{M_3, M_2*M_3}} (Mersenne monoid elements)")
print(f"  28 = 2^2 * 7 = 2^2 * M_3 = the 2nd perfect number, built from the forbidden-H Mersenne prime 7!")

print("\n=== PARITY exclusion: even perfect numbers are EVEN => never an H-value (Redei: H odd) ===")
print(f"  6,28,496,8128 all even => not in the H-spectrum (which is a sub-monoid of ODD numbers).")
print(f"  ODD perfect numbers: NONE known (famous open). If one exists & odd, is it in the H-monoid? -> a bridge to the converse-of-Redei.")

print("\n=== the H-spectrum = multiplicative MONOID of strong (strongly-connected) H-values ===")
import itertools
def make_beats(n,bits):
    b=[[False]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: b[i][j]=True
            else: b[j][i]=True
            idx+=1
    return b
def ham(n,beats):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1<<n):
        for last in range(n):
            c=dp[m][last]
            if not c: continue
            for nx in range(n):
                if not (m>>nx)&1 and beats[last][nx]: dp[m|(1<<nx)][nx]+=c
    return sum(dp[full][v] for v in range(n))
def strong(n,beats):
    reach=[[i==j or beats[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]: reach[i][j]=True
    return all(reach[i][j] and reach[j][i] for i in range(n) for j in range(n))
atoms={1}
for n in range(2,7):
    m=n*(n-1)//2
    for bits in range(1<<m):
        b=make_beats(n,bits)
        if strong(n,b): atoms.add(ham(n,b))
atoms=sorted(a for a in atoms if a<=50)
# monoid closure
mono={1}; ch=True
while ch:
    ch=False
    for v in list(mono):
        for a in atoms:
            if v*a<=50 and v*a not in mono: mono.add(v*a); ch=True
forb=[h for h in range(1,51,2) if h not in mono]
print(f"  strong (atom) H-values <=50: {atoms}")
print(f"  multiplicative monoid <=50: {sorted(mono)}")
print(f"  FORBIDDEN odd (not in monoid) <=50: {forb}  -> {{7, 21}} = the permanent gaps = {{M_3, M_2 M_3}}")
