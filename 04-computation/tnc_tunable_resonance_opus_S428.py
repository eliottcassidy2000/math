import sympy as sp
from itertools import combinations_with_replacement as cwr
u,a=sp.symbols('u a')
def CT(R,N,mm):
    return sp.Poly(sp.expand(R**mm),u).coeff_monomial(u**(N*mm))
def charges(N,c,M):  # charge set {-N, c, M}
    return [-N,c,M]
def min_reps(chs):
    for m in range(1,20):
        reps=[r for r in cwr(chs,m) if sum(r)==0 and any(x!=chs[1] for x in r)]  # use middle nontrivially? no, count all
        reps=[r for r in cwr(chs,m) if sum(r)==0]
        if reps: return m,len(reps),reps
    return None,0,[]
print("ANGULAR UNIFORM: CHARACTERISE the TUNABLE charge-triples {-N, c, M} (non-unique minimal")
print("rep of 0). Unique-minimal => THM-1655 (positivity, no resultant). Tunable => the thin")
print("subset needing the finite-place certificate. Find the arithmetic law.")
print()
print("   N   c   M   min-m  #reps  tunable?   condition guess")
tun=[]
for N in range(1,6):
    for M in range(1,8):
        for c in range(-N+1,M):
            if c==0: continue
            chs=[-N,c,M]
            m0,nr,reps=min_reps(chs)
            if m0 is None: continue
            tunable = nr>=2
            if tunable:
                tun.append((N,c,M,m0,reps))
                print(f"   {N:2d} {c:3d} {M:2d}   {m0:4d}   {nr:3d}    YES       reps={reps}")
print()
print(f"tunable triples found (N<=5, M<=7): {len(tun)}")
print()
print("LOOK FOR THE LAW: in each tunable triple, is there a linear relation among -N,c,M")
print("(a 'sub-resonance') that creates the second rep?")
for N,c,M,m0,reps in tun:
    # the two reps differ by a vector in the kernel; find the primitive relation
    from sympy import Matrix, gcd
    # relation: alpha*(-N)+beta*c+gamma*M=0 with small alpha,beta,gamma
    rel=None
    for al in range(-4,5):
        for be in range(-4,5):
            for ga in range(-4,5):
                if (al,be,ga)==(0,0,0): continue
                if al*(-N)+be*c+ga*M==0 and abs(al)+abs(be)+abs(ga)<=6:
                    rel=(al,be,ga); break
            if rel: break
        if rel: break
    print(f"   {{-{N},{c},{M}}}: primitive relation {rel[0]}*(-{N})+{rel[1]}*({c})+{rel[2]}*{M}=0")
