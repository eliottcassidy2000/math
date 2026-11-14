# Root-change corridors: J acts on the root triple (H, E, lam) by [[1,-1,1],[1,-1,0],[1,0,1]] (char. poly x^3-x^2-x-1).
# (1) The (H,lam)-block of J^k is rank one exactly at k = 1 and 14 (k <= 400): J^14 block = (47,56)^T (47,56).
# (2) After an O letter (E = 0), J^k forces gcd(J^k[0][0], J^k[0][2]) | H.  (3) Check on random root-change words.
from math import gcd
import random
def matmul(A,B): return [[sum(A[i][t]*B[t][j] for t in range(3)) for j in range(3)] for i in range(3)]
J=[[1,-1,1],[1,-1,0],[1,0,1]]
P=[[1,0,0],[0,1,0],[0,0,1]]; rank1=[]; forced=[]
for k in range(1,401):
    P=matmul(P,J)
    if P[0][0]*P[2][2]-P[0][2]*P[2][0]==0: rank1.append(k)
    g=gcd(P[0][0],P[0][2])
    if g>1 and k<=60: forced.append((k,g))
print("rank-one (H,lam)-blocks of J^k, k<=400:",rank1)
print("forced factors of H after O J^k, k<=60:",forced)
def step(s,m):
    H,E,L,h=s
    if m=='a': return (H,0,2*L,h)
    if m=='b': return (H+L-E,h,2*L-h,h)
    return (H+L-E,H-E,H+L,H)
def run(w,s=(1,1,1,1)):
    for m in w: s=step(s,m)
    return s
random.seed(7); ok=0
for _ in range(1000):
    x=''.join(random.choice('abc') for _ in range(random.randint(1,14)))+'a'
    s=run(x); t=run('c'*14,s); m=47*s[0]+56*s[2]
    ok+= (t[0]==47*m and t[2]==56*m)
print("words x O J^14 with H = 47(47H(x O)+56 lam(x O)) and lam = 56(...):",ok,"of 1000")
# (4) det of the (H,lam)-block of J^k equals -T(-k-2), T the two-sided tribonacci with T(0)=T(1)=0, T(2)=1
#     (wave 16's corridor determinant det K_q = -T(-q-3) at q = k-1); its zeros k=1,14 come from T(-3)=T(-16)=0.
T={0:0,1:0,2:1}
for n in range(-1,-500,-1): T[n]=T[n+3]-T[n+2]-T[n+1]
P=[[1,0,0],[0,1,0],[0,0,1]]; bad=0
for k in range(1,301):
    P=matmul(P,J)
    if P[0][0]*P[2][2]-P[0][2]*P[2][0] != -T[-k-2]: bad+=1
print("det(H,lam-block of J^k) == -T(-k-2) for k=1..300: mismatches",bad)
