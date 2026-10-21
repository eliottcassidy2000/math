import json
from itertools import product

def old(i,a,b,c):
    if b==0:
        if c==0:return None
        if c==1:return (1,2,2,a+1) if i==0 else (0,1,0,a+3)
        return (0,1,a+3,c-2) if i==0 else (1,2,a+3,c-2)
    if c==0:return (1,a+1,0,b-1) if i==0 else (0,a+2,0,b-1)
    return (i,a+3,b-1,c-1)
def enc(s):
    if s is None:return None
    i,a,b,c=s;return (2*a+i,b,c)
def new(a,b,c):
    if b==0:
        if c==0:return None
        k,i=divmod(a,2)
        if c==1:return (5,2,k+1) if i==0 else (2,0,k+3)
        return (2,k+3,c-2) if i==0 else (5,k+3,c-2)
    if c==0:return (a+3,0,b-1)
    return (a+6,b-1,c-1)
count=0
for s in product(range(2),range(31),range(31),range(31)):
    assert enc(old(*s))==new(*enc(s));count+=1
# Matrix checks using exact integers.
def mul(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
I=((1,0),(0,1));B=((0,-1),(1,1));J=((1,1),(0,-1))
powers=[I]
for k in range(6): powers.append(mul(powers[-1],B))
assert powers[3]==((-1,0),(0,-1)) and powers[6]==I
assert mul(J,J)==I and mul(mul(J,B),J)==powers[5]
G=set(powers[:6])|set(mul(g,J) for g in powers[:6])
assert len(G)==12
assert sum(mul(g,B)==mul(B,g) for g in G)==6
out={'A9_reparameterization_cases':count,'A9_phase_forgetting_counterexample':{'input_phase0':[0,0,0,1],'input_phase1':[1,0,0,1],'output_phase0':old(0,0,0,1),'output_phase1':old(1,0,0,1)},'six_vector_group_order':len(G),'centralizer_of_B_in_group':6,'scope':'Finite controls of algebraic proofs; not a raw Turing-machine audit or proof of halting.'}
print(json.dumps(out,indent=2))
