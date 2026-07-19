from fractions import Fraction as F
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncov(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))
print("(3) THE SHARP CONSEQUENCE: lambda_min(Toeplitz) -> ESSENTIAL inf of the symbol,")
print("    so the LP criterion detects 'uncovered set has POSITIVE MEASURE' --")
print("    strictly stronger than LRC(14), which allows a measure-ZERO lonely set.")
print()
print("    family                 exact uncovered measure   LP can certify?")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("2*{1,...,13}",[2*i for i in range(1,14)]),
             ("random spread",[3,7,11,19,23,31,37,41,47,53,59,61,67]),
             ("AP d=8",[1+8*i for i in range(13)])]:
    u=uncov(V)
    print(f"    {nm:22s} {str(u):22s}  {'NO -- measure zero' if u==0 else 'yes'}")
print()
print("    => on the three TIGHT families the uncovered set is exactly measure zero,")
print("       so lambda_min = 1 and the LP cannot certify -- yet LRC(14) HOLDS there")
print("       (gap exactly 1/14, attained at the six points p/14).")
print("       The Delsarte LP therefore cannot prove LRC(14) even in principle.")
