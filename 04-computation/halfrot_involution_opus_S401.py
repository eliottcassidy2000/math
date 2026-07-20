from fractions import Fraction as F
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
print("THE HALF-ROTATION s: t -> t+1/2 is the FREE Z/2 on R/Z (BU's antipodal map).")
print("The reflection t -> -t FIXES 0 and 1/2, so it is NOT BU's antipodal map.\n")
print("(a) How does ||vt|| transform under s?  v even: unchanged.  v odd: ||vt|| -> 1/2-||vt||")
bad=0
for _ in range(4000):
    v=random.randint(1,60); t=F(random.randint(0,997),998)
    L=fn((t+F(1,2))*v)
    R= fn(t*v) if v%2==0 else F(1,2)-fn(t*v)
    if L!=R: bad+=1
print(f"    checked 4000 random (v,t): mismatches = {bad}   -> law holds exactly\n")
print("(b) CONSEQUENCE: for v ODD, f_v(t) := ||vt|| - 1/4 is ODD under s:")
print("       f_v(t+1/2) = (1/2-||vt||) - 1/4 = 1/4-||vt|| = -f_v(t).")
bad=0
for _ in range(4000):
    v=2*random.randint(0,30)+1; t=F(random.randint(0,997),998)
    if fn((t+F(1,2))*v)-F(1,4) != -(fn(t*v)-F(1,4)): bad+=1
print(f"    checked 4000 random odd v: mismatches = {bad}   -> genuinely ODD equivariant\n")
print("(c) Is g_V itself s-invariant?  Only if ALL speeds are even (impossible for primitive V).")
for nm,V in [("{1,...,13}",list(range(1,14))),("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    diffs=[(t,fn_g) for t in [F(p,14) for p in range(1,14)]
           for fn_g in [ (min(fn(t*v) for v in V), min(fn((t+F(1,2))*v) for v in V)) ]]
    viol=[(str(t),str(a),str(b)) for t,(a,b) in diffs if a!=b]
    print(f"    {nm}: g(t) vs g(t+1/2) differ at {len(viol)}/13 sampled points -> NOT s-invariant")
print()
print("  => The free Z/2 that BU needs is s (half-rotation), and g is NOT s-invariant,")
print("     so Argmax is not an s-set and s does not organise the witnesses.")
print("     The reflection DOES organise them but is not free on the circle.")
print("     BU's two hypotheses (free action + odd map) are satisfied by DIFFERENT")
print("     involutions here.  That mismatch -- not evenness -- is the real obstruction.")
