"""
Upper bound via NEAREST-NEIGHBORS of 0 (three-gap). For the runners {vt:v=1..n-1}:
  d+ = nearest runner above 0 (index p), d- = nearest below (index p').
  merged gap G = d+ + d- (the runner-gap containing 0).
TEST relations: (a) p+p' >= n (three-gap)? (b) at optimal t, ||p't||=(1-2c)/n? (c) the bound reduces to
  ||p't|| <= (1-2c)/n  when c is in the merged gap centered.
Also: is c ALWAYS in the merged gap (for the MAXIMIZING t)? If yes, the proof localizes to 0's neighborhood.
"""
from fractions import Fraction
def frac(x): x=x%1; return min(x,1-x)
def nn0(t,n):
    # nearest runner above and below 0
    above=[]; below=[]
    for v in range(1,n):
        x=(Fraction(v)*t)%1
        if x==0: return None
        if x<=Fraction(1,2): above.append((x,v))
        else: below.append((1-x,v))
    dp,p=min(above) if above else (Fraction(1),0)
    dm,pp=min(below) if below else (Fraction(1),0)
    return dp,p,dm,pp
n=9; AP=list(range(1,n)); Qmax=7*n
print(f"n={n}: at the MAXIMIZING t for each c -- nearest-nbr structure of 0:")
print(f"  {'c':>6} {'M_c':>8} {'t*':>8} {'p(above)':>8} {'p~(below)':>9} {'p+p~':>5} {'||p~t||':>8} {'(1-2c)/n':>9} {'c in merged?':>12}")
for cn,cd in [(1,18),(1,9),(1,6),(2,9),(1,4),(5,18),(1,3)]:
    c=Fraction(cn,cd)
    best=Fraction(-1); bt=None
    for q in range(1,Qmax+1):
        for a in range(q):
            t=Fraction(a,q); m=min(frac(Fraction(v)*t-c) for v in AP)
            if m>best: best=m; bt=t
    r=nn0(bt,n)
    if r is None: print(f"  c={c}: t* hits a runner at 0"); continue
    dp,p,dm,pp=r
    # is c in the merged gap (-dm, dp)?
    cmod=c%1; inmerged = (cmod<dp) or (cmod>1-dm)
    print(f"  {str(c):>6} {str(best):>8} {str(bt):>8} {p:>8} {pp:>9} {p+pp:>5} {str(dm):>8} {str(Fraction(1,n)-2*c/n):>9} {str(inmerged):>12}")
print()
print("Now NON-optimal t: does p+p'>=n ALWAYS hold (three-gap), and is min_v||vt-c|| <= G/2 = (d++d-)/2?")
import random; random.seed(5); cnt=0; ppmin=99; badG=0
for _ in range(4000):
    q=random.randint(3,150); a=random.randint(1,q-1); t=Fraction(a,q)
    r=nn0(t,n)
    if r is None: continue
    dp,p,dm,pp=r; cnt+=1; ppmin=min(ppmin,p+pp)
    c=Fraction(random.randint(0,n),2*n)
    m=min(frac(Fraction(v)*t-c) for v in AP)
    if m > (dp+dm)/2 + Fraction(1,10**9): badG+=1
print(f"  over {cnt} t: min(p+p') = {ppmin} (>=n={n}? {ppmin>=n}); min_v||vt-c|| > G/2 count = {badG}")
print("  => if p+p'>=n always, the merged gap G=||pt||+||p't|| with p+p'>=n is the three-gap LARGEST-gap key.")
