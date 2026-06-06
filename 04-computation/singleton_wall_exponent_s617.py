from fractions import Fraction as F
import math

def p0_grid(speeds, delta, res):
    # speeds: list of Fraction; p0 = fraction of grid t with all ||s t|| >= delta  (lonely measure)
    cnt=0
    for j in range(res):
        t=F(j,res); ok=True
        for s in speeds:
            x=(s*t)%1
            if min(x,1-x)<delta: ok=False; break
        if ok: cnt+=1
    return cnt/res

# Perturb the TOP speed of a rank-r collapse wall by epsilon and measure how the lonely measure p0(eps) opens.
# Recursive-log prediction: exponent of the opening tracks the relation RANK (one averaging level per relation).
def opening_exponent(base, top_index, label, rank):
    n=len(base); delta=F(1,n+1)
    print(f"\n{label}: base wall {base}, rank(additive relations)={rank}, delta={delta}")
    print(f"   {'eps':>10} {'p0(eps)':>12} {'p0/eps':>10} {'p0/eps^2':>10}")
    eps_list=[F(1,d) for d in [50,100,200,400,800,1600]]
    xs=[]; ys=[]
    for eps in eps_list:
        sp=[F(b) for b in base]; sp[top_index]=F(base[top_index])+eps
        res= max(60000, int(8/float(eps)))   # resolution finer than eps
        p0=p0_grid(sp, delta, res)
        if p0>0:
            xs.append(math.log(float(eps))); ys.append(math.log(p0))
        print(f"   {float(eps):>10.5f} {p0:>12.6f} {p0/float(eps):>10.3f} {p0/float(eps)**2:>10.1f}")
    if len(xs)>=3:
        k=len(xs); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(x*y for x,y in zip(xs,ys))
        a=(k*sxy-sx*sy)/(k*sxx-sx*sx)
        print(f"   --> fitted opening exponent  p0 ~ eps^a   a = {a:.3f}")

# rank-1 singleton wall: {1,2} (single relation 1+1=2)
opening_exponent([1,2], 1, "SINGLETON WALL {1,2}", 1)
# rank-2 sporadic chain {1,3,4,7} -- perturb top 7 (breaks the relation 3+4=7)
opening_exponent([1,3,4,7], 3, "SPORADIC CHAIN {1,3,4,7}", 2)
# AP {1,2,3,4} -- many relations, perturb top
opening_exponent([1,2,3,4], 3, "AP WALL {1,2,3,4}", 4)
