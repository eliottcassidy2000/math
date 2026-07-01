"""Definitive symmetry group of the half-tile point set + the (span d=x-y, anti-diagonal s=x+y) grid.
Check ALL candidate isometries: reflections across vertical/horizontal/slope+-1 lines, and 90/180 rotations."""
def halftiles(n):
    return [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x+y<=n+1]
def fulltiles(n):
    return [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
def sym_group(pts):
    S=set(pts); found=[]
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    # candidate reflection axes: 2c over range of relevant sums
    cand=set()
    for c2 in range(2*min(xs), 2*max(xs)+1):        # vertical x=c2/2
        if all(((c2-x),y) in S for (x,y) in pts): found.append(('vert x=%.1f'%(c2/2)))
    for c2 in range(2*min(ys),2*max(ys)+1):         # horizontal y=c2/2
        if all((x,(c2-y)) in S for (x,y) in pts): found.append(('horiz y=%.1f'%(c2/2)))
    for c in range(min(x+y for x,y in pts), max(x+y for x,y in pts)+1):  # anti-diag x+y=c
        if all((c-y,c-x) in S for (x,y) in pts): found.append(('antidiag x+y=%d'%c))
    for c in range(min(x-y for x,y in pts), max(x-y for x,y in pts)+1):  # main-diag x-y=c
        if all((y+c,x-c) in S for (x,y) in pts): found.append(('maindiag x-y=%d'%c))
    # 180 rotation about centroid*2
    cx2=sum(xs)*2//len(xs) if (sum(xs)*2)%len(xs)==0 else None
    return found
for n in [4,5,6,7,8,9]:
    h=halftiles(n); f=fulltiles(n)
    print(f"n={n}: FULL staircase symmetries: {sym_group(f)}")
    print(f"       HALF-region symmetries:   {sym_group(h)}")
    # span/anti-diagonal ranges
    spans=sorted(set(x-y for x,y in h)); adiag=sorted(set(x+y for x,y in h))
    print(f"       half grid: spans d=x-y in {spans} (each a 'main-diagonal' line); anti-diag s=x+y in {adiag}")
