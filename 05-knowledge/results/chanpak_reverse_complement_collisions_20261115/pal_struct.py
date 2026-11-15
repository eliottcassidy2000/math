from collections import Counter
def p(u):
    if u<=1: return 0
    return p(u>>1) if u&1 else p(u>>1)+(u>>1)+((u>>1)&1)
def rev(x,a):
    r=0
    for i in range(a): r=(r<<1)|((x>>i)&1)
    return r
for a in range(3,13):
    pairs=[]
    for x in range(1<<(a-1),1<<a):
        for y in range(0,1<<(a-1)):
            if (x^y)&1==0: continue
            if x-y==rev(x,a)-rev(y,a): pairs.append((x,y))
    diffs=Counter((p(x)-p(y))//(x-y) if (p(x)-p(y))%(x-y)==0 else 'frac' for x,y in pairs)
    eq=[(x,y) for x,y in pairs if p(x)==p(y)]
    # structure of equal-p pairs: is y = x - (something)?  check if x,y share middle bits pattern
    sym=Counter()
    for x,y in eq:
        xs=format(x,'0%db'%a); ys=format(y,'0%db'%a)
        sym[(xs[0]+xs[-1], ys[0]+ys[-1])]+=1
    print('a=%2d pal pairs %5d  (p(x)-p(y))/a distribution %s  equal-p %d'%(a,len(pairs),dict(sorted(diffs.items(),key=str)),len(eq)))
    if a<=6:
        for x,y in eq: print('     x=%s y=%s'%(format(x,'0%db'%a),format(y,'0%db'%a)))
