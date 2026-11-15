from random import Random
def h(w):
    s=1
    while w:
        w,b=divmod(w,2)
        if not b:s+=w+(w&1)
    return s
def need(x):
    if not x:raise RuntimeError("odd contrast")
rng=Random(20261115)
for s in range(9,47):
    m=(1<<s)-3
    A={0:1,2*s-3:-1,2*s:-1}; C=(0,4,6)
    low={p+q:-z for p,z in A.items() for q in C}
    top={p+4*s:z for p,z in A.items()}
    need(len(low)==9 and len(top)==3 and not low.keys()&top.keys())
    coef=low|top
    old=sum(1<<i for i,z in coef.items() if z==-1)
    mask=sum(1<<i for i in coef)
    width=6*s+1
    def guarded(w):return (1<<(width+1))+(w<<1)
    def contrast(w):return h(guarded(w^mask))-h(guarded(w))
    base=contrast(old)
    # The formula for unchanged bits depends only on distances;
    # adding a common guard does not change it.
    j=4*s-2
    c=(pow(2,-j-1,m)*sum(z*pow(2,i,m) for i,z in low.items()))%m
    need(sum(z*pow(2,i,m) for i,z in coef.items())%m==0)
    need(sum(z*pow(2,-i,m) for i,z in coef.items())%m==0)
    for _ in range(40):
        target=rng.randrange(m)
        binary=((target-base)*pow(c,-1,m))%m
        ext=sum((binary>>r&1)<<(j-r) for r in range(s))
        need(not ext&mask and contrast(old|ext)%m==target)
print("1520 exact contrasts at moving moduli, through s=46")
