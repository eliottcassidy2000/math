from collections import Counter

def kids(x):
    return [2*x]+([(2*x-1)//3] if x%3==2 else [])

def trace(y):
    j=1
    sterile=[]
    while True:
        children=kids(y)
        live=[z for z in children if z%3]
        if len(live)==2:return j,tuple(sterile)
        if len(children)==2:sterile.append(j)
        y=live[0];j+=1

D={}
for x in range(3,1000):
    live=[y for y in kids(x) if y%3]
    if len(live)!=2:continue
    T=tuple(sorted(trace(y) for y in live))
    D.setdefault(x%9,Counter())[T]+=1
assert set(D[2])=={((2,()),(2,())),((2,()),(4,(2,)))}
assert set(D[8])=={((1,()),(4,(2,))),((3,(1,)),(4,(2,)))}
print('weighted skeleton and all sterile-arm rules independently checked:',D)
