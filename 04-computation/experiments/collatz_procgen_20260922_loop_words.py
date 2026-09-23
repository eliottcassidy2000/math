# find minimal-K reverse E-loops at 1 for s<=22 with explicit words (BFS with parent pointers, values <= XMAX)
import math
XMAX=3_000_000
SMAX=22
layer={1:(0,None)}  # value -> (K, (prev_value, k))
hist=[layer]
for s in range(1,SMAX+1):
    nxt={}
    for x,(K,_) in layer.items():
        if x%3==0: continue
        k0=0 if x%3==1 else 1
        k=k0
        while True:
            t=(2**k)*x
            if t>3*XMAX+1: break
            y=(t-1)//3
            if y>=1 and y%3!=0:
                if y not in nxt or K+k<nxt[y][0]: nxt[y]=(K+k,(x,k))
            k+=2
    layer=nxt; hist.append(layer)
    if 1 in layer:
        # backtrack
        word=[]; v=1; ss=s
        while ss>0:
            K,(pv,k)=hist[ss][v]; word.append((k,v)); v=pv; ss-=1
        word.reverse()
        K=layer[1][0]
        vals=[1]+[w[1] for w in word]
        print(f"s={s:2d} K={K:2d} ratio={2**K/3**s:.3f} ks={[w[0] for w in word]} values={vals}")
