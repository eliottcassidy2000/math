import itertools
from collections import Counter
def H2(x,y,e):
    # same DP as twohub.c
    st={}
    for o in (0,1): st[(o,0,0,0)]=1
    for xi,yi in zip(x,y):
        nt={}
        for (o,s0,s1,s2),c in st.items():
            h1b= xi if o==0 else yi; h2b= yi if o==0 else xi
            k=(o,2 if h1b else 1,s1,s2); nt[k]=nt.get(k,0)+c
            if s1 or h1b: k=(o,s0,2 if h2b else 1,s2); nt[k]=nt.get(k,0)+c
            if s2 or h2b: k=(o,s0,s1,1); nt[k]=nt.get(k,0)+c
        st=nt
    H=0
    for (o,s0,s1,s2),c in st.items():
        h12= e if o==0 else 1-e
        if s0==2 or s1==2 or (s1==0 and not h12): continue
        H+=c
    return H
L=7; stat=Counter(); stat_pre=Counter()
for m in range(0,L+1):
    for x in itertools.product((0,1),repeat=m):
        for y in itertools.product((0,1),repeat=m):
            for e in (0,1):
                H=H2(x,y,e)
                app={H2(x+(a,),y+(b,),f) for a in (0,1) for b in (0,1) for f in (0,1)}
                pre={H2((a,)+x,(b,)+y,f) for a in (0,1) for b in (0,1) for f in (0,1)}
                key=tuple(int(3*H+d in app) for d in (0,2,4)); stat[key]+=1
                key2=tuple(int(3*H+d in pre) for d in (0,2,4)); stat_pre[key2]+=1
tot=sum(stat.values())
print("append one spine vertex: which of (3H, 3H+2, 3H+4) occur, over",tot,"configs:",dict(stat))
print("prepend one spine vertex:",dict(stat_pre))
