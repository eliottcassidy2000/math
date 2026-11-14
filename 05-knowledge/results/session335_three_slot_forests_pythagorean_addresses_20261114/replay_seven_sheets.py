from collections import Counter
from replay_forest import orders

coefs=[(-5,24,18),(-9,25,20),(-9,5,40),(-11,24,24)]
for d,(a,b,c) in enumerate(coefs):
    vals=[]
    for k in (1,2,3):
        w=[1,3]+[d]*k+[0,0]
        ps=orders(w)
        assert len(ps)==a+b*2**k+c*3**k
        chosen=Counter()
        for p in ps:
            runs=[[]]
            for v in p:
                if v<len(w):runs[-1].append(v)
                else:runs.append([])
            q=tuple(sorted(tuple(g) for g in runs if g))
            if len(q)==2 and {f[0] for f in q}=={0,1} and {f[-1] for f in q}=={k+2,k+3}:
                chosen[q]+=1
        assert len(chosen)==2**(k+1) and set(chosen.values())=={7}
        vals.append(len(ps))
    print(d,vals)
