from math import gcd

def cycle(p, y):
    return [(2**k*y)%(2*p*p) for k in range(1,p*(p-1)+1)]
for p in (3,5):
    outputs=[]
    for y in range(1,2*p*p):
        if y%p==0: continue
        states=cycle(p,y)
        assert set(states)=={x for x in range(2,2*p*p,2) if gcd(p,x)==1}
        exits=[(z-1)//p%p for z in states if z%(2*p)==p+1]
        assert sorted(exits)==list(range(p))
        outputs.append(sum(x!=0 for x in exits))
    print("prime",p,"clock",p*(p-1),"live",sorted(set(outputs)))

def forward(v):
    return 3*v+1 if v%2 else v//2

level={8}
for m in range(1,9):
    nxt=set()
    for y in level:
        outcomes=[(2**k*y-1)//3 for k in range(1,7)
                  if (2**k*y)%18 in (4,16)]
        assert len(outcomes)==2 and not nxt.intersection(outcomes)
        for v in outcomes:
            a=v
            for _ in range(7):
                a=forward(a)
                if a==y:break
            assert a==y and v%3
        nxt.update(outcomes)
    assert len(nxt)==2**m
    level=nxt
print("first eight macro levels independently replayed")
