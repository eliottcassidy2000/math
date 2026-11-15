from replay_target_carry import branch,p

def need(t,*details):
    if not t: raise AssertionError(details)

def choose(s):
    D=64
    lo,hi,c,_=branch(s,D)
    for b in range(lo,hi):
        t=c-b; free=b^t
        for j in range(1,4):
            if free>>j & 1 and sum(free>>k & 1
                    for k in range(max(0,j-2),j+3) if k!=j)%2:
                return b,t,j
    raise AssertionError(s)

for s in range(64):
    b,t,j=choose(s)
    for u in (0,1,2,7,19,93849):
        for v in (0,1,2,15):
            x,y=64*u+b,64*v+t
            X,Y=x^(1<<j),y^(1<<j)
            need(X+Y==x+y and X^Y==x^y and X%2==x%2,s,u,v)
            need(((p(X)+p(Y)-p(x)-p(y))//2)%2==1,s,u,v)
print("64 local target carry cases; 1536 literal high-prefix toggle controls")
