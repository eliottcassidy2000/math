from random import Random
B,Q,T,S,h=2855,5760,4034,3914,45
codes=((1<<12,2787,2875,12),(1<<36,46767159475,48234217403,36))
def check(x,*d):
    if not x:raise AssertionError(d)
def p(x):
    y=0
    while x:
        x,b=divmod(x,2)
        if not b:y+=x+(x&1)
    return y
def F(x,y):
    return B+Q*x*y+T*x+S*y+h*(p(x)+p(y))
def rev(x,r):
    return int(f'{x:0{r}b}'[::-1],2)
def suffix(x,side):
    present=[i for i,code in enumerate(codes) if x%code[0]==code[side]]
    check(len(present)<=1)
    return present[0] if present else None
def canonical(x,y):
    while (i:=suffix(x,1)) is not None:
        a,u,v,_=codes[i]
        x,y=x//a,a*y+v
    word=[]
    while (i:=suffix(y,2)) is not None:
        word.append(i)
        y//=codes[i][0]
    return x,y,tuple(word)
rng=Random(45)
num=0
for k in range(29):
    for trial in range(10):
        x,z=[rng.getrandbits(4*k) for _ in range(2)]
        while suffix(x,1) is not None:x+=1
        while suffix(z,2) is not None:z+=1
        word=tuple(rng.randrange(2) for _ in range(k))
        y=z
        H=F(x,y)
        for i in reversed(word):
            A,U,V,r=codes[i]
            C=A-1-rev(V,r)
            H+=((Q*x+S)*(A-1)+45*C)*y+(Q*x+S)*V+45*p(V)
            y=A*y+V
            check(F(x,y)==H)
        left,right=x,y
        for j in range(k+1):
            check(canonical(left,right)==(x,z,word),x,z,word,j)
            check(F(left,right)==H-2*j,x,z,word,j)
            if j<k:
                A,U,V,_=codes[word[j]]
                left,right=A*left+U,right//A
            num+=1
print('mixed canonical orbit and consecutive output controls:',num)
for A,U,V,r in codes: print(r,'first repair from',F(0,V),'to',F(U,0),'radix',A)
