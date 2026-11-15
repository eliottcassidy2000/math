from itertools import product

def paths(rows):
    n=len(rows)
    h=[1]+[0]*((1<<n)-1)
    d=[[0]*n for _ in h]
    for s in range(1,len(h)):
        for i in range(n):
            if s>>i&1:
                z=s^(1<<i)
                d[s][i]=sum(d[z][j] for j in range(n) if (rows[i]&z)>>j&1) if z else 1
        h[s]=sum(d[s])
    return h,d

def reverse(rows):
    return tuple(sum((row>>i&1)<<j for j,row in enumerate(rows)) for i in range(len(rows)))

def apex(rows,mask):
    n=len(rows)
    return tuple(row|(1<<n if not mask>>i&1 else 0) for i,row in enumerate(rows))+(mask,)

def weights(rows,mask):
    h,d=paths(rows)
    f=len(h)-1
    v=[h[f^s]*sum(d[s][i] for i in range(len(rows)) if mask>>i&1) for s in range(1,f+1)]
    return sum(v),sum(z if s.bit_count()%2 else -z for s,z in enumerate(v,1))

def sig(rows,mask):
    h,_=paths(rows)
    f=len(h)-1
    return h[-1],sum(h[s]*h[f^s] for s in range(f+1)),*weights(rows,mask),*weights(reverse(rows),f^mask)

def tt(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))

def join(X,Y):
    a,b=len(X),len(Y)
    return tuple(x|(((1<<b)-1)<<a) for x in X)+tuple(y<<a for y in Y)

def build(x, core, m, y):
    a=x.bit_length()
    raw=bin(x)[2:] if x else ''
    left_mark=sum((v=='1')<<i for i,v in enumerate(raw))
    b=y.bit_length()
    right_mark=sum((not (y>>i&1))<<i for i in range(b))
    whole=join(join(tt(a),core),tt(b))
    return apex(whole,left_mark|(m<<a)|(right_mark<<(a+len(core))))

def p(z):
    if not z:return 0
    a,b=divmod(z,2)
    return p(a)+(0 if b else a+a%2)

def strong(rows):
    for i in range(len(rows)):
        r=1<<i
        for k in rows:
            for j,x in enumerate(rows):
                if r>>j&1:r|=x
        if r!=(1<<len(rows))-1:return False
    return True

def main():
    r=(6,24,18,5,9)
    q=(0,2,1,4,3)
    if any(r[q[i]]!=((1<<q[(i+1)%5])|(1<<q[(i+2)%5])) for i in range(5)):
        raise AssertionError(q)
    one=[(r,5,41),(r,3,43)]
    r=(26,44,17,36,10,21)
    two=[((6,28,40,17,37,11),6,135),(r,10,137),(r,6,139)]
    for family, expected in [(one,(15,120,42,6,63,9)),(two,(45,432,129,15,258,30))]:
        for rows,m,v in family:
            if sig(rows,m)!=expected or not strong(apex(rows,m)) or paths(apex(rows,m))[0][-1]!=v:
                raise AssertionError((rows,m,sig(rows,m)))
        print('exact common signature and apex counts',expected,[n for _,_,n in family])
    tests=0
    nine=[]
    for family,(h,q,s,e,t,f) in [(one,(15,120,42,6,63,9)),(two,(45,432,129,15,258,30))]:
        for (rows,m,v),(x,y) in product(family,product(range(7),repeat=2)):
            got=paths(build(x,rows,m,y))[0][-1]
            want=v+h*(p(x)+p(y))+t*x+s*y+q*x*y+f*(x%2)+e*(y%2)
            if got!=want:raise AssertionError((rows,m,x,y,got,want))
            tests+=1
            if family==two and (x,y)==(1,1):nine.append(got)
    if nine != [999,1001,1003]:raise AssertionError(nine)
    print('literal two-tail path counts',tests)
    print('literal order-nine triple',nine)
if __name__=='__main__':main()
