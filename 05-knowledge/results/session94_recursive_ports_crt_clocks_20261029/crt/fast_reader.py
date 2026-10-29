"""Exact inverse odd-period golden reader: 9 carry candidates, no word census."""
def fib(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a

def mul(x,y):
    a,b=x;c,d=y
    return a*c+b*d,a*d+b*c+b*d

def sign(x):
    u,v=x;A=2*u+v;B=v
    if B==0:return (A>0)-(A<0)
    if A>=0 and B>0:return 1
    if A<=0 and B<0:return -1
    z=A*A-5*B*B
    return ((z>0)-(z<0))*(1 if A>0 else -1)

def inverse_reader(n,target):
    if n<1 or n%2!=1:raise ValueError('positive odd period required')
    a,b=fib(n-1)-1,fib(n);L=abs(a*a+a*b-b*b)
    u,v=target;A=u*(a+b)-v*b;B=-u*b+v*a
    ca,cb=A//L,B//L
    powers=[(1,0)]
    for i in range(1,n):powers.append(mul(powers[-1],(0,1)))
    solutions=set()
    for bb in range(cb,cb+3):
        for aa in range(ca-1,ca+2):
            x,y=mul((a,b),(aa,bb));x+=u;y+=v
            if sign((x,y))<0 or sign((x-a,y-b))>0:continue
            w=0
            for i in range(n-1,-1,-1):
                p,q=powers[i]
                if sign((x-p,y-q))>=0:
                    w|=1<<i;x-=p;y-=q
            if (x,y)!=(0,0):continue
            if w&(w>>1) or ((w&1) and (w>>(n-1)&1)):continue
            solutions.add(w)
    if len(solutions)!=1:raise RuntimeError(('decoder not unique',n,target,solutions))
    return solutions.pop()

def reflection_residue(n,x,N,r):
    w=inverse_reader(n,(x,0))
    return sum(pow(r,(-i)%n,N) for i in range(n) if w>>i&1)%N

if __name__=='__main__':
    import json
    from check_crt import words,setup
    checks=0
    for n in (1,3,5,9,15):
        for w in words(n):
            u=v=0;p,q=1,0
            for i in range(n):
                if w>>i&1:u+=p;v+=q
                p,q=q,p+q
            if inverse_reader(n,(u,v))!=w:raise RuntimeError('source decode mismatch')
            checks+=1
    L,g,M,r,keys,R,T,U=setup(25)
    ws=list(words(25))
    for i in range(0,L,997):
        if inverse_reader(25,(keys[i][1],0))!=ws[i]:raise RuntimeError('integer decode mismatch')
        checks+=1
    n=55;a,b=fib(n-1)-1,fib(n);N=abs(a*a+a*b-b*b);r=(-a*pow(b,-1,N))%N
    factors=[(11,2),(199,1),(331,1),(39161,1)];divs=[1]
    for p,e in factors:divs=[d*p**j for d in divs for j in range(e+1)]
    certificate=[]
    for d in sorted(set(divs)-{1,N}):
        for x in range(100):
            y=x+d;rx=reflection_residue(n,x,N,r);ry=reflection_residue(n,y,N,r)
            if (rx-ry)%d:
                certificate.append({'d':d,'x':x,'y':y,'word_x':inverse_reader(n,(x,0)),'word_y':inverse_reader(n,(y,0)),'Rx':rx,'Ry':ry});break
        else:raise RuntimeError(('no witness found',d))
    print(json.dumps({'decoder_source_checks':checks,'period':n,'N':N,'root':r,'factorization':factors,'certificate':certificate},indent=2))
