from math import isqrt,gcd

def insist(ok, msg):
    if not ok: raise RuntimeError(msg)
# Algebraically independent direct scan over all odd legs, forward and reverse.
forward=[]; reverse=[]
for a in range(1,1000000,2):
    z=3*a+1
    while z%2==0:z//=2
    if z>a:
        b2=z*z-a*a;b=isqrt(b2)
        if b*b==b2:forward.append((a,b,z))
    # Need not search all c: q possible while 2^q*a<=3c+1; use bound c derived from q.
    for q in range(2,20):
        d=(1<<q)*a-1
        if d%3:continue
        c=d//3
        if c<=a:continue
        b2=c*c-a*a;b=isqrt(b2)
        if b*b==b2:reverse.append((a,b,c,q))
insist(forward==[(3,4,5),(1155,1292,1733),(372099,416020,558149)], str(forward))
# Independent norm recurrence iterates are included in direct reverse census.
x,y=381,144
for i in range(3):
    a=(x+4)//7;b=y//3;c=(4*a-1)//3
    if a<1000000:insist((a,b,c,2) in reverse,'missing reverse recurrence')
    x,y=127*x+336*y,48*x+127*y
print('forward, odd a<10^6:',forward)
print('reverse, q in [2,19], odd a<10^6:',len(reverse),'including',reverse[:10])
