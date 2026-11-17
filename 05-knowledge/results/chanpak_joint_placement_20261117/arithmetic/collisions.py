"""Exact same-count seam fibers and distinct literal continuations."""
import json

def F(n):
    if n==0:return 1
    if n&1:return F(n//2)
    return F(n//2)+n//2+((n//2)&1)

def H(u,w):
    return 81+142*(u+w)+256*u*w+5*(F(u)+F(w))

def reverse(r,k):
    return int(f'{r:0{k}b}'[::-1],2)

def family(G,t=0,shift=0):
    if G<4 or G&(G-1) or t<0 or shift%2:raise ValueError()
    k=G.bit_length()-1;r=3;s=G-4
    cr=G-1-reverse(r,k);cs=G-1-reverse(s,k)
    A=16384*G*G-5
    B=5120*G*G-1024*G*(s-r)+10-20*(cr-cs)
    z0=next(z for z in range(1,16) if (shift+A*z+B)%15==0)
    y0=(shift+A*z0+B)//15
    z=z0+15*t;y=y0+A*t;x=y+z
    u=G*(16*x+5)+r;w=G*(16*y+4)+s
    return u,w,u+4*G,w-4*G

if __name__=='__main__':
    tests=0
    for k in range(2,15):
        G=1<<k
        for t in list(range(100))+[2**j+1 for j in (20,50,100,200)]:
            for shift in (0,-4,2,4):
                u,w,v,z=family(G,t,shift)
                assert u%2==v%2==1 and w%4==z%4==0
                assert H(v,z)-H(u,w)==shift
                assert H(v,z+2)-H(u,w+2)==shift+2053*G
                assert u.bit_length()==v.bit_length() and w.bit_length()==z.bit_length()
                tests+=1
    u,w,v,z=family(4)
    print(json.dumps({'checks':tests,'example':{'G':4,'addresses':[u,w,v,z], 'same_count':H(u,w),'first_continuation':H(u,w+2),'second_continuation':H(v,z+2),'continuation_difference':2053*4}},indent=2))
