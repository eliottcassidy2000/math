"""Integer certificates for same-row folds of the physical two-tape map."""

Q, S = 5760, 3914
def require(ok):
    if not ok:
        raise AssertionError("row-fold identity failed")

def p(y):
    out = 0
    while y:
        y, b = divmod(y, 2)
        if not b:
            out += y + (y & 1)
    return out

def rev(b,k):
    return int(format(b, f"0{k}b")[::-1], 2)

def eq(x,y,z):
    return (Q*x+S)*(z-y) == 45*(p(y)-p(z))

def values(k,b,d):
    q=1<<k
    c=rev(b+d,k)-rev(b,k)
    return [p(q*e+b)-p(q*e+b+d)-c*e for e in (0,1)]

tops=[]
for b in range(19):
    d=rev(b+45,6)-rev(b,6)
    if d==45:
        tops.append(b)
        # u mod 2 controls the entire value mod 90.
        for e,k in enumerate(values(6,b,45)):
            require((45*e+k)%90 != S%90)
    else:
        require(d<=39)
require(tops==[0,2,16,18])
require(max(rev(b-19,6)-rev(b,6) for b in range(19,64))==43)
require(rev(67,7)-rev(22,7)==45 and values(7,22,45)[1]==-1)
require(rev(77,7)-rev(32,7)==87 and values(7,32,45)[0]==8)
for t in range(1000):
    x,u=t,128*t+87
    require(eq(x,128*u+22,128*u+67))
    x,u=15+29*t,1038+1920*t
    y=128*u+32
    require(eq(x,y,y+45))
    require(87*y == 128*(Q*x+S)+1760)
print("6-bit carry certificate, 7-bit fold identities, and 2000 exact row folds passed")
maxima=[]
for b in range(128-45):
    c=rev(b+45,7)-rev(b,7)
    V=values(7,b,45)
    if c == 90:
        require(all(k%90!=S%90 for k in V))
    elif c == 87:
        for e,k in enumerate(V):
            if (S-k-87*e)%3==0:
                maxima.append((b,e,k,87*b-128*k))
    else:
        require(c<=78 and max(V)<=92)
require(maxima==[(32,0,8,1760),(34,1,23,14),(48,0,20,1616),(50,1,35,-130)])
for t in range(1000):
    x,v=2+29*t,177+1920*t
    y=128*v+50
    require(eq(x,y,y+45) and 87*y == 128*(Q*x+S)-130)
print("Sharp lower-envelope arithmetic certificate and 1000 exact equality cases passed")
