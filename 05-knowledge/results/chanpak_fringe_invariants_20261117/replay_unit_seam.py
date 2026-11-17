from replay_joint_phase import direct

def word(F):return bin(F)[2:]
def revcomp(word):return ''.join(str(1-int(x)) for x in reversed(word))
def count(word):
    return 1+sum((1<<(j-i-1))+(j==i+1) for i,b in enumerate(word) if b=='1'
                 for j in range(i+1,len(word)) if word[j]=='0')
def literal(word):
    n=len(word)
    a=[((1<<n)-(2<<i))|((word[i]=='0')<<n) for i in range(n)]
    a.append(sum((b=='1')<<i for i,b in enumerate(word)))
    f=(1<<(n+1))-1
    def paths(end,m):
        if not m:return 1
        return sum(paths(i,m^(1<<i)) for i in range(n+1) if m>>i&1 and a[i]>>end&1)
    return sum(paths(i,f^(1<<i)) for i in range(n+1))
last=0
for F in range(1,1001):
    for d in (1,2,4,15):
        P,S=word(F),revcomp(word(F+2*d))
        w0,w1=P+'0'+S,P+'1'+S
        lo,hi=count(w0),count(w1)
        predicted=2*F*(F+2*d)+F+count(P)+count(S)-1+(F%2)
        if lo!=predicted or hi-lo!=2*d:raise ValueError((F,d,lo,hi))
        if d==1:
            if lo<=last or hi>2*F*F+7*F+6:raise ValueError(F)
        if F<5 and d in (1,2) and (literal(w0),literal(w1))!=(lo,hi):
            raise ValueError('literal control')
    last=count(word(F)+'1'+revcomp(word(F+2)))
print("4000 exact windows; 8 independent small literal controls")

def inv_dist(w):
    out=[0]*len(w)
    for i,b in enumerate(w):
        if b=='1':
            for j in range(i+1,len(w)):
                if w[j]=='0':out[j-i]+=1
    return out
for F in range(1,256):
    for j in range(1,6):
        if F>>j&1:continue
        P,S=word(F),revcomp(word(F+(1<<j)))
        lo,hi=P+'0'+S,P+'1'+S
        delta=[y-x for x,y in zip(inv_dist(lo),inv_dist(hi))]
        if sum(delta)!=1 or delta[j+1]!=1 or sum(map(abs,delta))!=1:
            raise ValueError((F,j,delta))
        if j==1:
            a= len(lo)+2
            h0,p0=direct(lo,1<<a); h1,p1=direct(hi,1<<a)
            if ((h1-h0)%(1<<a),(p1-p0)%(1<<a))!=(2,1<<len(lo)):
                raise ValueError((F,j,h0,p0,h1,p1))
print("carry-free switches preserve all inversion distances except the selected one")

def slots(word,q):
    L=len(word)
    return q**(L+1)+sum(k*q**(L-d)*((q+1)**(d-1)+(q-1)**(d-1))
                         for d,k in enumerate(inv_dist(word)) if d)
for u in range(300):
    P,S=word(20+64*u),revcomp(word(32+64*u))
    lo,hi=P+'0'+S,P+'1'+S
    if count(hi)-count(lo)!=12 or slots(hi,2)!=slots(lo,2) or slots(hi,3)-slots(lo,3)!=-100*3**(len(lo)-5):
        raise ValueError(('exact full coordinate',u))
print("300 growing switches: ΔH=12, exact ΔP2=0, nonzero ΔP3")

lams=[(3**s+1)<<(6-s) for s in range(7)]
for F in range(128):
    for dif in (5,6,7):
        Z=F+dif
        if Z>=128:continue
        if dif+Z%2-F%2==6 and sum(((Z>>s&1)-(F>>s&1))*x for s,x in enumerate(lams))==0:
            raise ValueError(('step six',F,Z))
print("sharpness: no +6 remains after carry-depth bound")
