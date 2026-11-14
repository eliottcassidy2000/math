# Verify: F(u) = 908u + 9g(u) = 917u - 9Psi(u) - 9u0 + 9;  F(2^m V + W) = 917*2^m V - 9Psi(V) + 9 + Phi_V^m(W);
# Phi_V^m(b 2^(m-1) + W) = Phi_{2V+b}^{m-1}(W) + b (917*2^(m-1) - 9V - 9 v0);  and cross collisions vanish for large V.
import random
def g(U):
    if U==0: return 1
    k=U.bit_length(); s=1
    for t in range(k):
        b=(U>>t)&1
        if not b: s+=U>>(t+1)
        if t+1<k and ((U>>(t+1))&1) and not b: s+=1
    return s
def Psi(U): return 1+U-(U&1)-g(U)
def rev(W,m): return int(format(W,'0%db'%m)[::-1],2) if m>0 else 0
def Phi(V,W,m):
    if m==0: return 0
    return 917*W-9*Psi(W)-9*V*rev(W,m)-9*((V&1)*((W>>(m-1))&1))-9*(W&1)
random.seed(3); bad1=bad2=bad3=0
for _ in range(5000):
    u=random.randrange(1,1<<30)
    if 908*u+9*g(u)!=917*u-9*Psi(u)-9*(u&1)+9: bad1+=1
    m=random.randint(1,min(20,u.bit_length())); V=u>>m; W=u&((1<<m)-1)
    if 908*u+9*g(u)!=917*(1<<m)*V-9*Psi(V)+9+Phi(V,W,m): bad2+=1
    if m>=2:
        b=(W>>(m-1))&1; Wl=W&((1<<(m-1))-1)
        if Phi(V,W,m)!=Phi(2*V+b,Wl,m-1)+b*(917*(1<<(m-1))-9*V-9*(V&1)): bad3+=1
print("F identity mismatches",bad1,"| level split mismatches",bad2,"| top-bit recursion mismatches",bad3)
# cross collisions at level m: Phi_{2V}(W) = Phi_{2V+1}(W') + D_V; find the largest V (relative to 2^m) with any
for m in range(3,11):
    M=1<<(m-1); best=0
    for V in range(0,400*(1<<m),max(1,(1<<m)//8)):
        A={}
        for W in range(M): A.setdefault(Phi(2*V,W,m-1),[]).append(W)
        D=917*M-9*V-9*(V&1)
        hit=any((Phi(2*V+1,W2,m-1)+D) in A for W2 in range(M))
        if hit: best=max(best,V)
    print("m=%d: largest sampled V with a cross collision %d = %.1f * 2^m"%(m,best,best/(1<<m)))
