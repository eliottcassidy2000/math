# independent check of Q2 bad classes at level r via explicit DFS on residues mod 3^(r+1)
import math,sys
L2,L3=math.log(2),math.log(3)
r=int(sys.argv[1]); M=3**(r+1)
def good(c):
    # DFS over reverse moves using exact rational path on representative c (valid for class while depth<=r)
    stack=[(c,0,0.0)]  # (value mod 3^(prec), depth, cost)
    while stack:
        x,d,cost=stack.pop()
        if d>=r: continue
        prec=r+1-d; mod=3**prec
        if x%3==0: continue
        k0=0 if x%3==1 else 1
        k=k0
        while True:
            step=k*L2-L3
            if cost+step-(r-d-1)*L3>=0 and cost+step>=0: 
                # even best future (all k=0 moves, -L3 each) cannot bring below 0
                break
            y=(pow(2,k,mod)*x-1)%mod
            y//=3
            if y%3!=0:
                if cost+step<0: return True
                stack.append((y%(3**(prec-1)),d+1,cost+step))
            k+=2
    return False
bad=[c for c in range(M) if c%3 and not good(c)]
print("r=",r,"bad count:",len(bad)); print(bad[:40])
