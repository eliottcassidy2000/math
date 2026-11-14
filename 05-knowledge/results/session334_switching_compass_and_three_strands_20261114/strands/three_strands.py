from collections import Counter
from itertools import product, permutations

def bit(d,h): return d>>(1-h)&1

def transfer(digits):
    """Return H_r for phase r=sum assigned run indices mod 3."""
    result=[0,0,0]
    for a,b in ((0,1),(1,0)):
        state=Counter({(0,0,0,0):1})
        for c in digits:
            nxt=Counter()
            for key,count in state.items():
                s, t=key[:3],key[3]
                for r in range(3):
                    if s[r]==0 and r and not bit(c,(a,b)[r-1]): continue
                    z=2 if r==2 or not bit(c,(a,b)[r]) else 1
                    nxt[s[:r]+(z,)+s[r+1:]+((t+r)%3,)]+=count
            state=nxt
        for (s,t,u,k),count in state.items():
            if s!=1 and (t==2 or (t==0 and (a,b)==(0,1))):result[k]+=count
    return tuple(result)

def direct(digits):
    n=len(digits)
    total=[0,0,0]
    for a,b in ((0,1),(1,0)):
        for rows in product(range(3),repeat=n):
            runs=[[v for v in range(n) if rows[v]==r] for r in range(3)]
            R,M,S=runs
            if ((R and bit(digits[R[-1]],a)) or
                (M and (not bit(digits[M[0]],a) or bit(digits[M[-1]],b))) or (
                not M and a==1) or (S and not bit(digits[S[0]],b))):
                continue
            total[sum(rows)%3]+=1
    return tuple(total)

def main():
    for n in range(6):
        for w in product(range(4),repeat=n):
            got=transfer(w); want=direct(w)
            if got!=want:raise ValueError((w,got,want))
    print('three-phase exact transfer equals literal enumeration on all 1365 words of length <=5')
    for k in range(1,8):
        print(k,[tuple(v-q[2] for v in q) for l in range(1,8) if (q:=transfer((1,1)+(0,)*k+(2,)*l+(3,1,1)))])
if __name__=='__main__':main()
