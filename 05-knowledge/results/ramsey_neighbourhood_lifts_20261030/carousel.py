"""Constructive all-size carousel parity chart. See CIRCULAR.md."""
import time

def chart(n,epsilon,mask):
    odd=[i for i in range(n) if mask>>i&1];m=len(odd);out=[-1]*n
    if not m:return out
    if epsilon==0:
        word=[i%2 for i in range(m)]
        if m>1 and m%2:word[-1]=2
    elif m<=2:word=[2]*m
    else:
        r=n//2
        gaps=[(odd[(i+1)%m]-odd[i])%n for i in range(m)]
        forced=[i for i in range(m) if gaps[i-1]+gaps[i]>=r+2]
        if m==3:
            word=[2]*3
            free=next((i for i in range(m) if i not in forced),None)
            if free is not None:word[free]=0
        else:
            if len(forced)>2:raise RuntimeError('too many isolated marks')
            if len(forced)==2:
                start=forced[0] if (forced[0]+1)%m==forced[1] else forced[1]
                if (start+1)%m not in forced:raise RuntimeError('isolated marks not adjacent')
            else:start=forced[0] if forced else 0
            odd=odd[start:]+odd[:start]
            word=[2,2]
            if m%2:
                for j in range((m-1)//2):
                    word.append(j%2)
                    if len(word)<m:word.append(2)
            else:
                word.extend([0,1])
                for j in range(2,m//2):word.extend([2,j%2])
    if len(word)!=m:raise RuntimeError('chart length')
    for i,c in zip(odd,word):out[i]=c
    return out

def check(n,mask,epsilon,c):
    parts=[sum(1<<i for i,x in enumerate(c) if x==j) for j in range(3)]
    for start in range(n):
        row=sum(1<<((start+i)%n) for i in range(n//2+1));s=(row&mask).bit_count()
        for j in range(3):
            if (row&parts[j]).bit_count()>(s-epsilon+1)//2+epsilon*(j==2):
                raise RuntimeError((n,mask,epsilon,c,start,j))

if __name__=='__main__':
    for n in range(3,16,2):
        t=time.time()
        for mask in range(1<<n):
            for epsilon in (0,1):check(n,mask,epsilon,chart(n,epsilon,mask))
        print(n,2**(n+1),'charts passed',round(time.time()-t,3),flush=True)
