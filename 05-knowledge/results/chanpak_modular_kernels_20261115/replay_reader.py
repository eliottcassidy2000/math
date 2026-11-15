"""Fixed-core mask readers on the two-tape fold; pure Python."""
OUT=(4,1,130,7,143,159,63,75)

def count(mark, extras):
    # Guards 1,0, protected eight-vertex core, optional terminal 0 singleton.
    n=11+abs(extras)
    ap=n-1
    rows=[sum(1<<j for j in range(i+1,ap)) for i in range(ap)]
    start=1+int(extras<0)
    mask_core=255<<start
    rows[ap-1] |= 1<<ap
    rows[0] &= ~(1<<ap)
    marks=(1<<start)-1+(mark<<start)
    for i,r in enumerate(OUT,start):
        rows[i]=(rows[i]&~mask_core)|(r<<start)
        if not((marks>>i)&1): rows[i]|=1<<ap
    if extras>0: rows[ap-2]|=1<<ap
    rows.append(marks)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for m in range(1,2**n-1):
        left=(2**n-1)^m
        for i,val in enumerate(dp[m]):
            w=rows[i]&left
            while w and val:
                b=w&-w; w^=b
                dp[m|b][b.bit_length()-1]+=val
    return sum(dp[-1])

def main():
    B0=count(35,0); S0=count(35,1)-B0
    if (B0,S0)!=(2855,3914):raise AssertionError((B0,S0))
    data=[]
    for mark in range(256):
        B=count(mark,0)
        S=count(mark,1)-B
        T=count(mark,-1)-B
        data.append((B,S,mark,T))
    near=sorted((abs(S-S0),B,S,mark) for B,S,mark,_ in data if S!=S0)
    print("smallest nonzero carriage changes:",near[:10])
    from math import gcd
    g=0
    for _,S,_,_ in data:g=gcd(g,abs(S-S0))
    if g!=2 or any(not(2131<=B<=4671 and 2880<=T<=5670) for B,S,mask,T in data):
        raise AssertionError("mask spectrum bounds")
    assert {m:(B,S,T) for B,S,m,T in data if m in (35,38,152)} == {35:(2855,3914,4034),38:(3141,3934,4414),152:(3485,3896,4976)}
    print("gcd of all carriage changes:",g)
    print("ranges B,T:", [(min(d[i] for d in data),max(d[i] for d in data)) for i in [0,3]])
    for B,S,mask,T in data:
        if mask in (35,152,38): print(mask, B,S,T)

if __name__=="__main__":main()
