import subprocess, itertools
def parse(s,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if s[k]=='1': out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def H_P2(adj,n):
    full=(1<<n)-1; cnt=[0]*(1<<n); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            c=dp[S][v]
            if not c: continue
            cnt[S]+=c; nb=adj[v]&~S
            while nb:
                u=(nb&-nb).bit_length()-1; nb&=nb-1; dp[S|1<<u][u]+=c
    return cnt[full], sum(cnt[S]*cnt[full^S] for S in range(1,full))
minratio=None; tot=0
for n in range(1,9):
    lines=[l for l in subprocess.run(['gentourng',str(n)],capture_output=True,text=True).stdout.split() if not l.startswith('>')] if n>1 else ['']
    for s in lines:
        adj=parse(s,n) if n>1 else [0]
        h,P=H_P2(adj,n); Pi=P+2*h; tot+=1
        r=Pi/h
        if minratio is None or r<minratio[0]: minratio=(r,n,h,Pi)
print('classes checked',tot,'; minimum of Pi_2/H =',minratio[0],'at order',minratio[1],'(H, Pi_2) =',minratio[2:])
print('so 3 Pi_2^2 - H^2 >= (3*%.3f^2 - 1) H^2 on all of them; isotropic ratio 1/sqrt3 = 0.577 is excluded'%minratio[0])
