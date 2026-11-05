# seeds: all tournaments of order 5 with hc in {1,2}, as canonical digraph6 (gentourng + labelg)
import subprocess, itertools
n=5
out=subprocess.run(['/opt/homebrew/bin/gentourng','-q',str(n)],capture_output=True,text=True).stdout.split()
def d6(rows,n):
    bits=[(rows[i]>>j)&1 for i in range(n) for j in range(n)]
    s='&'+chr(n+63)
    for c in range(0,len(bits),6):
        chunk=bits[c:c+6]+[0]*(6-len(bits[c:c+6])); v=0
        for b in chunk: v=(v<<1)|b
        s+=chr(v+63)
    return s
def hc(rows,n):
    # cycles through 0: paths from 0
    N=1<<n; dp=[[0]*n for _ in range(N)]; dp[1][0]=1
    for S in range(1,N,2):
        for v in range(n):
            c=dp[S][v]
            if not c: continue
            nx=rows[v]&~S
            while nx:
                u=(nx&-nx).bit_length()-1; nx&=nx-1; dp[S|(1<<u)][u]+=c
    return sum(dp[N-1][v] for v in range(1,n) if rows[v]&1)
lines=[]
for t in out:
    rows=[0]*n; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if t[idx]=='1': rows[i]|=1<<j
            else: rows[j]|=1<<i
            idx+=1
    h=hc(rows,n)
    if h in (1,2): lines.append((h,d6(rows,n)))
can=subprocess.run(['/opt/homebrew/bin/labelg','-q'],input='\n'.join(s for _,s in lines)+'\n',capture_output=True,text=True).stdout.split()
seen={}
for (h,_),c in zip(lines,can): seen[c]=h
with open('cls_5.txt','w') as f:
    for c,h in sorted(seen.items(), key=lambda x:(x[1],x[0])): f.write('%d %s\n'%(h,c))
print('order 5 seeds:', sum(1 for h in seen.values() if h==1), 'one-cycle,', sum(1 for h in seen.values() if h==2), 'two-cycle')
