#!/usr/bin/env python3
"""Build the literal carrier tournaments for solver hits (hostsolve output) and count their Hamiltonian
paths by subset DP (hp2), independently of the closed form.  usage: bf_check.py small_hosts.txt"""
import re, subprocess, sys
import re, subprocess, sys
# build carrier tournaments T_C(X,Y) for the 3-cycle core with one marked vertex; compute H and P(2)

def p(u):
    if u<=1: return 0
    return p(u>>1) if u&1 else p(u>>1)+(u>>1)+((u>>1)&1)
def carrier(core_adj, core_mark, X, Y):
    # spine order: g1 (marked), X-tape (digits of X below leading 1, MSB first, mark=digit), core, Y-tape (j-th unmarked iff bit j-1 of Y is 1), g2 (unmarked)
    spine=[]; marks=[]
    spine.append(('g1',)); marks.append(1)
    xb=bin(X)[3:]
    for ch in xb: spine.append(('x',)); marks.append(int(ch))
    k=len(core_adj); base=len(spine)
    for i in range(k): spine.append(('c',i)); marks.append(core_mark[i])
    ylen=Y.bit_length()-1
    for j in range(1,ylen+1):
        bit=(Y>>(j-1))&1
        spine.append(('y',)); marks.append(0 if bit==1 else 1)
    spine.append(('g2',)); marks.append(0)
    m=len(spine); n=m+1; apex=m
    out=[0]*n
    for i in range(m):
        for j in range(i+1,m):
            if spine[i][0]=='c' and spine[j][0]=='c':
                a,b=spine[i][1],spine[j][1]
                if core_adj[a][b]: out[i]|=1<<j
                else: out[j]|=1<<i
            else: out[i]|=1<<j        # ordinal sum: earlier beats later
    for i in range(m):
        if marks[i]: out[apex]|=1<<i   # apex beats marked
        else: out[i]|=1<<apex
    return n,out

lines=open(sys.argv[1] if len(sys.argv)>1 else 'small_hosts.txt').read().strip().split('\n'); inp=[]; meta=[]
for ln in lines:
    m=re.match(r'N=(\d+) core\(order (\d+); n,V,out\)=\[(.*?)\] X=(\d+) Y=(\d+)',ln)
    N,o,core,X,Y=int(m.group(1)),int(m.group(2)),list(map(int,m.group(3).split())),int(m.group(4)),int(m.group(5))
    n,V,outs=core[0],core[1],core[2:]
    adj=[[1 if (outs[i]>>j)&1 else 0 for j in range(n)] for i in range(n)]
    mark=[(V>>i)&1 for i in range(n)]
    nn,out=carrier(adj,mark,X,Y); inp.append('%d %s'%(nn,' '.join(map(str,out)))); meta.append((N,nn))
res=subprocess.run([sys.argv[2] if len(sys.argv)>2 else './hp2'],input='\n'.join(inp)+'\n',capture_output=True,text=True).stdout.split('\n')
for (N,nn),r in zip(meta,res): print('target',N,'vertices',nn,'brute-force H =',r.split()[0],'OK' if int(r.split()[0])==N else 'MISMATCH')
