"""Independent (numpy) re-implementation of the closure used in THM-5055, plus explicit witness graphs.
Sets per order n (values are odd k; index (k-1)//2):
  E  even graph on n vertices                W  ... with no universal vertex
  P  ... with a pendant vertex               S  ... with a pendant whose neighbour is non-universal
  S2 ... with >= 2 pendants whose neighbours are non-universal
  S3 ... with a vertex of closed-neighbourhood size 3
  G  good = S2 and S3 simultaneously
Rules (lemmas in PROOF.md): pad; cone (W_{n-1} + 2^{n-1} -> E_n); E2-join (E_{n-2} + 3*2^{n-2} - 2 -> E_n, W->W);
  cone-with-leaf (E_{n-2} + 2^{n-1} -> P_n); R (S_{n-2} + 3*2^{n-3} -> S_n; S2_{n-2} + 3*2^{n-3} -> G_n);
  PH32 (S3_{n-3} + 7*2^{n-4} -> G_n); A (E_{n-4} + 2^{n-2} -> S2_n); B (E_{n-3} + 3*2^{n-2} - 2 -> S3_n).
Base: every even graph on <= NB vertices (unions of connected even census classes), with explicit edges."""
import sys, numpy as np, itertools, json
from collections import defaultdict
CENSUS, NB, NM = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
# ---------- explicit base graphs ----------
def parse(line):
    p=line.split(); n=int(p[0]); d=int(p[2]); par=p[3]; e=p[5] if len(p)>5 else ''
    adj=[0]*n
    for i in range(0,len(e),2):
        u=int(e[i],16); v=int(e[i+1],16); adj[u]|=1<<v; adj[v]|=1<<u
    return n,d,par,adj
def flags(adj):
    n=len(adj); deg=[bin(a).count('1') for a in adj]
    univ=[deg[v]==n-1 for v in range(n)]
    pend=[v for v in range(n) if deg[v]==1]
    pnu=[v for v in pend if not univ[(adj[v]&-adj[v]).bit_length()-1]]
    return dict(U=any(univ) if n>1 else True, P=len(pend)>0, S=len(pnu)>0, S2=len(pnu)>=2, S3=any(x==2 for x in deg))
def union(a,b):
    na=len(a); return a+[x<<na for x in b]
comps=defaultdict(list)   # n -> list of (d, adj) connected even
for line in open(CENSUS):
    n,d,par,adj=parse(line)
    if par=='E' and n<=NB: comps[n].append((d,adj))
# all unions with total order <= NB, components in nondecreasing (n, index) order to avoid repeats
base={}   # (n, setname, d) -> adj
def record(adj,d):
    n=len(adj); f=flags(adj)
    names=['E'] + (['W'] if not f['U'] else []) + (['P'] if f['P'] else []) + (['S'] if f['S'] else []) \
          + (['S2'] if f['S2'] else []) + (['S3'] if f['S3'] else []) + (['G'] if f['S2'] and f['S3'] else [])
    for s in names:
        base.setdefault((n,s,d),adj)
clist=[(n,i,d,adj) for n in sorted(comps) for i,(d,adj) in enumerate(comps[n])]
def rec(start,adj,d,tot):
    if adj: record(adj,d)
    for idx in range(start,len(clist)):
        n,i,dc,ac=clist[idx]
        if tot+n>NB: break
        rec(idx, union(adj,ac), d*dc, tot+n)
clist.sort(key=lambda t:(t[0],t[1]))
rec(0,[],1,0)
sys.stderr.write(f"base representatives: {len(base)}\n")
# ---------- numpy closure ----------
names=['E','W','P','S','S2','S3','G']
A={s:{} for s in names}
for n in range(1,NM+1):
    for s in names: A[s][n]=np.zeros(1<<(n-1),dtype=bool)
for (n,s,d) in base:
    if n<=NM: A[s][n][(d-1)//2]=True
def pad_into(s,n):
    A[s][n][:len(A[s][n-1])] |= A[s][n-1]
def shift_into(dst,n,src,m,off):   # dst_n |= src_m + off  (off even)
    h=off//2; L=len(A[src][m]); D=len(A[dst][n])
    if h>=D: return
    seg=A[src][m][:max(0,min(L,D-h))]
    A[dst][n][h:h+len(seg)] |= seg
for n in range(2,NM+1):
    for s in names: pad_into(s,n)
    A['W'][n][:len(A['E'][n-1])] |= A['E'][n-1]          # padding removes universal vertices
    A['S'][n][:len(A['P'][n-1])] |= A['P'][n-1]          # padding makes a pendant's neighbour non-universal
    if n>NB:
        shift_into('S',n,'S',n-2,3<<(n-3))                  # R
        A['P'][n] |= A['S'][n]; shift_into('P',n,'E',n-2,1<<(n-1))   # cone over (G u K1)
        shift_into('W',n,'W',n-2,(3<<(n-2))-2); A['W'][n] |= A['S'][n]
        A['E'][n] |= A['W'][n]; shift_into('E',n,'W',n-1,1<<(n-1)); shift_into('E',n,'E',n-2,(3<<(n-2))-2)
    # good / seeds (valid at every order)
    if n>=3: shift_into('G',n,'S2',n-2,3<<(n-3))
    if n>=4: shift_into('G',n,'S3',n-3,7<<(n-4))
    if n>=6: shift_into('S2',n,'E',n-4,1<<(n-2))
    if n>=5: shift_into('S3',n,'E',n-3,(3<<(n-2))-2)
    A['S2'][n] |= A['G'][n]; A['S3'][n] |= A['G'][n]
def has(s,n,k): return n>=1 and n<=NM and 1<=k<(1<<n) and k%2==1 and bool(A[s][n][(k-1)//2])
# ---------- checks ----------
def M(n): return 2**n-n if n%2 else 2**n-n-1
def formula(n,k):
    if k<1 or k%2==0 or k>M(n) or k in (3,7,13): return False
    if n%2==1 and n>=5 and k==2**n-n-2: return False
    return True
bad=0
for n in range(1,NM+1):
    for k in range(1,1<<n,2):
        if formula(n,k) and not has('E',n,k): bad+=1
        if has('E',n,k) and not formula(n,k): bad+=1; print('EXTRA',n,k)
print('E(n) = formula for n <=',NM,': mismatches',bad)
for c in (1,3,5,7,11,13):
    print('Gamma(+%d):'%c, [n for n in range(3,NM+1) if has('G',n,2**(n-1)+c)])
def L(k): return k.bit_length()
seedbad=0; worst=0; n0s={}
for j in range(1,1<<12,2):
    n0=None
    for n in range(3,NM):
        if 2**(n-1)>j and has('G',n,2**(n-1)-j) and has('G',n+1,2**n-j): n0=n; break
    if n0 is None or n0>L(j)+8: seedbad+=1
    else: n0s[j]=n0; worst=max(worst,n0-L(j))
print('seed base: odd j < 2^12 lacking two consecutive Gamma orders <= L(j)+8:',seedbad,' max n0-L(j)=',worst)
# ---------- explicit witnesses ----------
sys.setrecursionlimit(100000)
def pad(adj,k=1): return adj+[0]*k
def cone(adj):
    n=len(adj); return [a|(1<<n) for a in adj]+[(1<<n)-1]
def e2join(adj):
    n=len(adj); full=(1<<n)-1; return [a|(3<<n) for a in adj]+[full,full]
def ph(adj,u,l):
    n=len(adj); h=n; new=list(adj); hm=0
    for v in range(n):
        if v!=u: new[v]|=1<<h; hm|=1<<v
    for i in range(l): hm|=1<<(n+1+i)
    return new+[hm]+[1<<h]*l
def pick(adj,kind):
    n=len(adj); deg=[bin(a).count('1') for a in adj]
    univ=[deg[v]==n-1 for v in range(n)]
    if kind=='pnu':
        for v in range(n):
            if deg[v]==1 and not univ[(adj[v]&-adj[v]).bit_length()-1]: return v
    if kind=='r3':
        for v in range(n):
            if deg[v]==2: return v
    raise Exception('pick failed')
memo={}
def wit(s,n,k):
    key=(s,n,k)
    if key in memo: return memo[key]
    assert has(s,n,k), key
    r=None
    if (n,s,k) in base: r=base[(n,s,k)]
    elif s in ('E','W','P','S','S2','S3','G') and has(s,n-1,k): r=pad(wit(s,n-1,k))
    elif s=='S' and has('P',n-1,k): r=pad(wit('P',n-1,k))
    elif s=='W' and has('E',n-1,k): r=pad(wit('E',n-1,k))
    else:
        if s=='E':
            if has('W',n,k): r=wit('W',n,k)
            elif has('W',n-1,k-2**(n-1)): r=cone(wit('W',n-1,k-2**(n-1)))
            elif has('E',n-2,k-3*2**(n-2)+2): r=e2join(wit('E',n-2,k-3*2**(n-2)+2))
        elif s=='W':
            if has('W',n-2,k-3*2**(n-2)+2): r=e2join(wit('W',n-2,k-3*2**(n-2)+2))
            elif has('S',n,k): r=wit('S',n,k)
        elif s=='P':
            if has('S',n,k): r=wit('S',n,k)
            elif has('E',n-2,k-2**(n-1)): r=cone(pad(wit('E',n-2,k-2**(n-1))))
        elif s=='S':
            if has('S',n-2,k-3*2**(n-3)):
                g=wit('S',n-2,k-3*2**(n-3)); r=ph(g,pick(g,'pnu'),1)
        elif s in ('S2','S3') and has('G',n,k): r=wit('G',n,k)
        elif s=='S2':
            if has('E',n-4,k-2**(n-2)): r=pad(cone(pad(wit('E',n-4,k-2**(n-2)),2)))
        elif s=='S3':
            if has('E',n-3,k-3*2**(n-2)+2): r=e2join(pad(wit('E',n-3,k-3*2**(n-2)+2)))
        if r is None and s=='G':
            if has('S2',n-2,k-3*2**(n-3)):
                g=wit('S2',n-2,k-3*2**(n-3)); r=ph(g,pick(g,'pnu'),1)
            elif n>=4 and has('S3',n-3,k-7*2**(n-4)):
                g=wit('S3',n-3,k-7*2**(n-4)); r=ph(g,pick(g,'r3'),2)
    if r is None: raise Exception('no rule for %s'%(key,))
    assert len(r)==n, (key,len(r))
    memo[key]=r; return r
out=open(sys.argv[4],'w')
def emit(tag,s,n,k):
    g=wit(s,n,k)
    edges=' '.join(f'{u}-{v}' for u in range(n) for v in range(u+1,n) if (g[u]>>v)&1)
    out.write(f'{tag} {s} {n} {k} | {edges}\n')
WMAX=int(sys.argv[5])
cnt=0
for n in range(1,WMAX+1):
    for k in range(1,1<<n,2):
        if formula(n,k) and not has('E',n-1,k): emit('Ebase','E',n,k); cnt+=1
for c in (1,3,5,7,11,13):
    for n in (8,9): emit('off%+d'%c,'G',n,2**(n-1)+c); cnt+=1
for j,n0 in n0s.items():
    emit('seed-%d'%j,'G',n0,2**(n0-1)-j); emit('seed-%d'%j,'G',n0+1,2**n0-j); cnt+=2
out.close(); print('witnesses written:',cnt)
