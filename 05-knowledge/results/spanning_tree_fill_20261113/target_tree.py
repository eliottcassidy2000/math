"""Exact odd-target macro tree; graph witnesses only for the five bypasses."""
from pathlib import Path
from itertools import permutations
from subprocess import check_output
import json

HERE=Path(__file__).resolve().parent
FORBIDDEN={7,21}
SOURCE_REF='dfaf7c7e88dc'


def need(condition,message):
    if not condition:raise RuntimeError(message)


def parent(h):
    need(h>1 and h%2==1,'odd nonroot')
    if h%4==3:return (h-1)//2
    if h%12==5:return (h-2)//3
    if h%12==1:return (h-4)//3
    return (2*h-3)//3


def repaired_parent(h):
    p=parent(h)
    while p in FORBIDDEN:p=parent(p)
    return p


def near_transitive(n):
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    rows[0]^=1<<(n-1);rows[-1]|=1
    return rows


def ordinal(a,b):
    n=len(a);tail=((1<<len(b))-1)<<n
    return [row|tail for row in a]+[row<<n for row in b]


def hp(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            rem=rows[v]&~mask
            while value and rem:
                bit=rem&-rem;rem-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])


def literal(rows):
    return sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(len(rows))))


def main():
    cap=100001;depth={1:0};bypasses=[]
    for h in range(3,cap+1,2):
        p=parent(h)
        need(0<p<h and p%2==1,'strict odd descent')
        need(3*(p+1)<=2*(h+1),'contraction')
        if h in FORBIDDEN:continue
        q=repaired_parent(h);need(q not in FORBIDDEN,'avoids holes')
        depth[h]=depth[q]+1
        if p!=q:bypasses.append([q,h])
    need(sorted(bypasses)==[[3,15],[3,25],[13,33],[13,43],[13,65]],'complete bypass head')
    witnesses={15:ordinal(near_transitive(3),near_transitive(4)),
               25:ordinal(near_transitive(4),near_transitive(4)),
               33:near_transitive(7),65:near_transitive(8)}
    source='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n7_representatives.txt'
    bank=check_output(['git','show',SOURCE_REF+':'+source],cwd=HERE,text=True)
    for line in bank.splitlines():
        rows=list(map(int,line.split()))
        if hp(rows)==43:witnesses[43]=rows;break
    need(43 in witnesses,'inherited43 witness')
    controls=[]
    for target,rows in sorted(witnesses.items()):
        need(hp(rows)==literal(rows)==target,'independent path counts')
        controls.append(dict(target=target,order=len(rows),rows=rows,subset_count=target,permutation_count=target))
    need(159252697%4==1 and 3*159252697+2==477758093,'blocked next interface')
    result=dict(status='PASS',head=cap,allowed_targets=len(depth),maximum_depth=max(depth.values()),
                bypasses=bypasses,base_witnesses=controls,inherited43_bank=source,inherited_git_ref=SOURCE_REF,
                scope='Arithmetic tree and five actual bases only; universal tournament lift OPEN')
    (HERE/'target_tree_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(depth),'arithmetic targets;',len(controls),'literal base witnesses')


if __name__=='__main__':main()
