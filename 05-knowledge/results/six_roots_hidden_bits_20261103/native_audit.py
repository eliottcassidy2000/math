#!/usr/bin/env python3
"""Audit all one-birth outputs via permutation-derived response coefficients."""
import itertools,json
from pathlib import Path
P=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def coefficients(rows):
    n=len(rows);head=[0]*n;tail=[0]*n;middle=[[0]*n for _ in range(n)]
    for w in itertools.permutations(range(n)):
        bad=[i for i,(a,b) in enumerate(zip(w,w[1:])) if not rows[a]>>b&1]
        if len(bad)>1:continue
        cuts=range(n+1) if not bad else (bad[0]+1,)
        for t in cuts:
            if t==0:head[w[0]]+=1
            elif t==n:tail[w[-1]]+=1
            else:middle[w[t-1]][w[t]]+=1
    return head,tail,middle
def isstrong(r):
    for v in range(len(r)):
        seen={v};todo=[v]
        while todo:
            a=todo.pop()
            for b in range(len(r)):
                if r[a]>>b&1 and b not in seen:seen.add(b);todo.append(b)
        if len(seen)!=len(r):return False
    return True
def main():
    data=json.loads((P/'native_probe.json').read_text());base=data['source_rows'];totals={'195':0,'197':0,'199':0};tested=0
    for c in data['census']:
        t=c['address'];w=list(range(8)) if t is None else [i+1 for i in range(6) if t>>i&1]+[7,0]+[i+1 for i in range(6) if not t>>i&1]
        keep={(a,b) for a,b in zip(w,w[1:])};r=[0]*8
        for a in range(8):
            for b in range(8):
                if a==b:continue
                if (a,b) in keep or ((b,a) not in keep and base[b]>>a&1):r[a]|=1<<b
        head,tail,middle=coefficients(r);need(sum(head)==c['complement_H'],'intermediate H')
        local={h:0 for h in totals}
        for mask in range(256):
            bits=[mask>>i&1 for i in range(8)]
            h=sum(head[i]*(1-bits[i])+tail[i]*bits[i] for i in range(8))+sum(middle[i][j]*bits[i]*(1-bits[j]) for i in range(8) for j in range(8))
            born=[r[i]|(bits[i]<<8) for i in range(8)]+[sum((1-bits[i])<<i for i in range(8))]
            if str(h) in local and isstrong(born):local[str(h)]+=1
            tested+=1
        need(local==c['strong_target_hits'],('all cuts',t,local));
        for h in totals:totals[h]+=local[h]
    need(totals==data['total_strong_hits'],'total')
    # Rooted pair counts on each selected literal target by direct permutations/cuts.
    from audit import brute,count
    for h,w in data['best_port_witnesses'].items():
        got,k=brute(w['target_rows']);need(got==int(h) and [a-got for a in k]==w['J_all'],'best ports')
        need(w['J_all'][6]==0,'common physical zero port6')
    trib=[1,1,3]
    for _ in range(9):trib.append(sum(trib[-3:]))
    for n in range(1,13):
        q=[sum(1<<j for j in range(n) if j==i+1 or j<i-1) for i in range(n)]
        need(count(q)==trib[n-1],'inherited tribonacci recurrence control')
    for n in range(4,51):
        path=[n-2,n-1]+list(range(n-2));pairs={frozenset((a,b)) for a,b in zip(path,path[1:])}
        def before(i,j):return i<j if {i,j}!={0,n-1} else i==n-1
        for i,j in itertools.combinations(range(n),2):
            actual=before(i,j) if frozenset((i,j)) in pairs else not before(i,j)
            expected=(j==i+1) if j<n-2 else (i==n-2 and j==n-1)
            need(actual==expected,('native tribonacci map',n,i,j))
    out={'status':'PASS independent full one-birth coefficient audit and literal best-port permutation cuts','paths':65,'birth_outputs':tested,'strong_target_hits':totals,'common_zero_port':6,'native_tribonacci_pair_checks_through':50,'tribonacci_counts_through':12}
    (P/'native_audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: all16640 native/birth counts, strong filters and selected rooted profiles')
if __name__=='__main__':main()
