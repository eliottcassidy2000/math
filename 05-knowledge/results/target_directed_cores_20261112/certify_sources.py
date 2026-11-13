"""Retain literal source addresses and certify the selective bridge at the hole.

Large counts use proved full-cover compiler identities. Only sources of order
at most13 are additionally counted by literal subset DP. No large census is
claimed by the witness certificates.
"""
from pathlib import Path
import json
from build_core_bank import SINGLE,advance
HERE=Path(__file__).resolve().parent
def need(p,msg):
    if not p:raise RuntimeError(msg)
def word(x):return tuple((x>>i)&1 for i in range(x.bit_length()))
def funnel(rows,root,w):
    n=len(rows);t=len(w);x=n+t;out=list(rows)+[0]*(t+1)
    for i in range(n):
        for j in range(n,n+t):out[i]|=1<<j
        if i==root:out[x]|=1<<i
        else:out[i]|=1<<x
    for i,b in enumerate(w):
        v=n+i
        for j in range(i+1,t):out[v]|=1<<(n+j)
        if b:out[v]|=1<<x
        else:out[x]|=1<<v
    return tuple(out),x
def birth(rows,root):
    n=len(rows);out=list(rows)+[0];out[root]|=1<<n
    out[n]=((1<<n)-1)^(1<<root)
    return tuple(out),n
def constructor(address):
    rows=(0,);root=0;state=SINGLE
    for text in address:
        w=tuple(map(int,text));need(bool(w) or len(rows)>1,'strong empty-word guard')
        rows,root=funnel(rows,root,w);state=advance(state,w)
    parent=rows;pred=root;rows,root=birth(rows,root)
    h=state[0]+state[4]-1;q=state[0]
    return dict(rows=rows,root=root,predecessor=pred,order=len(rows),parent_rows=parent,
                parent_state=state,state=(h,q-1,q,h+1),address=address)
def hp(rows,start=None):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        if start is None or start==v:dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,val in enumerate(dp[mask]):
            if not val:continue
            rem=rows[v]&~mask
            while rem:
                bit=rem&-rem;rem-=bit;dp[mask|bit][bit.bit_length()-1]+=val
    return sum(dp[-1])
def validate(rows):
    n=len(rows)
    for i in range(n):
        need(not(rows[i]>>i&1),'no loops')
        for j in range(i):need(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1,'tournament pair')
def main():
    bank=json.loads((HERE/'core_bank.json').read_text())['cores']
    lookup={(r['order'],r['rank']):r for r in bank}
    search=json.loads((HERE/'target_inverse_receipt.json').read_text())['targets']
    witnesses=[];literal=0
    for r in search:
        if r.get('found') is False:continue
        base=lookup[r['core_order'],r['core_rank']]
        address=base['address']+[''.join(map(str,word(r['word_integer'])))]
        cert=constructor(address);validate(cert['rows'])
        need(cert['state'][0]==r['target'],'full compiler target')
        cert.update(k=r['k'],target=r['target'],core_order=base['order'],core_rank=base['rank'])
        cert['literal_checked']=cert['order']<=13
        if cert['literal_checked']:
            need(hp(cert['rows'])==r['target'],'literal total')
            need(hp(cert['parent_rows'],cert['predecessor'])==1,'literal unique start')
            need(hp(cert['parent_rows'])==cert['state'][2],'literal deletion')
            literal+=1
        witnesses.append(cert)
    (HERE/'source_certificates.json').write_text(json.dumps(dict(status='PASS',scope='Declared positive targets only; full compiler proofs and small literal controls',literal_controls=literal,witnesses=witnesses),indent=2)+'\n')
    # The raw 405 ray at k17, whose count is inherited from MODULE_RETURN.
    k=17;D=1<<k
    parent,p=funnel((2,4,1),0,(1,)*k+(0,0,0,1))
    rows,r=birth(parent,p)
    z=(135*D+23,27*D+26,27*D+27,135*D+24)
    for _ in range(2):
        h,s,q,F=z;z=(h+F-s,h-s,h,h+F);rows,r=birth(rows,r)
    need(z[0]==405*D+71,'raw target')
    old=rows;oldroot=r;oldpred=r-1;module=list(range(3,3+k))
    target=constructor([''.join(map(str,word(45157534)))])
    need(target['state'][0]==3*z[0]+2,'next ternary target')
    new=target['rows'];newroot=target['root'];newpred=target['predecessor']
    mapping={v:j+1 for j,v in enumerate(module)}
    mapping.update({oldroot:newroot,oldpred:newpred})
    available=[v for v in range(len(new)) if v not in mapping.values()]
    rest=[v for v in range(len(old)) if v not in mapping]
    mapping.update(zip(rest,available))
    need(len(mapping)==len(old),'all old vertices retained')
    for a in range(len(old)):
        if a!=oldroot:need(((old[oldroot]>>a)&1)==((new[newroot]>>mapping[a])&1),'entire root star retained')
    for a in module:
        for b in module:
            need(((old[a]>>b)&1)==((new[mapping[a]]>>mapping[b])&1),'module internal arcs retained')
    flips=[(a,b) for a in range(len(old)) for b in range(a+1,len(old)) if ((old[a]>>b)&1)!=((new[mapping[a]]>>mapping[b])&1)]
    bridge=dict(status='PASS',scope='One explicit direct target-changing macro; individual flips are not H-neutral',source_rows=old,source_state=z,source_root=oldroot,source_predecessor=oldpred,module=module,target=target,retained_map=sorted(mapping.items()),new_vertices=[v for v in range(len(new)) if v not in mapping.values()],changed_old_pairs=flips,selective_predecessor_pairs=[pair for pair in flips if oldpred in pair])
    (HERE/'gap_bridge.json').write_text(json.dumps(bridge,indent=2)+'\n')
    print('PASS',len(witnesses),'source certificates;',literal,'literal controls; bridge retains',len(old),'vertices and changes',len(flips),'old pairs')
if __name__=='__main__':main()
