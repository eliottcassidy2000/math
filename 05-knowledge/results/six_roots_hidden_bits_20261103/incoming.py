#!/usr/bin/env python3
"""Independent phase census and literal source checks for the incoming wing chart."""
import collections,itertools,json
from pathlib import Path
from verify import profile
P=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def wing(m,l):
    # Integer formula first; phases0 are represented by positive exponent3.
    m=m or 3;l=l or 3
    return (2**m*(2**l-1)-2**(l-1)+2)%7,(2**(l-1)-1)%7
def joint(z):
    m,l,n,k=z;return wing(m,l),wing(n,k)
def obs(z):
    (x,b),(y,c)=joint(z);return (2*(x*y+b*c)-1)%7
def source(m,l,k,n):
    q=(2,4,9,19,39,15);sizes=(m,l,1,1,k,n)
    blocks=[]
    for i,s in enumerate(sizes):blocks += [i]*s
    rows=[0]*len(blocks)
    for a,b in itertools.combinations(range(len(blocks)),2):
        forward=blocks[a]==blocks[b] or q[blocks[a]]>>blocks[b]&1
        if forward:rows[a]|=1<<b
        else:rows[b]|=1<<a
    return rows
def main():
    phases=list(itertools.product(range(3),repeat=4))
    signatures=[tuple(obs(tuple((a+b)%3 for a,b in zip(z,t))) for t in phases) for z in phases]
    need(len(set(signatures))==81,'all future phases separate')
    joint_sizes=collections.Counter(joint(z) for z in phases)
    h_sizes=collections.Counter(obs(z) for z in phases)
    need(len(joint_sizes)==49 and max(joint_sizes.values())==9,'joint current fibres')
    need(dict(sorted(h_sizes.items()))=={0:18,1:12,2:11,3:8,4:13,5:13,6:6},'H fibres')
    controls=[]
    for m in (1,2,3):
        a=source(m,1,2,2);b=source(m,1,3,1)
        differences=sum(bool((a[i]^b[i])>>j&1) for i in range(len(a)) for j in range(i+1,len(a)))
        ha=profile(a)[0];hb=profile(b)[0]
        ca=profile(source(m,2,2,2))[0];cb=profile(source(m,2,3,1))[0]
        need(differences==1 and ha==hb==24*2**m+23,'hidden incidence')
        need((ca,cb)==(72*2**m+1,72*2**m+5),'literal revealing clone')
        controls.append({'m':m,'H':ha,'clone_H':[ca,cb],'changed_old_edges':differences})
    out={'incoming_revision':'2374ba81ef','status':'FINITE-EXACT complete phase universe and literal controls; no new universal return','phase_states':81,'future_tests_per_state':81,'current_joint_observations':49,'largest_joint_fibre':9,'joint_tag_bits':4,'Hmod7_fibre_sizes':dict(sorted(h_sizes.items())),'Hmod7_tag_bits':5,'one_bit_literal_controls':controls}
    (P/'incoming.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: 81 distinct phase futures; fibre sizes9/18 need4/5 bits; literal one-bit clone switch')
if __name__=='__main__':main()
