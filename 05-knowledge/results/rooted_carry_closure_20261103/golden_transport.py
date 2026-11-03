"""Keep the inherited Q4 golden phase separate from the new balanced-pair field clock."""
from functools import lru_cache
import hashlib,json
from observer import HERE,old,load,observer,root_profile,need,records
from returns import chosen
graded=load('inherited_graded','05-knowledge/results/source_selection_closure_20261103/graded_carrier.py')
golden=load('inherited_golden','05-knowledge/results/golden_inverse_defects_20261103/inverse.py')

@lru_cache(None)
def op(command):return graded.operator(2,command)

def obs(V,slot):
    h=e=s2=0
    for (sk,m),value in V.items():
        c=len(sk)+m
        if c==1:h+=value;e+=value*any(b==slot for a,b in sk)
        if c in (1,2):s2+=value*any(a==slot for a,b in sk)
    need((s2-e)%2==0,'integral golden extraction')
    return h,(s2-e)//2,(s2+e)//2

def phase(v):
    h,x,y=v;return [4*(h+2*x-y)%7,y%7,(h-y)%7]

def apply(T,ports,V,command):
    T,ports=graded.birth(T,ports,*command)
    return T,ports,graded.apply(V,2,op(command))

def sealed(T,V,k):
    ports=(7,8);slot=0 if k%2==0 else 1
    for letter in 'R'*k+'NR':
        command=(0,1<<slot,1-slot);T,ports,V=apply(T,ports,V,command)
        if letter=='N':slot=1-slot
    return T,ports,V,slot

def virtual_cycle(T,ports,V,slot,k):
    n=len(T);commands=[]
    def do(command):
        nonlocal T,ports,V
        T,ports,V=apply(T,ports,V,command);commands.append(command)
    do((0,1<<slot,1-slot));slot=1-slot # N, root r
    do((0,1<<slot,1-slot)) # designated oldest z first in the virtual history
    for _ in range(k-1):do((0,1<<slot,2)) # retain r and z
    slot=1-slot # select z
    do((0,1<<slot,1-slot));slot=1-slot # N at z
    do((0,1<<slot,1-slot)) # R at y
    # Virtual-name to actual-name map after the real native oldest-port selection.
    names=list(range(n))+[n,n+k]+list(range(n+1,n+k))+[n+k+1,n+k+2]
    U=[0]*len(T)
    for a,row in enumerate(T):
        for b in range(len(T)):
            if row>>b&1:U[names[a]]|=1<<names[b]
    return tuple(U),tuple(names[p] for p in ports),V,slot,commands,names

def main():
    w,u=old.chart(3,4);T=old.source(w,u,False);V=observer(T,(7,8));rows=[]
    receipts=json.loads((HERE/'returns.json').read_text())
    literal={x['seal_k']:x['cycle'] for x in receipts['C5_entries']}
    for k in range(1,19):
        S,ports,W,slot=sealed(T,V,k);b,t=chosen(k)
        need(root_profile(W,slot)==(b,t,t,b),'full observer retains exact sealed state')
        if k<=2:need(W==observer(S,ports),'independent literal complete sealed observer')
        U,pp,Z,ss,commands,names=virtual_cycle(S,ports,W,slot,5)
        need(root_profile(Z,ss)==(421*b-2*t,33*b,33*b,421*b-2*t),'full observer C5 return')
        if k in literal:
            need(U==tuple(literal[k]['final_rows']) and pp==tuple(literal[k]['final_ports']),'actual native endpoint equals named virtual history')
        g=obs(Z,ss);q4={name:phase(golden.step(g,golden.BANK[name])) for name in ('q0','q1','q2','q4','q5')}
        a,c,charge=phase(g)
        need(q4['q2']==[2*c%7,(a+4*c)%7,charge],'inherited golden field transport')
        rows.append({'k':k,'sealed_golden':obs(W,slot),'C5_golden':g,'C5_golden_phase':phase(g),
                     'all_five_Q4_phases':q4,'physical_ports':pp,'root_slot':ss,'virtual_commands':commands,
                     'virtual_to_actual_names':names,'full_observer_sha256':hashlib.sha256(json.dumps(records(Z),sort_keys=True).encode()).hexdigest()})
    # Supplied ready/converse theorem: exact formula, its phase clock, and no identification of charts.
    ready=[]
    for k in range(1,19):
        h=3*2**k-1;x=2**(k+1)-2;y=3**k+3*2**(k-1)-2
        direct=(h,x,x+1);opposite=(h,y,y+1)
        need((phase(direct)==phase(opposite))==(k%6==2),'incoming ready-family phase equality')
        if k<=6:ready.append({'k_mod6':k%6,'direct_phase':phase(direct),'converse_phase':phase(opposite)})
    (HERE/'golden_transport.json').write_text(json.dumps({'status':'PROVED_TRANSPORT_WITH_FINITE_EXACT_PHASE_EXAMPLES',
        'scope':'golden Q4 phase and balanced-pair beta phase are DIFFERENT sourced observables',
        'examples':rows,'incoming_ready_clock':ready},indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: full eighteen-polynomial return transport, actual native name matching, all five Q4 phases retained')
if __name__=='__main__':main()
