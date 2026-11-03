"""Literal H497 reset, live-port sealing and balanced source reentry certificates."""
from math import gcd,isqrt
import json
from observer import HERE,old,junction,load,need,seal
sel=load('inherited_selector','05-knowledge/results/source_selection_closure_20261103/native_port_selector.py')

def chosen(k):
    return (1659*2**k+93,553*2**k-45) if k%2==0 else (1401*2**k+157,467*2**k+35)

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def replay_reset(w,u):
    receipt=old.reset(w,u);T=tuple(receipt['source']);n=len(T)
    V=junction.ear(T,((1<<n)-1)^receipt['birth_out'])
    for path in receipt['paths']:V=junction.native(V,path)
    V=tuple(row&((1<<n)-1) for row in V[:n])
    need(V==tuple(receipt['target'])==old.source(w,u,False),'independent native reset replay')
    need(junction.strong(T) and junction.strong(V),'strong reset endpoints')
    return V,receipt

def cycle(T,root,ports,b,t,k):
    """Generic supplied source; no assumption of ancestry from the three-cycle."""
    need(k in (3,5) and root in ports,'cycle source and mark')
    initial=T;ports=list(ports);events=[];stages=[];trace=[]
    def add(root,retain_new,kind,count):
        nonlocal T,ports
        need(root in ports,'birth root remains live')
        child=len(T);slot=ports.index(root);before=tuple(ports)
        T=junction.ear(T,1<<root);ports[1-slot]=child
        events.append({'kind':'birth','root':root,'child':child,'before_ports':before,'after_ports':tuple(ports),'retained_root':child if retain_new else root})
        stages.append({'kind':kind,'rows':T,'H':count})
        return child
    r=add(root,True,'N',2*b-t);leaves=[]
    for i in range(k):leaves.append(add(r,False,'root-block',b*(2**(i+1)+1)))
    z=leaves[-1];order=tuple(reversed(leaves))
    for i in range(k-1):
        expected,new_order,receipt=sel.adjacent_native_swap(T,order,i,f'cycle{k}_helper{i}')
        n=len(T);V=junction.ear(T,receipt['temporary_incoming_mask'])
        for path in (receipt['first_path'],receipt['second_path']):V=junction.native(V,path)
        V=tuple(row&((1<<n)-1) for row in V[:n])
        need(V==expected,'independent module-native replay')
        need(set(ports)=={r,z},'both actual selected ports retained during reset')
        T=V;order=new_order;receipt['ports']=tuple(ports);trace.append(receipt)
        events.append({'kind':'native-swap','receipt_index':len(trace)-1,'ports':tuple(ports)})
    need(order[-1]==z,'selected physical leaf is oldest')
    oldest=sel.block_profile(sel.N((b,t,t,b)),k-1,0)
    need(oldest[0]==b*(2**k+1),'oldest profile transfer')
    exit_h=2*b*(2**k+3**(k-1))-t
    y=add(z,True,'oldest-exit',exit_h)
    bp=(3*2**k+4*3**(k-1)+1)*b-2*t;tp=(2**k+1)*b
    add(y,False,'balanced-return',bp)
    need(sel.R(sel.N(oldest))==(bp,tp,tp,bp),'balanced return formula')
    # Full graph replay using only the recorded source, named births and native paths.
    U=initial;live=list(events[0]['before_ports'])
    for event in events:
        if event['kind']=='birth':
            need(tuple(live)==event['before_ports'] and event['root'] in live and len(U)==event['child'],'replayed birth clock and live root')
            U=junction.ear(U,1<<event['root']);live=list(event['after_ports'])
        else:
            receipt=trace[event['receipt_index']];n=len(U)
            need(tuple(live)==event['ports'],'native port continuity')
            U=junction.ear(U,receipt['temporary_incoming_mask'])
            U=junction.native(junction.native(U,receipt['first_path']),receipt['second_path'])
            U=tuple(row&((1<<n)-1) for row in U[:n])
    need(U==T and tuple(live)==tuple(ports) and junction.strong(T),'complete actual source replay')
    literal=[]
    for stage in stages:
        if len(stage['rows'])<=15:
            need(old.hp(stage['rows'])==stage['H'],'independent subset-DP stage count');literal.append(len(stage['rows']))
    return T,y,tuple(ports),bp,tp,{'k':k,'source_rows':initial,'initial_root':root,'source_pair':[b,t],
        'events':events,'native_trace':trace,'stages':stages,'final_rows':T,'final_root':y,'final_ports':ports,
        'return_pair':[bp,tp],'oldest_profile':oldest,'independently_counted_stage_orders':literal,
        'count_scope':'literal stage DP through order15; later counts by proved marked block transfer',
        'clocks':{'permanent_births':k+3,'native_moves':2*(k-1),'temporary_births':k-1,'temporary_deletions':k-1}}

def main():
    w,u=old.chart(3,4);T,reset=replay_reset(w,u);profiles=[junction.profile(T,r) for r in (7,8)]
    need(profiles==[(497,58,69,622),(497,56,61,528)],'H497 repaired rooted profiles')
    arithmetic=[]
    for k in range(1,181):
        A=seal(profiles[0],k);B=seal(profiles[1],k)
        need(A==(1659*2**k+93,553*2**k-45,553*2**k-45,1659*2**k+93),'closed first-port formula')
        need(B==(1401*2**k+157,467*2**k+35,467*2**k+35,1401*2**k+157),'closed second-port formula')
        need(gcd(A[0],A[1])==(19 if k%18==5 else 1),'exact first-port divisor clock')
        need(gcd(B[0],B[1])==(13 if k%12==8 else 1),'exact second-port divisor clock')
        b,t=chosen(k);bp,tp=421*b-2*t,33*b
        need(gcd(b,t)==gcd(bp,tp)==1 and bp%3 and tp%3==0 and tp%11==0 and bp%11,'direct primitive C5 controller entry')
        if k<=18:arithmetic.append({'k':k,'root':7 if k%2==0 else 8,'pair':[b,t],'C5_pair':[bp,tp],'C5_oldest_exit':226*b-t})
    for k in range(1,31):
        need(tuple(a-5*b+4*c for a,b,c in zip(chosen(k+4),chosen(k+2),chosen(k)))==(0,0),'order-four parity-merged recurrence')
    examples=[];seals={}
    for k in (1,2,3,6):
        root=7 if k%2==0 else 8
        U,r,ports,steps=junction.literal_seal(T,root,k,(7,8));b,t=chosen(k)
        if k<=3:need(junction.profile(U,r)==(b,t,t,b),'independent literal sealed profile')
        V=T;live=(7,8)
        for step in steps:
            need(step['root'] in live and step['newborn']==len(V),'seal live-root and birth clock')
            V=junction.ear(V,1<<step['root']);live=tuple(step['ports'])
        need(V==U and live==ports and junction.strong(U),'literal sealed source program')
        seals[k]={'k':k,'root_choice':root,'steps':steps,'rows':U,'root':r,'ports':ports,'pair':[b,t],
                  'profile_method':'independent subset DP' if k<=3 else 'proved R/N transfer'}
        result=cycle(U,r,ports,b,t,5);examples.append({'seal_k':k,'cycle':result[-1]})
    U=tuple(seals[1]['rows']);r=seals[1]['root'];ports=seals[1]['ports'];b,t=chosen(1)
    C3=cycle(U,r,ports,b,t,3)
    need(len(C3[0])==18 and old.hp(C3[0])==178561 and prime(178561),'literal 18-vertex prime return')
    C3[-1]['final_count_method']='independent full subset DP at order18'
    exit6=next(stage for stage in examples[-1]['cycle']['stages'] if stage['kind']=='oldest-exit')
    need(len(exit6['rows'])==24 and exit6['H']==23981447==7*3425921 and prime(3425921),'transfer-certified seven-prime exit')
    for k in range(2,32,2):
        b,t=chosen(k);need(226*b-t==7*(53483*2**k+3009),'all-even seven-multiple exit formula')
    out={'status':'PROVED_SCOPED_WITH_FINITE_EXACT_CONTROLS','reset':reset,'repaired_profiles':profiles,
         'arithmetic_period_controls':180,'first_selected_returns':arithmetic,'seals':list(seals.values()),'C5_entries':examples,
         'prime_C3_example':C3[-1],'seven_prime_exit':{'order':24,'H':23981447,'prime_quotient':3425921,
         'rows':exit6['rows'],'count_method':'proved actual marked block transfer; not order24 subset DP'},
         'total_C5_entry_clock':{'native_moves':10,'temporary_births':5,'temporary_deletions':5,'permanent_births':'k+10','final_order':'k+19'}}
    (HERE/'returns.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: named H497 reset, all-k parity entry, literal prime178561 and transferred7*3425921; full native traces')
if __name__=='__main__':main()
