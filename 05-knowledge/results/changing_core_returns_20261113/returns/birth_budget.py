"""Exact H71 finite head for the all-order one-birth ray obstruction."""
from pathlib import Path
import hashlib,json
import inverse as inv
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def oddpart(n):
    while n%2==0:n//=2
    return n

def main():
    census=ROOT/'session212/cycle-spectrum/tournaments8.txt'
    presentations=[tuple(map(int,line.split())) for line in census.read_text().splitlines() if line.strip()]
    inv.need(len(presentations)==6880,'inherited complete order8 census size')
    cores=set();h71=0
    for T in presentations:
        H,_=inv.f.a.tables(T)
        if H[-1]==71:
            h71+=1
            for component in inv.f.a.ordered_components(T):
                if len(component)>1:cores.add(inv.f.a.induced(T,component))
        inv.f.a.tables.cache_clear()
    rows=[];checks=0;equality=[]
    for T in sorted(cores,key=lambda T:(len(T),T)):
        H,_=inv.f.a.tables(T);full=(1<<len(T))-1
        inv.need(inv.f.a.strong(T) and H[-1]==71,'extracted strong H71 component')
        fields=[];P2=None
        for root in range(len(T)):
            p,e,s=inv.f.a.bin_values(T,root,2);P2=p[2]
            q=H[full^(1<<root)];b=s[1];F=e[2]
            d=sum(bool(T[v]>>root&1) for v in range(len(T)))
            delta=F-71-q+b
            inv.need(F>=71 and q>=b and delta>=d*(d-1),'inherited cut-budget inequalities')
            inv.need(oddpart(F)!=405,'no core end-cover response with odd part405')
            if F==71:
                inv.need(d==1 and q==b,'sharp equality classification')
                equality.append(dict(rows=T,root=root,state=(71,b,q,F)))
            fields.append(dict(root=root,start=b,deletion=q,F=F,indegree=d,delta=delta));checks+=1
        inv.need(oddpart(P2)!=405,'no core P2 response with odd part405')
        rows.append(dict(rows=T,order=len(T),H=71,P2=P2,P2_oddpart=oddpart(P2),roots=fields))
        inv.f.a.tables.cache_clear()
    inv.need(len(cores)==11 and h71==18,'finite-head census controls')
    triangle=(2,4,1);p,e,s=inv.f.a.bin_values(triangle,0,2)
    inv.need((p[1],s[1],e[2])==(3,1,3),'strong equality positive control')
    witness=next((row,root['root']) for row in rows for root in row['roots'] if root['F']==404)
    core,mark=witness[0]['rows'],witness[1]
    parent=tuple(core)+((1<<len(core))-1,)
    p,e,s=inv.f.a.bin_values(parent,mark,2)
    inv.need((p[1],s[1],e[2])==(71,0,808),'literal neighboring-ray precursor')
    T,root=inv.f.a.birth(parent,mark);neighbor=[]
    for k in range(1,5):
        H=404*(1<<k)+71
        inv.need(inv.f.a.strong(T),'neighboring actual strong source')
        inv.need(inv.direct(T)==H and inv.direct(T,root)==71,'literal neighboring full H/start')
        deletion=inv.f.a.induced(T,[i for i in range(len(T)) if i!=root])
        inv.need(inv.direct(deletion)==71,'literal neighboring deletion')
        neighbor.append(dict(k=k,rows=T,root=root,state=(H,71,71,H)))
        T,_=inv.f.a.birth(T,root)
    rays={}
    for row in rows:
        for value,extra in [(row['P2'],0)]+[(r['F'],1) for r in row['roots']]:
            odd=oddpart(value);exponent=(value//odd).bit_length()-1+extra
            rays[odd]=min(rays.get(odd,exponent),exponent)
    inv.need(405 not in rays and rays[101]<=3,'exact missing and neighboring dyadic rays')
    data=dict(status='PASS',source=dict(path='session212/cycle-spectrum/tournaments8.txt',sha256=hashlib.sha256(census.read_bytes()).hexdigest(),order8_presentations=len(presentations),H71_presentations=h71),scope='Inherited complete census; independently counted selected profiles, not regenerated isomorphism classes',strong_core_presentations=rows,root_controls=checks,equality_controls=equality,P2_values=sorted({r['P2'] for r in rows}),maximum_F=max(x['F'] for r in rows for x in r['roots']),excluded_oddpart=405,balanced_t71_ray_support=[dict(oddpart=o,minimum_exponent=rays[o]) for o in sorted(rays)],neighboring_ray=dict(core_rows=core,core_root=mark,parent_rows=parent,literal_levels=neighbor))
    suffix='' if __debug__ else '.optimized'
    (HERE/('birth_budget_receipt'+suffix+'.json')).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS H71 finite head:',len(cores),'strong cores;',checks,'root profiles; no odd part405')

if __name__=='__main__':main()
