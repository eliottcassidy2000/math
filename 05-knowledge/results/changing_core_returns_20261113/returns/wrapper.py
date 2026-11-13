"""Independent ordered-bin controls for the arbitrary-source inverse interface."""
from pathlib import Path
import json,math
import inverse as inv
from growing import labelled
HERE=Path(__file__).resolve().parent

def carrier(B,m,l):
    order=list(range(1,l+1))+[0]+list(range(l+1,m));x=m;offset=m+1
    T=[0]*(m+1+len(B))
    for i,u in enumerate(order):
        T[u]|=1<<x
        for v in order[i+1:]:T[u]|=1<<v
    T[x]=sum(1<<(offset+i) for i in range(len(B)))
    for i,row in enumerate(B):T[offset+i]=(row<<offset)|((1<<m)-1)
    return tuple(T)

def profile(h,p,t,m,l):
    r=m-1-l;A=2**l*3**r;B=2**r;C=int(l==0)
    H=(2**m-1)*(p-h)+2*h
    P=2*((3**m-2**m)*t+(3*2**m-3**m-1)*p+(1-2**m)*h)
    S=(2*A-B)*p+(4*B-2*A-C)*h
    E=2**l*(p-h)-(p-2*h if r==0 else 0)
    F=2*(2**r*3**l-2**l)*(t-p)+3*2**l*p-2*2**l*h-(p-h if r==0 else 0)
    return dict(H=H,E=E,S2=S,P2=P,F=F,start=2**r*h)

def main():
    contexts=0
    for n in range(1,5):
        for B in labelled(n):
            bp,_,_=inv.f.a.bin_values(B,0,3)
            for m in range(2,5):
                for l in range(m):
                    T=carrier(B,m,l);p,e,s=inv.f.a.bin_values(T,0,2)
                    expected=profile(*bp[1:4],m,l)
                    actual=dict(H=p[1],E=e[1],S2=s[2],P2=p[2],F=e[2],start=s[1])
                    inv.need(expected==actual,('all-source inverse interface',B,m,l,actual,expected))
                    inv.need(inv.f.a.strong(T),'literal strong carrier')
                    contexts+=1;inv.f.a.tables.cache_clear()
    # Independent small literal control of the sibling H199 preparation.
    B=(6,12,8,1);last=carrier(B,3,2);prepared=carrier(B,3,1)
    inv.need(inv.direct(last)==199 and inv.direct(prepared)==199,'H199 neutral preparation')
    inv.need(inv.direct(last,0)==5 and inv.direct(prepared,0)==10,'physical prepared root')
    T=prepared;root=0
    for move in (False,True,False):
        T,new=inv.f.a.birth(T,root)
        if move:root=new
    inv.need(inv.direct(T)==4267 and inv.direct(T,root)==1041,'literal RNR full H and start')
    deleted=inv.f.a.induced(T,[i for i in range(len(T)) if i!=root])
    inv.need(inv.direct(deleted)==1041,'literal balanced deletion')
    miss_checks=0
    for m in range(2,21):
        for l in range(m):
            z=profile(5,32,117,m,l);G=z['S2']+z['P2']
            ns,_=inv.inverse(53084231,G,z['H'],z['E'],G-z['E']-z['F']-z['start'])
            inv.need(not ns,'declared strong4 wrapper target miss');miss_checks+=1
    rec=dict(status='PASS',contexts=contexts,strong4_k17_miss=dict(widths=[2,20],root_phases=miss_checks,arbitrary_outer_tail=True),H199=dict(last_rows=last,prepared_rows=prepared,returned_rows=T,root=root,state=(4267,1041,1041,4267)))
    suffix='' if __debug__ else '.optimized'
    (HERE/('wrapper_receipt'+suffix+'.json')).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS',contexts,'arbitrary-source interfaces; literal H199 return')

if __name__=='__main__':main()
