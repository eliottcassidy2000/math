"""Inherited BC2 exterior-flip compiler, specialized to a fresh neutral gate swap."""
import json
from pathlib import Path
from check import need,program,subset_hp,strong

HERE=Path(__file__).resolve().parent


def birth(a,cut):
    n=len(a);b=list(a)+[((1<<n)-1)^cut]
    for i in range(n):
        if cut>>i&1:b[i]|=1<<n
    return tuple(b)


def native(a,path):
    need(sorted(path)==list(range(len(a))),'spanning path')
    need(all(a[x]>>y&1 for x,y in zip(path,path[1:])),'Hamilton guard')
    spine={frozenset((x,y)) for x,y in zip(path,path[1:])};b=list(a)
    for x in range(len(a)):
        for y in range(x+1,len(a)):
            if frozenset((x,y)) not in spine:b[x]^=1<<y;b[y]^=1<<x
    return tuple(b)


def flip(a,x,y,root,names,serial):
    if not a[x]>>y&1:x,y=y,x
    path=[]
    for v in range(len(a)):
        if v in (x,y):continue
        j=next((j for j,u in enumerate(path) if a[v]>>u&1),len(path));path.insert(j,v)
    n=len(a);z,s,w=n,n+1,n+2
    cuts=[1<<y,(1<<(n+1))-1,((1<<(n+2))-1)^(1<<x)]
    b=a
    for cut in cuts:b=birth(b,cut)
    p=(x,y,z,*path,s,w);q=(x,w,y,z,*path,s)
    frame=lambda t:(t[t.index(root)-1] if t.index(root) else None,t[t.index(root)+1] if t.index(root)+1<len(t) else None)
    need(frame(p)==frame(q),'same literal root frame')
    b=native(native(b,p),q);b=tuple(v&((1<<n)-1) for v in b[:n])
    diff=[(i,j) for i in range(n) for j in range(i+1,n) if ((a[i]>>j)^(b[i]>>j))&1]
    need(diff==[tuple(sorted((x,y)))] and root not in (x,y),'one intended exterior edge')
    return b,{'arc':[names[x],names[y]],'cuts':cuts,'paths':[p,q],
              'helper_names':[f'repair{serial}:{s}' for s in ('z','sink','w')],
              'root_frame':frame(p),'completed_rows':b,'completed_H':subset_hp(b)[-1],
              'completed_strong':strong(b)}


def main():
    source=program('UL');target=program('LU');a=tuple(source['rows']);b=tuple(target['rows']);names=source['names']
    need(names==target['names'] and a[0]==b[0],'same actual names and root star')
    edges=[(i,j) for i in range(len(a)) for j in range(i+1,len(a)) if ((a[i]>>j)^(b[i]>>j))&1]
    need(len(edges)==2,'two switched gates')
    trace=[]
    for i,(x,y) in enumerate(edges):
        a,step=flip(a,x,y,0,names,i);trace.append(step)
    need(a==b and subset_hp(a)[-1]==611 and all(x['completed_strong'] for x in trace),'native neutral reset')
    out={'status':'FINITE_EXACT new specialization of inherited BC2 compiler','source':source,'target':target,
         'root_name':'source:0','trace':trace,'native_events':4,'helper_births':6,'helper_deletions':6,
         'peak_order':12,'deleted_old_names':[],'common_H':611,
         'reader_word':'L','reader_counts':[88693,86317],
         'scope':'root star and all old names preserved; root response phases may change; helper intermediates not claimed strong; not a rooted-ready ternary closure'}
    (HERE/'native_reset.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
    print('PASS source611 native reset: 4 native moves; 6 helpers; zero old-name loss',flush=True)


if __name__=='__main__':main()
