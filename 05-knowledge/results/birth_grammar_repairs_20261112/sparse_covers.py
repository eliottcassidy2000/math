"""Independent endpoint-table/convolution compiler, storing only reachable runs."""
from collections import defaultdict

def ends(rows,ports):
    fixed=[v for v in range(len(rows)) if v not in ports];order=fixed+list(ports)
    n=len(rows);k=len(fixed);t=len(ports);bound=k+2;zero=(0,)*t
    states={};out=[{} for _ in range(1<<k)]
    for v in range(n):
        c=[0]*t
        if v>=k:c[v-k]=1
        states[(1<<v if v<k else 0,tuple(c),v)]=1
    while states:
        nxt=defaultdict(int)
        for (mask,c,v),number in states.items():
            vals=out[mask].setdefault(c,[0]*n);vals[v]+=number
            for w in range(n):
                if not rows[order[v]]>>order[w]&1:continue
                if w<k:
                    if mask>>w&1:continue
                    nxt[(mask|(1<<w),c,w)]+=number
                elif c[w-k]<bound:
                    cc=list(c);cc[w-k]+=1;nxt[(mask,tuple(cc),w)]+=number
        states=nxt
    return order,out

def kernel(rows,ports,selected):
    order,end=ends(rows,ports);n=len(rows);full=len(end)-1
    reverse=[sum(1<<b for b in range(n) if rows[b]>>a&1) for a in range(n)]
    _,start=ends(reverse,ports);sel=[i for i,v in enumerate(order) if v in selected]
    D=defaultdict(int);U=defaultdict(int)
    for c,vals in end[full].items():D[c]+=sum(vals)
    for c,vals in start[full].items():U[c]+=sum(vals[i] for i in sel)
    for mask in range(full+1):
        left=[(a,sum(vals),sum(vals[i] for i in sel)) for a,vals in end[mask].items()]
        right=[(b,sum(vals[i] for i in sel)) for b,vals in start[full^mask].items()]
        for a,h,e in left:
            for b,s in right:
                if not s:continue
                c=tuple(x+y for x,y in zip(a,b));D[c]-=e*s;U[c]+=h*s
    return ({c:v for c,v in D.items() if all(c) and v},{c:v for c,v in U.items() if all(c) and v})
