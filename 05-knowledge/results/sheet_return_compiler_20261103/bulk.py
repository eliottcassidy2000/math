"""Exact paired-boundary transfer for the literal native escape lifts."""
from collections import defaultdict, deque
import json
from pathlib import Path
from escape_lifts import paths, native, birth, tables, strong
HERE=Path(__file__).resolve().parent


def need(ok,message):
    if not ok:raise ValueError(message)


def target(m,n,kind):
    T,ad,pp=paths(m,n)
    p=dict(pp)[kind]
    incoming={(0,0),(0,m-1),(2,0),(5,n-2)} if kind=='straight' else {(0,1),(1,0),(2,0)}
    mask=sum(1<<ad.index(x) for x in incoming)
    return birth(native.native(T,p),mask),ad+[('new',0)]


def graph(m,n,kind,probe=None):
    T,ad=target(m,n,kind)
    if probe:
        from return_program import golden
        block=golden.literal((probe,));d=len(block);root=ad.index((0,2))
        rows=list(block)+[(row<<d)|((1<<d)-1) for row in T]
        rows[d-1]|=1<<(d+root);rows[d+root]&=~(1<<(d-1))
        T=tuple(rows);ad=[('probe',i) for i in range(d)]+ad
    N=len(T)
    # A universal source and sink turn all H paths into one paired query.
    T=tuple(row|(1<<(N+1)) for row in T)+(((1<<N)-1)|(1<<(N+1)),0)
    ad=ad+[('source',0),('sink',0)]
    variable=[i for i,a in enumerate(ad) if a[0]==5]
    outside=[i for i in range(len(T)) if i not in variable]
    # Index one is a generic interior vertex at n>=4; the penultimate
    # vertex has exceptional newborn incidence on the straight branch.
    middle=variable[1]
    left=[i for i in outside if T[i]>>middle&1]
    right=[i for i in outside if not T[i]>>middle&1]
    order=left+list(reversed(variable))+right
    # Each fixed vertex gets a stable name; 100 is the moving bulk endpoint.
    names={i:j for j,i in enumerate(outside)}
    names.update({variable[n-1]:101,variable[n-2]:100,variable[1]:102,variable[0]:103})
    return T,ad,order,names,len(left)


def advance(states,v,back,active,start=False,end=False):
    out=defaultdict(int)
    for (done,components),weight in states.items():
        lefts=[None] if start else [None]+[i for i,(h,t) in enumerate(components) if t==-1 or t>=0 and t not in back]
        rights=[None] if end else [None]+[i for i,(h,t) in enumerate(components) if h>=0 and h in back]
        for a in lefts:
          for b in rights:
            if a is not None and a==b:continue
            h=components[a][0] if a is not None else -2 if start else v
            t=components[b][1] if b is not None else -2 if end else v
            pairs=[pair for i,pair in enumerate(components) if i not in (a,b)]
            finished=done
            if h<0 and t<-1:
                if h!=t or done:continue
                finished=1
            else:pairs.append((h,t))
            if any(h>=0 and h not in active for h,t in pairs):continue
            pairs=[(h,-1 if t>=0 and t not in active else t) for h,t in pairs]
            out[finished,tuple(sorted(pairs))]+=weight
    return dict(out)


def blueprint(m,kind,probe=None):
    T,ad,order,names,L=graph(m,4,kind,probe)
    steps=[]
    for j,old in enumerate(order):
        previous=order[:j];future=order[j+1:]
        back={names[x] for x in previous if T[old]>>x&1}
        active={names[x] for x in order[:j+1] if any(T[y]>>x&1 for y in future)}
        steps.append((names[old],back,active,ad[old][0]=='source',ad[old][0]=='sink'))
    initial={(0,()):1}
    for args in steps[:L+2]:initial=advance(initial,*args)
    persistent=steps[L+1][2]-{100}
    def move(state):
        moved=advance({state:1},104,{100},persistent|{104})
        return {(done,tuple(sorted((100 if h==104 else h,100 if t==104 else t) for h,t in comp))):w
                for (done,comp),w in moved.items()}
    def finish(state):
        values={state:1}
        for args in steps[L+2:]:values=advance(values,*args)
        return values.get((1,()),0)
    states=list(initial);known=set(states);pending=deque(states);edges={}
    while pending:
        s=pending.popleft();edges[s]=move(s)
        for t in edges[s]:
            if t not in known:known.add(t);states.append(t);pending.append(t)
    index={s:i for i,s in enumerate(states)};d=len(states)
    B=[[0]*d for _ in states]
    for s,targets in edges.items():
        for t,w in targets.items():B[index[s]][index[t]]=w
    return dict(m=m,kind=kind,probe=probe,states=states,steps=[(v,sorted(b),sorted(a),s,t) for v,b,a,s,t in steps],
                persistent=sorted(persistent),transfer=B,initial=[initial.get(s,0) for s in states],terminal=[finish(s) for s in states])


def values(data,count):
    B=data['transfer'];row=data['initial'];end=data['terminal'];out=[]
    for _ in range(count):
        out.append(sum(a*b for a,b in zip(row,end)))
        row=[sum(row[i]*B[i][j] for i in range(len(B))) for j in range(len(B))]
    return out


def audit(data):
    from sympy import Matrix
    B=Matrix(data['transfer']);a=Matrix([data['initial']]);b=Matrix(data['terminal'])
    d=B.rows
    powers=[b]
    for _ in range(3):powers.append(B*powers[-1])
    residual=powers[3]-powers[2]-powers[1]-powers[0]
    # Exact finite-dimensional certificate, NOT agreement of initial scalars.
    row=a
    for _ in range(d):
        need((row*residual)[0]==0,'full reachable-space cubic certificate')
        row=row*B
    seq=values(data,8)
    minor=int(Matrix([[seq[i+j] for j in range(3)] for i in range(3)]).det())
    need(minor!=0,'minimal cubic Hankel witness')
    data.update(values_from_n4=seq,hankel=minor,cubic_residual=list(map(int,residual)))
    for n in range(4,9):
        T,_=target(data['m'],n,data['kind'])
        need(tables(T)[0][-1]==seq[n-4],'independent literal subset count')
        need(strong(T),'strong target controls');tables.cache_clear()
    return data


def main():
    records=[]
    for kind in ('straight','revisit'):
        data=blueprint(3,kind)
        print(kind,'built',len(data['states']),'initial values',values(data,5),flush=True)
        data=audit(data);records.append(data)
        print(kind,'states',len(data['states']),'values',data['values_from_n4'],'Hankel',data['hankel'],flush=True)
    (HERE/'bulk.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':main()
