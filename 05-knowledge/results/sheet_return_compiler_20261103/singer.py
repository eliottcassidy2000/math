"""Transport the actual joint-response clock into the older Singer57 chart."""
from phase import HERE,need
import json
from math import gcd


def add(a,b):return tuple((x+y)%7 for x,y in zip(a,b))
def scale(c,a):return tuple(c*x%7 for x in a)
def mul(a,b):
    out=[0]*5
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    for k in (4,3):
        out[k-1]-=out[k];out[k-3]-=out[k];out[k]=0
    return tuple(x%7 for x in out[:3])
ONE=(1,0,0);ZERO=(0,0,0)
def power(a,n):
    v=ONE
    while n:
        if n&1:v=mul(v,a)
        a=mul(a,a);n//=2
    return v
def trace(a):
    b=add(add(a,power(a,7)),power(a,49))
    need(b[1:]==(0,0),'trace lies in base field')
    return b[0]


def run():
    primitive=(0,1,5);powers=[power(primitive,i) for i in range(342)]
    need(len(set(powers))==342,'inherited Singer generator is primitive')
    logs={v:i for i,v in enumerate(powers)}
    roots=[i for i,a in enumerate(powers) if add(add(add(power(a,3),scale(-2,power(a,2))),scale(-1,a)),scale(-2,ONE))==ZERO]
    exponent=min(roots);alpha=powers[exponent]
    need(gcd(exponent,57)==1 and power(alpha,57)==scale(2,ONE),'Singer generator with scalar return2')
    derivative=add(add(scale(3,power(alpha,2)),scale(-4,alpha)),scale(-1,ONE))
    gamma=power(derivative,341)
    aa=add(add(power(alpha,2),scale(-2,alpha)),scale(-1,ONE));bb=add(alpha,scale(-2,ONE))
    def encode(v):return add(add(scale(v[0],aa),scale(v[1],bb)),scale(v[2],ONE))
    for i in range(3):
        v=tuple(int(i==j) for j in range(3));nextv=(v[1],v[2],(2*v[0]+v[1]+2*v[2])%7)
        need(encode(nextv)==mul(alpha,encode(v)),'basis proof of field conjugacy')
        need(trace(mul(gamma,power(alpha,i)))==int(i==2),'basis proof of trace extraction')
    D=[i for i in range(57) if trace(powers[i])==0]
    need(D==[8,15,19,38,48,50,51,56],'literal older trace-zero line')
    inputs=json.loads((HERE/'phase_bulk.json').read_text())['records'];records=[]
    for record in inputs:
        v=tuple(x%7 for x in record['S_initial']);z=encode(v);offset=logs[mul(gamma,z)]%57
        zeros=[];classes=set()
        for r in range(171):
            need(encode(v)==mul(z,power(alpha,r)),'joint-state field intertwining')
            need(trace(mul(gamma,encode(v)))==v[0],'actual ordered-two-bin response is trace')
            classes.add(logs[encode(v)]%57)
            if r<57 and v[0]==0:zeros.append(r)
            v=(v[1],v[2],(2*v[0]+v[1]+2*v[2])%7)
        need(len(classes)==57 and len(zeros)==8,'all projective points and eight zero phases')
        need(zeros==[r for r in range(57) if (offset+exponent*r)%57 in D],'explicit old-line transport')
        records.append(dict(kind=record['kind'],encoded_initial=z,trace_offset=offset,zero_growth_residues_mod57=zeros,
                            zero_n_residues_mod57=sorted((r+4)%57 for r in zeros)))
    out=dict(field='F7[t]/(t^3+t^2+1)',primitive=primitive,new_root_exponents=roots,
             chosen_exponent=exponent,alpha=alpha,gamma=gamma,trace_zero_line=D,records=records)
    (HERE/'singer.json').write_text(json.dumps(out,indent=2)+'\n')
    print('literal Singer transport',out,flush=True)
    return out


if __name__=='__main__':run()
