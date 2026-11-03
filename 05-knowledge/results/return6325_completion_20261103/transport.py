"""Named joint ports and two distinct finite-field charts; no invented clock."""
from return6325 import HERE,need,observer,phases,values
from observer import records
import json
def mul(a,b):
    c=[0]*5
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    for i in (4,3):c[i-1]-=c[i];c[i-3]-=c[i]
    return tuple(x%7 for x in c[:3])
def power(a,k):
    out=(1,0,0)
    while k:
        if k&1:out=mul(out,a)
        a=mul(a,a);k//=2
    return out
def chart(g):
    h,x,y=g;z=tuple(v%7 for v in g);orbit=[z,power(z,7),power(z,49)]
    need(power(z,343)==z,'literal cubic-field Frobenius closes')
    a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    need(((b+c)%7,(a-4*c)%7,b)==z,'linear change from F49/charge coordinates')
    return dict(observer_mod7=z,F343_vector_codec=orbit,F343_codec_degree=next(i for i in (1,3) if power(z,7**i)==z),
        F49_charge=[a,b,c],F49_frobenius=[(a+4*b)%7,(-b)%7,c])
def main():
    out=[];cache={}
    for file in ('return6325.json','return6325_plus2.json'):
        d=json.loads((HERE/file).read_text());names=d.get('old_names',d['names'])
        for stage,T,ns,r in (('source',d['source'],names,d['root']),('target',d['target'],d['names'],d['root']),('next',d['next_rows'],d['next_names'],d['next_root'])):
            b=ns.index(6);T=tuple(T);key=(T,r,b)
            if key not in cache:
                p=phases(T,r);c=chart(p['golden']);need(c['F49_charge']==list(p['golden_phase']),'inherited chart agrees')
                cache[key]=dict(phases=p,charts=c,joint_ports=[ns[r],ns[b]],joint_response=records(observer(T,(r,b))))
            out.append(dict(branch=file,stage=stage,rows=T,names=ns,root=r,spare=b,**cache[key]));print('joint transport',file,stage,flush=True)
    # The old 233/223 staircase is a different named source, recovered literally.
    T=(252,169,242,164,202,144,42,64);square=[]
    for mask in range(4):
        U=list(T)
        for bit,(a,b) in enumerate(((0,3),(2,3))):
            if mask>>bit&1:U[a]^=1<<b;U[b]^=1<<a
        row=[]
        for r in range(8):row.append(phases(tuple(U),r)['golden'])
        square.append(row)
    mixed=[[square[3][r][j]-square[1][r][j]-square[2][r][j]+square[0][r][j] for j in range(3)] for r in range(8)]
    need(all(v[0]==42 for v in mixed),'inherited scalar mixed response')
    inherited=json.loads((HERE.parent/'adaptive_four_returns_20261103/audit.json').read_text());ancestral=[]
    for row in inherited['literal']:
        a=row['ancestry'];target=(4,6,1) if a['family']=='straight' else (6,4,1)
        need(mul((4,4,1),a['F343'])==target,'frozen external ancestral clock control')
        ancestral.append(dict(source=row['source'],ancestry=a))
    result=dict(status='FINITE_EXACT_NAMED_JOINT_TRANSPORT',states=out,
        ancestral_controls=ancestral,ancestral_scope='unchanged earlier 1457/1459/1463/1465 source history; not a claimed 6325 interval-growth ancestry',
        F343_scope='new F7-vector codec of the current three-coordinate observer; NOT inherited spectral-clock transport',
        F343_modulus=[1,0,1,1],F343_frobenius_columns=[power(z,7) for z in ((1,0,0),(0,1,0),(0,0,1))],
        recovered_staircase=dict(rows=T,edges=[[0,3],[2,3]],golden_by_state_and_root=square,mixed=mixed,mixed_mod7=[[x%7 for x in v] for v in mixed]))
    (HERE/'transport.json').write_text(json.dumps(result,indent=2)+'\n')
    print('recovered233/223 mixed golden',mixed,flush=True)
if __name__=='__main__':main()
