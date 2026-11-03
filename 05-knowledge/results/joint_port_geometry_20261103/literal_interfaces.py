"""Independent graph replay of replacement and persistent-spare cap operators."""
import json
import sympy as s
from return_probe import HERE,junction,observer,need
import effects as e
from returns import cycle
from persistent_port import macro as persistent

def vector(T,ports):
    V=observer(T,ports)
    return s.Matrix([[int(V.get(k,0)) for k in e.STATES]])

def main():
    witnesses=json.loads((HERE/'rank_transfer.json').read_text())['witnesses']
    source=next(x for x in witnesses if len(x['rows'])==5 and x['ports'][1]!=x['provenance'][3] and x['ports'][1]<3)
    T=tuple(source['rows']);r,b=source['ports'];n=len(T);v=vector(T,(r,b));h,t,q,F=junction.profile(T,r)
    need(t==q and F==h,'literal balanced sealed source')
    tests=[]
    for m in range(4):
        inc=((1<<r) if m&1 else 0)|((1<<b) if m&2 else 0)
        U=junction.ear(T,inc);U=junction.ear(U,1<<r);U=junction.ear(U,1<<(n+1))
        out=(n+1,b);M=e.dense((0,m,2))*e.macro(0)
        need(vector(U,out)==v*M,'literal J'+str(m)+' matches exact cap transport')
        need(junction.strong(U),'joint return endpoint strong')
        tests.append({'operation':'J'+str(m),'rows':U,'ports':out,'H':junction.H(U)})
    for j in (3,5):
        U,rr,pp,hh,tt,receipt=cycle(T,r,(r,b),h,t,j)
        replacement=(rr,next(p for p in pp if p!=rr));fixed=(rr,b)
        need(vector(U,replacement)==v*e.macro(j),'actual native replacement macro C'+str(j))
        need(vector(U,fixed)==v*persistent(j),'same physical endpoint with persistent second mark C'+str(j))
        need(b not in replacement and b in fixed,'different live-port policies are explicit')
        tests.append({'operation':'C'+str(j),'rows':U,'replacement_ports':replacement,'persistent_ports':fixed,
                      'root_profile':junction.profile(U,rr),'native_trace':receipt['native_trace'],'clocks':receipt['clocks']})
    (HERE/'literal_interfaces.json').write_text(json.dumps({'status':'FINITE_EXACT_LITERAL_NATIVE_REPLAY',
        'source':source,'tests':tests},indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: J0..J3 literal graphs and C3/C5 native replays match both live-port policies',flush=True)
if __name__=='__main__':main()
