"""Positive341-clock controls and one literal H25/C3 meeting."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess,json
from probe import g,rev,N,new_pair,qstate,construct

HERE=Path(__file__).resolve().parent
def D19(n):return 25*g(n)-3*g(19*n)+16*n+176
def F19(n):return D19(n)+22*(n&1)-18
def lower(n):return 25*g(n)+244*n+186
def prepare(p):
    qs=[128*p+27,128*p+28]
    q=next(q for q in qs if F19(q)%4==0)
    assert all(F19(t)<0 for t in qs)
    return q,-F19(q)//4
def complete(q,m):
    power=1024**(m+1)
    return q*power+3*(power-1)//1023
def old_pair(n):
    B=2*n+1
    U=147*(1<<B.bit_length())-25*rev(B)
    return [(lower(n)+j,25,U,25) for j in (0,2)]
def main():
    for p in range(128):
        assert D19(1024*p+3)==F19(p)
        assert F19(1024*p+3)==F19(p)+4
        a=25*g(p)-3*g(19*p+4)-1068*p+22*(p&1)-44
        b=25*g(p)-3*g(19*p+4)-1576*p+22*(p&1)-182
        assert F19(128*p+27)==a and F19(128*p+28)==b
        assert a<0 and b<0 and (a-b)%4==2
    completed=[]
    for p in range(4):
        q,m=prepare(p);n=complete(q,m)
        assert D19(n)==0 and lower(n)==N(19*n)
        for a,b in zip(old_pair(n),new_pair(19*n)):
            diff=[x-y for x,y in zip(qstate(a),qstate(b))]
            assert diff[0]==0 and min(diff[1:])>0
        completed.append({'prefix':p,'q':q,'actuator_repeats':m,'total_blocks':m+1,'address_bits':n.bit_length(),'target_bits':lower(n).bit_length()})
    n=890;B=2*n+1;C=19*n
    oldcores=((50,8,3,5,46,14),(50,40,35,5,46,8))
    oldrows=[construct(c,cut,B) for c,cut in zip(oldcores,(38,35))]
    newrows=[construct((2,4,1),cut,b) for cut,b in ((3,2*C+1),(1,4*C+3))]
    actual=[]
    with TemporaryDirectory(prefix='odd-tables-341-') as td:
        exe=Path(td)/'literal_state'
        subprocess.run(['clang++','-O3','-std=c++17',str(HERE/'literal_state.cpp'),'-o',str(exe)],check=True)
        for rows,root,mode in [(r,len(r)-1,1) for r in oldrows]+[(r,len(r)-2,0) for r in newrows]:
            data=' '.join(map(str,[len(rows),root,mode]+rows))+'\n'
            actual.append(tuple(map(int,subprocess.check_output([str(exe)],input=data,text=True).split())))
    assert actual==old_pair(n)+new_pair(C)
    receipt={'scope':'Exact341-actuator/preparation controls on128 inputs; four completed prefixes; literal endpoint DP atn890 on18/18/20/21 vertices. New U uses the proved one-out-neighbour identity, old U is counted independently.', 'identity_inputs':128,'completed_prefixes':completed,'literal_meeting':{'n':n,'old_B':B,'new_C':C,'rows':oldrows+newrows,'states':actual,'target_pair':[lower(n),lower(n)+2]}}
    (HERE/'composite_clock_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'pass','completed_prefixes':completed,'literal_states':actual},indent=2))
if __name__=='__main__':main()
