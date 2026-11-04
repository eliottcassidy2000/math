"""Named joint ports and current field codecs for the actual uniform3H family."""
from pathlib import Path
from hashlib import sha256
import json
from audit_moving import query,chart,bridge,need,birth,strong
HERE=Path(__file__).resolve().parent
def family(m):
    names=[31]+[f'clone_{j}' for j in range(1,m)]+[32,22,2,1,6,8,4,13,5,14,7]+[21]
    p=m+11;T=[0]*(p+1)
    for i in range(p):
        for j in range(i+1,p):T[i]|=1<<j
        if i<m or not (1243>>(i-m))&1:T[p]|=1<<i
        else:T[i]|=1<<p
    return tuple(T),names
def main():
    path=HERE.parent/'return18979_ports_20261103/grandparent_repair.json';old=json.loads(path.read_text());records=[]
    for m in range(1,8):
        T,names=family(m)
        if m==1:
            perm=[old['source_names'].index(a) for a in names];R=old['repaired_rows'];U=tuple(sum(((R[a]>>b)&1)<<j for j,b in enumerate(perm)) for a in perm)
            need(U==T,'literal inherited repaired source, including every name')
        cut=((1<<m)-1)|(1<<(m+1));V=birth(T,cut)
        for stage,U,labels,ports in [('source',T,names,[31,6,21]),('second_port',V,names+['feedback_q'],[31,6,21,'feedback_q'])]:
            need(strong(U),'actual strong family sample');pairs=[];g=None
            for i in range(len(ports)):
              for j in range(i+1,len(ports)):
                v=query(U,labels.index(ports[i]),labels.index(ports[j]));expected=(1243*(2**m)-1061)*(1 if stage=='source' else 3)
                need(v[0]==expected,'exact source clock and3H identity')
                if i==0 and j==1:
                    need((v[6]-v[2])%2==0,'integer golden vector');g=[v[0],(v[6]-v[2])//2,(v[6]+v[2])//2]
                pairs.append(dict(ports=[ports[i],ports[j]],P2=v[1],port_data=[v[2:7],v[7:12]],ordered_W2=v[12:16]))
            records.append(dict(m=m,stage=stage,rows=U,names=labels,ports=ports,birth_cut=cut if stage=='second_port' else None,H=expected,joint_pairs=pairs,golden=g,charts=chart(g),five_Q4=bridge.five_charts(g)))
    report=dict(status='FINITE_EXACT_NAMED_PORTS_AND_PHASES',inherited_source='return18979_ports_20261103/grandparent_repair.json',source_sha256=sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),m_range=[1,7],fields=['E1','S1','deletion_H','E2','S2'],states=records,
        scope='Named full relations plus all ordered pair responses. Current codecs are not replacements for joint tensors or ancestral clocks.')
    (HERE/'feedback_ports.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS named uniform3H ports and current phase codecs')
if __name__=='__main__':main()
