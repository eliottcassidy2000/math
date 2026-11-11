"""Complete address inverses for incoming paired carriers on the +4 target ray."""
from pathlib import Path
import importlib.util,json
import chart_inverse as inv
import audit_defect as audit

HERE=Path(__file__).resolve().parent
Q=(18,8,3,5,14)
MINUS=(50,8,3,5,46,14)
PLUS=(50,40,35,5,46,8)


def state(rows):
    n=len(rows);root=n-1;full=(1<<n)-1
    starts=audit.starts(rows)
    ends=audit.starts([sum(((rows[w]>>v)&1)<<w for w in range(n)) for v in range(n)])
    H=[1]+[sum(x) for x in starts[1:]]
    U=sum(starts[S][root]*H[full^S] for S in range(1,full+1) if S>>root&1)
    return [H[full],ends[full][root],U,H[full^(1<<root)]]


def main():
    src=inv.RESULTS/'session251_creative_completion_20261111/ternary/compiler.py'
    spec=importlib.util.spec_from_file_location('signed',src)
    oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    families=[('H13_width2',13,32,6,-1,0,None),
              ('H13_width4',13,136,6,-59,1,None),
              ('H25_width2',25,122,16,-31,0,1),
              ('H25_width2or4',25,426,16,-151,1,None)]
    reports=[];literal=[]
    for name,d,c,e,offset,A,parity in families:
        rows=[]
        for k in range(1,31):
            target=102*2**k+19
            bs=sorted(inv.linear_inverse(102*2**k+offset,d,c,e))
            if parity is not None:bs=[B for B in bs if B%2==parity]
            if name.startswith('H25'):bs=[B for B in bs if B>0]
            for B in bs:
                if name.startswith('H13'):
                    lower=oracle.construct(Q,6,A,B);upper=oracle.construct(Q,3,A,B)
                elif A==0:
                    lower=oracle.construct(MINUS,38,A,B);upper=oracle.construct(PLUS,35,A,B)
                else:
                    lower=oracle.construct(PLUS,35,A,B);upper=oracle.construct(MINUS,38,A,B)
                lo=state(lower);hi=state(upper)
                inv.require(hi[0]==target and hi[1:]==lo[1:],'literal pure-H pair')
                inv.require(hi[0]-lo[0] in (2,4),'allowed baseline gap')
                literal.append({'family':name,'k':k,'A':A,'B':B,'lower_state':lo,'upper_state':hi,
                                'lower_rows':lower,'upper_rows':upper,'retained_root':len(lower)-1})
            rows.append({'k':k,'target':target,'all_qualifying_B':bs})
        reports.append({'family':name,'linear_inverse_parameters':{'d':d,'c':c,'e':e,'rhs':'102*2^k'+str(offset)},
                        'left_address':A,'required_B_parity':parity,'rows':rows})
    inv.require(len(literal)==1 and (literal[0]['k'],literal[0]['B'])==(7,97),'reported finite intersection')
    report={'status':'FINITE-EXACT complete address intersections, k1..30; no same-H source operation or uniform return established',
            'scope':'Four inherited H13/H25 pure-gap families; no address or graph-order cutoff, using the proved exact linear-g inverse.',
            'families':reports,'literal_pair_controls':literal,
            'missing_operation':'An actual compatible replacement from the uniform-reset source of H13073 to the H25 lower witness; equal scalar H alone does not provide it.'}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'carry_ray_probe.json').write_text(out,encoding='utf-8')
    print(json.dumps({'status':report['status'],'complete_inverses':120,'literal_hits':len(literal)}))


if __name__=='__main__':main()
