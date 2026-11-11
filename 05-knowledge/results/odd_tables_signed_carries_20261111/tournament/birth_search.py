"""Small complete source/port-restoration control at the actual133/135 pair."""
from pathlib import Path
from tempfile import TemporaryDirectory
from collections import defaultdict
import subprocess,json
from probe import construct,state,end_dp

HERE=Path(__file__).resolve().parent

def child(rows,cut):
    n=len(rows)
    return [r|(0 if cut>>i&1 else 1<<n) for i,r in enumerate(rows)]+[cut]

def preserving_cuts(rows):
    n=len(rows);full=(1<<n)-1
    dp=end_dp(rows)
    converse=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
    sp=end_dp(converse);E=dp[full];S=sp[full];H=sum(E)
    W=[[0]*n for _ in range(n)]
    for A in range(1,full):
        B=full^A
        for i,x in enumerate(dp[A]):
            if x:
                for j,y in enumerate(sp[B]):
                    if y:W[i][j]+=x*y
    good=[]
    for O in range(full+1):
        I=[i for i in range(n) if not(O>>i&1)]
        J=[j for j in range(n) if O>>j&1]
        value=sum(E[i] for i in I)+sum(S[j] for j in J)+sum(W[i][j] for i in I for j in J)
        if value==H:good.append(O)
    return good

def shared(a,b):
    return [[i,j,s,t] for i,s in enumerate(a) for j,t in enumerate(b) if s[1:]==t[1:]]

def scc(rows):
    n=len(rows);reach=[x|(1<<i) for i,x in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    unseen=set(range(n));components=[]
    while unseen:
        i=min(unseen)
        c=[j for j in sorted(unseen) if reach[i]>>j&1 and reach[j]>>i&1]
        components.append(c);unseen.difference_update(c)
    return components

def main():
    sources=[construct((2,4,1),cut,B) for cut,B in ((3,19),(1,39))]
    roots=[len(r)-2 for r in sources]
    allroots=[[state(rows,j) for j in range(len(rows))] for rows in sources]
    preserved=[preserving_cuts(r) for r in sources]
    padded_matches=[]
    for a in preserved[0]:
        ra=child(sources[0],a);sa=[state(ra,j) for j in range(len(ra))]
        for b in preserved[1]:
            rb=child(sources[1],b);sb=[state(rb,j) for j in range(len(rb))]
            for x in shared(sa,sb):padded_matches.append({'cuts':[a,b],'match':x})
    with TemporaryDirectory(prefix='odd-tables-birth-') as td:
        exe=Path(td)/'birth_probe'
        subprocess.run(['clang++','-O3','-std=c++17',str(HERE/'birth_probe.cpp'),'-o',str(exe)],check=True)
        responses=[]
        for rows,r in zip(sources,roots):
            data=' '.join(map(str,[len(rows),r]+rows))+'\n'
            lines=subprocess.check_output([str(exe)],input=data,text=True).splitlines()
            responses.append([list(map(int,x.split())) for x in lines])
    lookup=defaultdict(list)
    for x in responses[1]:lookup[tuple(x[2:])].append(x)
    pairs=[]
    for lo in responses[0]:
        for hi in lookup.get(tuple(lo[2:]),[]):
            pairs.append({'lower_birth':lo,'upper_birth':hi,'signed_width':hi[1]-lo[1]})
    controls=[]
    for x in pairs:
        lo=x['lower_birth'];hi=x['upper_birth']
        a=child(sources[0],lo[0]);b=child(sources[1],hi[0])
        assert state(a,roots[0])==tuple(lo[1:])
        assert state(b,roots[1])==tuple(hi[1:])
        controls.append({'rows':[a,b],'roots':roots,'states':[lo[1:],hi[1:]],'strong_components':[scc(a),scc(b)]})
    receipt={'scope':'All roots of the two named sources; all H-preserving one-birth cuts and all roots afterward; all512/1024 retained-last-root one-birth cuts with unrestricted new H. No claim about other sources or multiple births.', 'source_rows':sources,'source_roots':roots,'source_all_root_states':allroots,'source_strong_components':[scc(x) for x in sources],'source_shared_ports':shared(*allroots),'H_preserving_cuts':preserved,'H_preserving_birth_shared_ports':padded_matches,'all_common_port_birth_pairs':pairs,'positive_pair_direct_controls':controls}
    (HERE/'birth_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'source_shared_ports':receipt['source_shared_ports'],'H_preserving_cuts':preserved,'all_common_port_birth_pairs':pairs},indent=2))

if __name__=='__main__':main()
