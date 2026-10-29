"""Actual C21 source decomposition and its native Frobenius normalizer."""
from itertools import product
from pathlib import Path
import json

from primary_clock_layers import Ring, inverse_mod, identity, local_smith_with_sources, paired_partition, common_partition
from primary_native_clocks import native_actions, induce, modular_kernel, mul, power, find_perfect, radical_profile
from three_run_glue import three_runs, mm_mod, rref_mod
from involution_root_extensions import run as extension_recurrence
from paley_pencil import local_layers

HERE=Path(__file__).resolve().parent
GATES=0


def need(ok,why):
    global GATES
    GATES+=1
    if not ok:
        raise ValueError(why)


def convolution(a,b,q):
    return [sum(x*y for i,x in enumerate(a) for j,y in enumerate(b) if i+j==r)%q
            for r in range(len(a)+len(b)-1)]


def sector_case(native,sector):
    slots=[i for i,x in enumerate(native['coordinate_names']) if x['sector']==sector]
    clocks={name:[[C[i][j] for j in slots] for i in slots] for name,C in native['clocks'].items()}
    for C in native['clocks'].values():
        need(all(C[i][j]==0 for i in slots for j in range(len(C)) if j not in slots),
             'every full native clock preserves this integral isotype')
    return {'prime':2,'exponents':[native['exponents'][i] for i in slots],
            'orders':[native['orders'][i] for i in slots],'clocks':clocks,'source_slots':slots}


def fixed_points(C,orders):
    q,n=max(orders),len(orders)
    k=q.bit_length()-1
    matrix=[[(q//orders[i])*(C[i][j]-int(i==j)) % q for j in range(n)] for i in range(n)]
    basis,_=modular_kernel(matrix,2,k)
    points={tuple([0]*n)}
    for b in basis:
        b=tuple(x%m for x,m in zip(b,orders))
        multiples={tuple(t*x % m for x,m in zip(b,orders)) for t in range(q)}
        need(len(points)*len(multiples)<=65536,'bounded exact invariant subgroup expansion')
        points={tuple((x+y)%m for x,y,m in zip(a,b,orders)) for a in points for b in multiples}
    need(all(tuple(sum(C[i][j]*v[j] for j in range(n))%orders[i] for i in range(n))==v for v in points),
         'every enumerated descent vector is exactly fixed')
    return sorted(points)


def element_order(v,orders):
    q=max(orders)
    return next(t for t in (2**j for j in range(q.bit_length())) if all(t*x%m==0 for x,m in zip(v,orders)))


def descend(case,degree,partition):
    orders,clocks=case['orders'],case['clocks']
    n,q=len(orders),max(orders)
    C=mul(clocks['translation'],clocks['run_cycle'],orders)
    H=mul(clocks['vertex_frobenius'],clocks['run_swap'],orders)
    gal=H if degree==6 else clocks['vertex_frobenius']
    points=fixed_points(gal,orders)
    es=sorted(partition)
    need(len(points)==2**sum(es),'Galois-fixed subgroup has the complete descended order')
    vectors=[[v for v in points if element_order(v,orders)==2**e] for e in es]
    chosen=None
    for bank in product(*vectors):
        columns=[]
        for v in bank:
            v=list(v)
            for _ in range(degree):
                columns.append(v)
                v=[sum(C[i][j]*v[j] for j in range(n))%orders[i] for i in range(n)]
        W=list(map(list,zip(*columns)))
        if len(rref_mod(W,2)[1])==n:
            chosen=(bank,W)
            break
    need(chosen is not None,'explicit fixed vectors give a full unramified module basis')
    bank,W=chosen
    Wi=inverse_mod(W,q,2)
    need(mul(W,Wi,orders)==identity(n),'descent chart has a named inverse')
    need(all(orders[j]*W[i][j]%orders[i]==0 for i,j in product(range(n),repeat=2)),
         'descent respects every cyclic order')
    converted={name:mul(mul(Wi,A,orders),W,orders) for name,A in clocks.items()}
    # Fixed coordinates are the constant coefficient of each O-module factor.
    positions=[j*degree for j in range(len(es))]
    S=converted['run_swap']
    need(all(S[i][j]==0 for j in positions for i in range(n) if i not in positions),
         'residual involution acts on the actual descent constants')
    small=[[S[i][j] for j in positions] for i in positions]
    return {'fixed_point_count':len(points),'fixed_basis':bank,'chart':W,'inverse':Wi,
            'descended_exponents':es,'descended_swap':small,'converted_clocks':converted}


def involution_embedding(plus,minus):
    a=plus['descended_swap'][0][0]
    es=minus['descended_exponents']
    orders=[2**e for e in es]
    B=minus['descended_swap']
    dual=[[orders[i]*B[j][i]//orders[j]%orders[i] for j in range(len(es))] for i in range(len(es))]
    matches=[]
    for x in product(*(range(q) for q in orders)):
        if element_order(x,orders)!=8:continue
        if all(sum(dual[i][j]*x[j] for j in range(len(es)))%orders[i]==a*x[i]%orders[i] for i in range(len(es))):
            matches.append(x)
    return {'plus_swap':a,'minus_orders':orders,'minus_dual_swap':dual,
            'order8_equivariant_characters':matches,'minimum_attained':bool(matches)}


def fixed_block_certificate(case):
    orders=case['orders']
    need(orders==[2,2,2,2,4,4],'actual trivial-C21 source type')
    F,S=case['clocks']['vertex_frobenius'],case['clocks']['run_swap']
    a,b=[0,0,0,0,1,0],[1,0,0,0,0,1]
    A=list(map(list,zip(a,b)))
    J4=[[0,1],[3,0]]
    swap=[[3,2],[2,1]]
    need(mul(F,A,orders)==A and mul(S,A,orders)==mul(A,swap,orders),
         'literal fixed basis and determinant-minus-one involution')
    need((3*1-2*2)%4==3,'orientation obstruction survives at modulus4')
    fixed=fixed_points(F,orders)
    need({tuple((x*a[i]+y*b[i])%orders[i] for i in range(6)) for x,y in product(range(4),repeat=2)}==set(fixed),
         'entire F-fixed subgroup is the displayed C4 squared')
    E0=[[sum(power(F,t,orders)[i][j] for t in range(3))*pow(3,-1,orders[i])%orders[i]
         for j in range(6)] for i in range(6)]
    E1=[[(int(i==j)-E0[i][j])%orders[i] for j in range(6)] for i in range(6)]
    basis=[]
    for col in zip(*E1):
        need(all(2*x%m==0 for x,m in zip(col,orders)),'nontrivial F-block is elementary')
        if len(rref_mod(basis+[list(col)],2)[1])>len(basis):basis.append(list(col))
    need(len(basis)==4,'four independent elementary complementary coordinates')
    B=list(map(list,zip(*basis)))
    points={tuple(sum(c*b[i] for c,b in zip(cs,basis))%orders[i] for i in range(6)):cs
            for cs in product(range(2),repeat=4)}
    need(len(points)==16,'elementary chart is injective')
    all_points={tuple((u[i]+v[i])%orders[i] for i in range(6)) for u,v in product(fixed,points)}
    need(len(all_points)==256,'fixed and nontrivial blocks recover the entire named source')
    radical=radical_profile(case)
    need(radical['common_radical_log_order']==radical['best_radical_log_order']==2,
         'exact small-block radical lower bound')
    J=radical['best_form']['form']
    JB=mm_mod(mm_mod(list(map(list,zip(*B))),J,4),B,4)
    need(all(x%2==0 for row in JB for x in row),'elementary pairing has the correct denominator')
    JB=[[x//2 for x in row] for row in JB]
    need(len(rref_mod(JB,2)[1])==4,'complementary elementary block has a perfect invariant pairing')
    # Check the canonical doubled fixed-block clock and its pairing literally.
    doubled=[[0]*4 for _ in range(4)]
    hyper=[[0]*4 for _ in range(4)]
    for i,j in product(range(2),repeat=2):
        doubled[i][j]=swap[i][j];doubled[i+2][j+2]=swap[j][i]
        hyper[i][j+2]=int(i==j);hyper[j+2][i]=-int(i==j)%4
    need(mm_mod(mm_mod(list(map(list,zip(*doubled))),hyper,4),doubled,4)==hyper,
         'canonical minimum fixed-block double retains the actual involution')
    return {'fixed_basis':A,'fixed_swap':swap,'fixed_order':16,'nontrivial_basis':B,
            'nontrivial_perfect_form_mod2':JB,'fixed_double_clock':doubled,
            'fixed_double_pairing_mod4':hyper,'minimum_log_index':4,
            'order64_exclusion':'C8^2 requires det(S)=1 mod4; C4^2+C2^2 forces a perfect restriction; C2^6 has no order4 elements',
            'radical_certificate':radical}


def main():
    p,k,q=2,6,64
    eta=0
    for h in range(1,k):
        eta=next(x for x in (eta,eta+2**h) if (x*x+x+2)%2**(h+1)==0)
    bar=(-1-eta)%q
    factors=[[-1,1],[1,1,1],[-1,bar,-eta,1],[-1,eta,-bar,1],
             [1,bar,-1,0,-1,eta,1],[1,eta,-1,0,-1,bar,1]]
    total=[1]
    for f in factors:
        total=convolution(total,f,q)
    need(total==[q-1]+[0]*20+[1],'complete integral X^21-1 factorization')
    need([pow(2,i,21) for i in range(1,7)]==[2,4,8,16,11,1],'primitive21 factors have degree6')
    rows=[sum(1<<((u+s)%7) for s in (1,2,4)) for u in range(7)]
    _,ports,M=three_runs(rows)
    actions=native_actions(ports)
    T,F,S,R=[actions[x] for x in ('translation','vertex_frobenius','run_swap','run_cycle')]
    C=[T[R[i]] for i in range(210)]
    H=[F[S[i]] for i in range(210)]
    need(all(H[C[i]]==C[C[H[i]]] for i in range(210)),'H C H^-1 = C^2 on actual ports')
    def perm_power(P,e):
        a=list(range(210))
        for _ in range(e):a=[P[i] for i in a]
        return a
    need(perm_power(C,21)==perm_power(H,6)==list(range(210)),'full cyclic and Frobenius orders')
    need(perm_power(C,15)==T and perm_power(C,7)==R and perm_power(H,4)==F and perm_power(H,3)==S,
         'two supplied generators recover all four native clocks')
    seen,orbits=set(),[]
    for i in range(210):
        if i in seen:continue
        orbit=[i]
        for _ in range(20):orbit.append(C[orbit[-1]])
        need(len(set(orbit))==21 and C[orbit[-1]]==i,'free joint-cycle source orbit')
        orbits.append(orbit);seen.update(orbit)
    need(len(orbits)==10,'ten complete regular C21 orbits')
    rings=[Ring(p,k,f) for f in factors]
    powers=[]
    for ring in rings:
        t=ring.one if ring.d==1 else (0,1)+(0,)*(ring.d-2)
        ps=[ring.one]
        for _ in range(21):ps.append(ring.mul(ps[-1],t))
        need(ps[-1]==ring.one,'factor clock has order dividing21')
        powers.append(ps)
    chart=[list(row) for ps in powers for row in zip(*ps[:21])]
    ci=inverse_mod(chart,q,p)
    need(mm_mod(chart,ci,q)==identity(21),'complete invertible source chart')
    operators=[[[r.zero]*10 for _ in range(10)] for r in rings]
    for i,j in product(range(10),repeat=2):
        B=[[M[u][v]%q for v in orbits[j]] for u in orbits[i]]
        expected=[[0]*21 for _ in range(21)]
        off=0
        for h,(r,ps) in enumerate(zip(rings,powers)):
            x=r.zero
            for a,c in enumerate(B[0]):x=r.add(x,r.mul(r.scalar(c),ps[-a%21]))
            operators[h][i][j]=x
            block=r.times(x)
            for a,b in product(range(r.d),repeat=2):expected[off+a][off+b]=block[a][b]
            off+=r.d
        need(mm_mod(mm_mod(chart,B,q),ci,q)==expected,'literal native21-block decomposition')
    records,parts,layers=[],[],[0]*k
    for r,A in zip(rings,operators):
        es,U,Ui=local_smith_with_sources(A,r)
        expanded=[[0]*(10*r.d) for _ in range(10*r.d)]
        for i,j in product(range(10),repeat=2):
            B=r.times(A[i][j])
            for a,b in product(range(r.d),repeat=2):expanded[i*r.d+a][j*r.d+b]=B[a][b]
        observed=local_layers(expanded,p,k)
        need(observed==[r.d*es.count(e) for e in range(k)],'independent integer layer check')
        layers=[a+b for a,b in zip(layers,observed)]
        part=sorted((e for e in es if e),reverse=True)
        parts.append(part)
        records.append({'factor_low_coefficients':[x%q for x in r.f],'degree':r.d,
                        'O_module_partition':part,'pivot_valuations':es,
                        'source_row_map_U':U,'source_row_map_inverse':Ui,'Zp_layers':observed})
    need(layers==[165,31,2,3,9,0],'all original dyadic classes recovered')
    targets=[paired_partition(parts[0]),parts[1],common_partition(parts[2],parts[3]),
             common_partition(parts[2],parts[3]),common_partition(parts[4],parts[5]),
             common_partition(parts[4],parts[5])]
    index=sum(r.d*(sum(b)-sum(a)) for r,a,b in zip(rings,parts,targets))
    case={'prime':p,'precision':k,'modulus':q,'source_orbit_chart':chart,
          'source_orbit_chart_inverse':ci,'factors':records}
    data={'translation_orbits':orbits,'ports':ports}
    native=induce(case,data,M,actions)
    sectors=[sector_case(native,i) for i in range(6)]
    small=[find_perfect(sectors[i]['exponents'],sectors[i]['clocks'],2) for i in (0,1)]
    descents={str(i):descend(sectors[i],rings[i].d,parts[i]) for i in (2,3,4,5)}
    embedding=involution_embedding(descents['2'],descents['3'])
    need(embedding['plus_swap']==3 and descents['3']['descended_swap']==[[1,1],[8,5]],
         'extension recurrence uses the exact descended native clocks')
    extension=extension_recurrence()
    need(extension['witness']['depth']==2 and extension['levels'][0]['complete_level'],
         'all smaller equivariant overgroups excluded; minimum cubic completion attained')
    # Simple closed coordinates for the recurrence's attaining group.
    Worders=[4,32];Wswap=[[3,1],[8,5]]
    source_swap=descents['3']['descended_swap']
    injection=[[2,0],[0,2]]
    need(mul(Wswap,injection,Worders)==mul(injection,source_swap,Worders),
         'source C2+C16 embeds with its original involution')
    need(power(Wswap,2,Worders)==identity(2),'larger cubic involution closes')
    need(all((6*(3*x+y)+8*x+5*y-3*(6*x+y))%8==0 for x,y in product(range(4),range(32))),
         'primitive character 6x+y preserves the scalar3 clock')
    small_fixed=fixed_block_certificate(sectors[0])
    need(small[1]['perfect'],'remaining self-reciprocal block needs no completion')
    out={'status':'FINITE-EXACT C21 native source modules','joint_cycle':C,'frobenius_normalizer':H,
         'joint_cycle_orbits':orbits,'source_chart':case,'native_clocks':native,
         'ambient_C21_partitions':targets,'minimum_C21_log_index':index,
         'small_self_dual_forms':small,'descents':descents,'cubic_embedding':embedding,
         'fixed_block_completion':small_fixed,
         'cubic_completion':{'source_orders':[2,16],'ambient_half_orders':Worders,
                             'ambient_half_swap':Wswap,'source_injection':injection,
                             'character_mod8':[6,1],'minimum_log_index':18},
         'full_native_minimum_log_index':40,'component_minimum_log_indices':[4,0,18,18],
         'gates':GATES}
    (HERE/'dyadic_joint_cycle.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'partitions':parts,'minimum_C21_log_index':index,
                      'small_forms':[r['perfect'] for r in small],
                      'cubic_embedding':embedding,'full_native_log_index':40,'gates':GATES}),flush=True)


if __name__=='__main__':main()
