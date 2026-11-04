#!/usr/bin/env python3
import pathlib,json,math,hashlib,sys
HERE=pathlib.Path(__file__).resolve().parent
S=pathlib.Path(sys.argv[1]).resolve()if len(sys.argv)>1 else HERE.parent/'session189-frontier-eight'
def need(x,msg):
    if not x:raise RuntimeError(msg)
def modules(r):
    full=(1<<len(r))-1;out=[]
    for m in range(1,full+1):
        first=(m&-m).bit_length()-1;outside=full^m
        if all(not((r[first]^r[v])&outside)for v in range(len(r))if m>>v&1):out.append(m)
    return out
def changed_nodes(r,t):
    ms=modules(r);mt=set(modules(t));ss=[a for a in ms if not any(a&b and a&~b and b&~a for b in ms)]
    need(set(ss)<=mt,'An original strong module is lost')
    parent={a:min((b for b in ss if a!=b and a&b==a),key=int.bit_count,default=None)for a in ss}
    records=[];changed=[]
    for m in ss:
        ch=sorted(a for a in ss if parent[a]==m)
        if not ch:continue
        v=[(a&-a).bit_length()-1 for a in ch]
        q=[[sum(1<<j for j,w in enumerate(v)if rows[u]>>w&1)for u in v]for rows in(r,t)]
        records.append({'node':m,'children':ch,'source_quotient':q[0],'target_quotient':q[1],'changed':q[0]!=q[1]})
        if q[0]!=q[1]:changed.append(m)
    high=[a for a in changed if not any(a!=b and a&b==a for b in changed)]
    return ms,records,high
def all_profiles(r):
    # Flattened vertex words with a cut marker. This independently counts
    # ordered covers rather than unordered set partitions times h!.
    n=len(r);need(math.factorial(n)*2**max(0,n-1)<2**64,'Packed integer bound')
    N=1<<n;B=1<<64;dp=[[0]*n for _ in range(N)];out=[[1]]
    for mask in range(1,N):
        words=0;vv=mask
        while vv:
            bit=vv&-vv;vv-=bit;v=bit.bit_length()-1;rest=mask^bit
            if not rest:count=1
            else:
                count=0;uu=rest
                while uu:
                    b=uu&-uu;uu-=b;u=b.bit_length()-1
                    count+=dp[rest][u]*(B+int(bool(r[u]>>v&1)))
            dp[mask][v]=count;words+=count
        profile=[0]
        for _ in range(mask.bit_count()):profile.append(words&(B-1));words>>=64
        need(words==0,'Polynomial degree spill');out.append(profile)
    return out
def ends(r):
    n=len(r);N=1<<n;E=[[0]*n for _ in range(N)];incoming=[sum(1<<u for u in range(n)if r[u]>>v&1)for v in range(n)]
    for mask in range(1,N):
        vv=mask
        while vv:
            bit=vv&-vv;vv-=bit;v=bit.bit_length()-1;rest=mask^bit
            if not rest:E[mask][v]=1;continue
            uu=rest&incoming[v];ans=0
            while uu:b=uu&-uu;uu-=b;ans+=E[rest][b.bit_length()-1]
            E[mask][v]=ans
    return E
def H(r):return sum(ends(r)[-1])
def joint(r):
    n=len(r);N=(1<<n)-1;E=ends(r);R=[sum(1<<u for u in range(n)if r[u]>>v&1)for v in range(n)];A=ends(R);C=[[0]*n for _ in range(n)]
    for mask in range(1,N):
        ee=[(j,x)for j,x in enumerate(E[mask])if x];aa=[(j,x)for j,x in enumerate(A[N^mask])if x]
        for j,x in ee:
            for k,y in aa:C[j][k]+=x*y
    return C
def augment(r,a,c):
    n=len(r);out=r[:]+[(1<<(n+1))|(((1<<n)-1)^a),((1<<n)-1)^c]
    for v in range(n):out[v]|=((a>>v&1)<<n)|((c>>v&1)<<(n+1))
    return out
def main():
    need((S/'normal.json').read_bytes()==(S/'optimized.json').read_bytes(),'Producer normal/-O differ')
    producer=json.loads((S/'normal.json').read_text());need(producer['status']=='PASS','Producer has not passed')
    for name,digest in producer['source_pins'].items():need(hashlib.sha256((S/name).read_bytes()).hexdigest()==digest,'Frozen integration source changed: '+name)
    for name,pin in producer['artifacts'].items():need(hashlib.sha256((S/name).read_bytes()).hexdigest()==pin['sha256'],'Frozen integration artifact changed: '+name)
    for pin in json.loads((S/'DEPENDENCIES.json').read_text()):
        path=S.parent/pin['path']
        if not path.exists():path=S.parent/pin['local_folder']/pin['path'].split('/',1)[1]
        need(hashlib.sha256(path.read_bytes()).hexdigest()==pin['sha256'],'Reviewed dependency changed: '+pin['path'])
    for pin in json.loads((S/'INHERITED.json').read_text()):
        data=(S/pin['local']).read_bytes();need(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==pin['git_blob'],'Inherited immutable Git blob changed')
    controls=json.loads((S/'controls.json').read_text());records=[];all_coordinates=0
    for z in controls:
        r=z['source'];v=z['intermediate'];u=z['target'];n=len(r)
        ms,nodes,highest=changed_nodes(r,u)
        need(ms==z['original_modules']and nodes==z['strong_node_records'],'Original module/tree records differ')
        need(set(highest)==set(z['highest_changed_nodes'])=={z['core'],z['other']},'Highest node antichain')
        need(set(z['family'])<=set(ms),'Selected family is not original')
        for i in range(n):
            for j in range(i+1,n):
                toggle=sum(1 for m in z['family']if m>>i&1 and m>>j&1)%2
                need((r[i]>>j&1)^toggle==(u[i]>>j&1),'Wrong original-module XOR')
                if toggle:need(any(m>>i&1 and m>>j&1 for m in highest),'Changed edge outside highest antichain')
        pr,pu,pv=all_profiles(r),all_profiles(u),all_profiles(v)
        need(pr==pu==pv==z['hereditary_profiles_source']==z['hereditary_profiles_target'],'Independent all-subset cut profiles disagree')
        all_coordinates+=sum(len(row)for row in pr)
        a,c=z['query'];vals=[H(augment(q,a,c))for q in(r,v,u)]
        need(vals==z['actual_response_values']and vals[1]==vals[2]and(vals[0]-vals[1])%4==2,'Physical separator/restoration values')
        M=z['core'];N=z['other'];full=(1<<n)-1;outside=full^M;root=(M&-M).bit_length()-1
        for j in range(n):
            if outside>>j&1:need((a>>j&1)==(c>>j&1)==(r[j]>>root&1),'Wrong ambient port incidence')
        need(a&N in(0,N)and c&N in(0,N),'Other node not uniform to physical ports')
        ja=(a&outside)|(M^(c&M));jc=(c&outside)|(M^(a&M))
        dual=[H(augment(q,ja,jc))for q in(r,v,u)]
        need(dual[1]==dual[2]and dual[0]-dual[1]==-(vals[0]-vals[1]),'Exact baseline ambient duality/restoration')
        records.append({'name':z['name'],'old_order':n,'all_subset_profile_coordinates_checked':sum(len(row)for row in pr),'actual_rows':{'source':r,'intermediate':v,'target':u},'highest_nodes':highest,'query':[a,c],'response_values':vals,'core_only_dual_query':[ja,jc],'dual_response_values':dual})
    banks=json.loads((S/'named_banks.json').read_text());basis_count=shared_count=zero_count=0
    for bank in banks:
        rows=bank['rows'];n=len(rows);masks=bank['physical_masks'];ix={a:i for i,a in enumerate(masks)};co=bank['coefficients'];vals=bank['values'];L=len(masks)
        need(masks==[a for a in range(1<<n)if a.bit_count()<=2]and len(co)==L and all(len(z)==L for z in co),'Incomplete full named bank')
        for a in masks:
            for c in masks:
                pred=sum(co[ix[x]][ix[y]]for x in masks if x&a==x for y in masks if y&c==y)
                need(pred==vals[ix[a]][ix[c]],'Full coefficient/value inversion')
                basis_count+=1
        rec=next(x for x in records if x['name']==bank['case']);role=['source','intermediate','target'].index(bank['role'])
        a,c=rec['core_only_dual_query'];pred=sum(co[ii][jj]for ii,x in enumerate(masks)if x&a==x for jj,y in enumerate(masks)if y&c==y)
        need(pred==rec['dual_response_values'][role],'Nonbasis full polynomial versus independently counted dual query')
        for a in masks:
            if a.bit_count()==2:need(co[ix[a]][ix[a]]==0,'Two-shared-name quartic must vanish');zero_count+=1
        for deleted in range(n):
            keep=[v for v in range(n)if v!=deleted];di={v:i for i,v in enumerate(keep)}
            subrows=[sum(1<<j for j,b in enumerate(keep)if rows[a]>>b&1)for a in keep];C=joint(subrows)
            for u in range(n):
                for v in range(u+1,n):
                    if deleted in(u,v):continue
                    a=(1<<deleted)|(1<<u);c=(1<<deleted)|(1<<v);value=C[di[u]][di[v]]+C[di[v]][di[u]]
                    need(co[ix[a]][ix[c]]==co[ix[c]][ix[a]]==value,'Independent shared-name quartic')
                    shared_count+=2
    need(basis_count==34800,'Full six-bank coordinate count')
    pins={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in(S/'PROOF.md',S/'controls.py',S/'controls.json',S/'INHERITED.json',S/'named_banks.json',S/'build_banks.py')}
    result={'status':'PASS_INDEPENDENT_HIGH_NODE_CONTROLS','method':'Backward endpoint polynomial recurrence in cut marker; no producer imports','all_subset_profile_coordinates_per_orientation':all_coordinates,'controls':records,'complete_named_bank_coordinates':basis_count,'independent_shared_quartics':shared_count,'independent_two_shared_zero_quartics':zero_count,'source_pins':pins,'producer_normal_optimized_byte_equal':True,'producer_receipt_sha256':hashlib.sha256((S/'normal.json').read_bytes()).hexdigest(),'reviewed_dependency_pins':producer['reviewed_dependency_pins']}
    (HERE/'frontier_controls_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: original module antichains, all induced-subset complete cover profiles, literal separator/restoration, and exact core-only dual queries.')
if __name__=='__main__':main()
