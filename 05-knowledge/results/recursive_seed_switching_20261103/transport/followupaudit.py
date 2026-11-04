"""Independent follow-up native interpreter and field-codec audit.

Imports only the earlier elementary audit, never a compiler or producer helper.
Hamilton counts and17 continuations are independently computed by the sibling
small-feedback engine; this interpreter checks their exact named graph binding.
"""
from pathlib import Path
from itertools import combinations
import hashlib,importlib.util,json
HERE=Path(__file__).resolve().parent;PACKET=HERE.parent;RESULTS=PACKET.parent
ELEMENTARY=RESULTS/'zero_ready_component_returns_20261103/component/audit_native.py'
spec=importlib.util.spec_from_file_location('elementary_followup_audit',ELEMENTARY);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)

def strong(rows):
    for r in range(len(rows)):
        seen=1<<r
        while True:
            nxt=seen
            for i in range(len(rows)):
                if seen>>i&1:nxt|=rows[i]
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<len(rows))-1:return False
    return True

def field_codec(g):
    def mul(a,b):
        terms=[0]*5
        for i,x in enumerate(a):
            for j,y in enumerate(b):terms[i+j]+=x*y
        for j in (4,3):terms[j-1]-=terms[j];terms[j-3]-=terms[j]
        return tuple(x%7 for x in terms[:3])
    def power(a,n):
        result=(1,0,0)
        for _ in range(n):result=mul(result,a)
        return result
    h,x,y=g;z=tuple(t%7 for t in g);orbit=[z,power(z,7),power(z,49)]
    A.need(power(z,343)==z,'current cubic-field codec closes')
    alpha,beta,charge=(x+4*(h-y))%7,y%7,(h-y)%7
    return dict(observer_mod7=list(z),F343_vector_codec=list(map(list,orbit)),
        F343_codec_degree=1 if power(z,7)==z else 3,F49_charge=[alpha,beta,charge],
        F49_frobenius=[(alpha+4*beta)%7,(-beta)%7,charge])

def main():
    inp=PACKET/'selector/followup_named.json';record=HERE/'followup.json';data=json.loads(inp.read_text());full=json.loads(record.read_text());source=data['source']
    repair=json.loads((PACKET/'selector/named_repair.json').read_text())['target']
    parent_audit=json.loads((HERE/'audit_native.json').read_text())
    A.need(source==repair and parent_audit['status']=='PASS' and parent_audit['target_H']==source['H'],'same actually audited repaired parent')
    A.need(full['status']=='FINITE_EXACT_NAMED_FOLLOWUP' and full['input_sha256']==hashlib.sha256(inp.read_bytes()).hexdigest(),'complete follow-up receipt binds exact input')
    A.need(full['source_names']==source['names'] and full['source_rows']==source['rows'],'explicit actual parent rows/names retained')
    A.need(len(full['programs'])==len(full['records'])==len(data['targets'])==2,'exactly two next graphs')
    targets={t['H']:t for t in data['targets']};counts={t['H']:t for t in full['records']};reports=[]
    A.need(set(targets)=={41506343,41506345}==set(counts),'declared two ternary children')
    for program in full['programs']:
        target=targets[program['target_H']];named=counts[target['H']];names=list(source['names']);rows=list(source['rows']);A.valid(rows)
        A.need(program['parent_H']==source['H']==13835447 and target['H']==3*source['H']+target['offset'],'actual common parent and numerical labels')
        source_digest=A.sha([names,rows]);old_names=set(names)
        A.need(target['names'][:len(names)]==names and target['new_names']==target['names'][len(names):],'all old physical names retained in order')
        A.need([z['name'] for z in program['permanent_births']]==target['new_names'],'exact permanent birth names/order')
        for z in program['permanent_births']:
            A.need(z['name'] not in names,'fresh permanent physical name');rows=A.add_vertex(rows,z['out_mask']);names.append(z['name'])
        A.need(names==target['names'],'literal permanent names')
        A.need(A.differences(rows,target['rows'],names)==target['old_arc_differences'],'complete old-arc discrepancy list after births')
        A.need([z['changed_pair'] for z in program['changes']]==target['old_arc_differences'],'every listed arc repaired once in order')
        guards=0;helpers=[]
        for z in program['changes']:
            A.need(z['helper'] not in names,'fresh temporary helper');helpers.append(z['helper']);local=names+[z['helper']]
            A.need(set(z['changed_pair'])<=old_names,'only an original-parent pair is retargeted')
            word=z['pair_guard_names'];where={v:i for i,v in enumerate(names)}
            A.need(sorted(word)==sorted(names),'constructed pair word retains every actual name')
            pair=set(z['changed_pair']);A.need(abs(word.index(z['changed_pair'][0])-word.index(z['changed_pair'][1]))==1,'chosen pair is adjacent in constructed word')
            A.need(all(set((u,v))==pair or rows[where[u]]>>where[v]&1 for u,v in zip(word,word[1:])),
                   'constructed word is legal in one of the two single-pair orientations')
            state=A.add_vertex(rows,z['helper_out_mask']);A.need(len(z['paths'])==2,'exactly two native guards per gadget')
            for word in z['paths']:state=A.native(state,local,word);guards+=1
            after=[row&((1<<len(names))-1) for row in state[:-1]]
            A.need(A.differences(rows,after,names)==[z['changed_pair']],'literal one-arc endpoint')
            A.need(hashlib.sha256(json.dumps(after,separators=(',',':')).encode()).hexdigest()==z['after_delete_sha256'],'literal intermediate matrix digest')
            A.need(not(set(z['changed_pair']) & {39,50,51}),'root and surviving two-prefix incidences unchanged at gadget endpoints')
            rows=after
        A.need(rows==target['rows']==named['rows'] and names==named['names'],'exact compiled and counted full named endpoint')
        A.valid(rows);A.need(strong(rows),'strong full endpoint')
        loc={v:i for i,v in enumerate(names)};arrow=lambda a,b:(rows[loc[a]]>>loc[b])&1
        bn=target['component_names'];B=[[arrow(a,b) for b in bn] for a in bn]
        component=[sum(bit<<j for j,bit in enumerate(row)) for row in B]
        A.need(component==target['component_rows'] and strong(component),'literal strong component')
        A.need(target['prefix']==[50,51] and target['mark']==4 and target['root']==39 and target['spare']==6,'actual active frame')
        A.need(arrow(50,51) and all(arrow(a,b) for a in (50,51) for b in bn),'declared prefix dominates component')
        A.need(arrow(4,39) and all(arrow(39,v) for v in names if v not in (4,39)),'outer root sole in-port4')
        A.need(31 in bn and len({39,4,31,6})==4,'old mark31 and spare6 retained as distinct names')
        A.need([p['ports'] for p in named['joint_pairs']]==[[39,6],[31,6],[4,6],[50,51]],'all four requested ordered pairs')
        charts={}
        for pair in named['joint_pairs']:
            A.need(pair['H']==named['H']==target['H'] and pair['P2']==named['P2'],'same full graph H/P2 in every pair')
            A.need(len(pair['observer17'])==17 and pair['observer17'][0]==target['H'],'17 continuation count baseline')
            A.need(pair['ordered_W2'][0]==pair['ordered_W2'][3]==0,'diagonal disjoint-path endpoints vanish')
            for v,fields in zip(pair['ports'],pair['port_data']):
                A.need(len(fields)==5 and all(x>=0 for x in fields),'endpoint/deletion data type')
                end,start,deletion,E2,start2=fields;A.need((start2-end)%2==0,'integral golden data')
                g=[target['H'],(start2-end)//2,(start2+end)//2]
                chart=dict(golden=g,charts=field_codec(g))
                if str(v) in charts:A.need(charts[str(v)]==chart,'shared spare has same chart in every ordered pair')
                charts[str(v)]=chart
        A.need(charts==named['named_port_charts'] and charts['39']['golden']==named['golden'],'independent named field-codec audit')
        A.need(named['root_fields'][3]==named['H'],'zero-ready root response')
        reports.append(dict(parent_H=source['H'],target_H=target['H'],order=len(names),names=names,
            actual_parent_graph_sha256=source_digest,actual_target_graph_sha256=A.sha([names,rows]),
            permanent_births=len(program['permanent_births']),native_complements=guards,temporary_births_and_deletions=len(helpers),
            deleted_old_names=[],peak_order=len(names)+1,root=39,active_mark=4,retained_old_mark=31,spare=6,prefix=[50,51],
            joint_pairs=[p['ports'] for p in named['joint_pairs']],named_port_charts=charts))
    A.need(sum(x['native_complements'] for x in reports)==86 and sum(x['permanent_births'] for x in reports)==7,'bounded follow-up totals')
    control_keys=set()
    for n,code,a,b,cut,p,q,word in full['direct_guard_controls']:
        pairs=list(combinations(range(n),2));rows=[0]*n
        for j,(u,v) in enumerate(pairs):
            x,y=(u,v) if code>>j&1 else (v,u);rows[x]|=1<<y
        A.valid(rows);A.need((n,code,a,b) not in control_keys,'distinct direct-guard control');control_keys.add((n,code,a,b))
        A.need(sorted(word)==list(range(n)) and abs(word.index(a)-word.index(b))==1,'small pair word type')
        A.need(all({u,v}=={a,b} or rows[u]>>v&1 for u,v in zip(word,word[1:])),'small constructed one-orientation word')
        state=A.add_vertex(rows,cut)
        state=A.native(state,list(range(n+1)),p);state=A.native(state,list(range(n+1)),q)
        after=[row&((1<<n)-1) for row in state[:n]]
        A.need(A.differences(rows,after,list(range(n)))==[[a,b]],'independent small direct-gadget endpoint')
    wanted={(n,code,a,b) for n in (2,3,4) for code in range(1<<(n*(n-1)//2)) for a,b in combinations(range(n),2)}
    A.need(control_keys==wanted and len(wanted)==410,'complete small direct-guard universe')
    A.need(A.native([6,4,0],[0,1,2],[0,1,2])==[2,4,1] and A.native([2,4,1],[0,1,2],[0,1,2])==[6,4,0],'positive and involution controls')
    for hostile in ([2,1,0],[0,1,1]):
        rejected=False
        try:A.native([6,4,0],[0,1,2],hostile)
        except RuntimeError:rejected=True
        A.need(rejected,'illegal/duplicate guard rejected')
    pins=[inp,record,PACKET/'selector/named_repair.json',HERE/'audit_native.json',ELEMENTARY,Path(__file__)]
    print(json.dumps(dict(status='PASS',scope='Independent native and named-field-codec audit of two next endpoints only. No compiler import, no endpoint H recount in this interpreter; exact small-feedback counts are bound by the sibling receipt.',
        actual_parent_H=source['H'],records=reports,native_complements=86,permanent_births=7,temporary_births_and_deletions=43,
        direct_guard_controls=410,positive_hostile_controls=4,source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}),indent=2,sort_keys=True))
if __name__=='__main__':main()
