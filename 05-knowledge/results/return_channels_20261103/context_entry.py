"""Exact counterexamples to right-cancellation of observer-state membership."""
from itertools import permutations
import json
from common import *

def main():
    tested=0;absent=[];rescues=[]
    for h in range(1,16,2):
        for x in range(2*h+1):
            for y in range(2*h+1):
                v=(h,x,y)
                if not physical_cone(v) or decode_state(v) is not None:continue
                absent.append(v)
                for l,m in MATS.items():
                    tested+=1;z=step(v,m);word=decode_state(z)
                    if word is not None:
                        need(vector(word)==z and decode_state(v) is None,'rescued state has no actual predecessor')
                        rescues.append({'formal_state':v,'letter':l,'actual_output_state':z,'actual_word':word})
    need(len(absent)==942 and tested==61230 and len(rescues)==49,'complete stated contextual entry universe')
    need(decode_state((3,0,3)) is None and step((3,0,3),MATS['a'])==vector(('q1',))==(3,3,3),'small cancellation hostile')
    def literal_state(rows,port):
        h=e=d=0
        for p in permutations(range(len(rows))):
            bad=[i for i,(u,v) in enumerate(zip(p,p[1:])) if not(rows[u]>>v&1)]
            if not bad:h+=1;e+=p[-1]==port
            elif len(bad)==1 and p[bad[0]+1]==port:d+=1
        C=h+d;need((C-e)%2==0,'physical parity')
        return (h,(C-e)//2,(C+e)//2)
    outside=(10,12,9,0);extended=(26,28,25,0,8);representative=actual_target(('q1',))
    need(literal_state(outside,3)==(3,0,3),'actual outside-state source')
    need(literal_state(extended,4)==literal_state(representative,3)==(3,3,3),'same observer, different source orders')
    queue=[((),(1,0,0))];fibre=[]
    for word,v in queue:
        if v==(3,3,3):fibre.append(word)
        for l,m in MATS.items():
            z=step(v,m)
            if z[0]<=3:queue.append((word+(l,),z))
    need(len(queue)==12 and fibre==[('q1',)],'all-order exact small observer fibre')
    result={'status':'PASS','universe':'odd h<=15, x,y<=2h, physical cone, initially absent',
            'absent_states':len(absent),'state_letter_pairs':tested,'rescues':rescues,
            'literal_clock_control':{'outside_rows':outside,'outside_port':3,'extended_rows':extended,
                'extended_port':4,'bank_representative_rows':representative,'bank_port':3,
                'orders':[4,5,4],'states':[(3,0,3),(3,3,3),(3,3,3)],'bank_state_fibre':fibre},
            'scope':'A context can enter the observer-state image; this does not realize its formal predecessor or identify the supplied graph with the returned representative.'}
    (HERE/'context_entry.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','absent_states':len(absent),'pairs':tested,'rescues':len(rescues)}))
if __name__=='__main__':main()
