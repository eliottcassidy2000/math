import json,hashlib,os
root=os.path.dirname(__file__)
primary=os.path.join(root,'..','session7-rational-cycles')
def read(p):return json.load(open(p))
def compare(left,right):
 a=read(left)['results'];b=read(right)['results']
 if [x['b'] for x in a]!=[x['b'] for x in b]:raise ValueError('denominators differ')
 if [x['cycles'] for x in a]!=[x['cycles'] for x in b]:raise ValueError('cycles differ')
 return dict(denominator_count=len(a),cycle_count=sum(len(x['cycles'])for x in a),periods=sorted({c['period'] for x in a for c in x['cycles']}),independent_triangle_count=sum(x['primitive_triangle_candidates'] for x in a))
receipts={
 'b400_unrestricted_smooth_pairs':compare(root+'/independent_results.json',primary+'/denominators_through_400.json'),
 'b400_normal_vs_optimized':compare(root+'/independent_results.json',root+'/independent_results_optimized.json'),
 'b1000_primitive_word_projected_pairs':compare(root+'/projected_results.json',primary+'/denominators_through_1000.json'),
 'b1000_normal_vs_optimized':compare(root+'/projected_results.json',root+'/projected_results_optimized.json'),
 'b1000_universal_pairs':compare(root+'/universal_pair_results.json',primary+'/denominators_through_1000.json'),
 'b1000_author_new_pair_implementation_replay':compare(root+'/author_pair_full_replay.json',primary+'/denominators_through_1000.json')}
if open(root+'/author_large_exponent_replay.json','rb').read()!=open(root+'/author_large_exponent_replay_optimized.json','rb').read():raise ValueError('large-exponent normal/optimized mismatch')
large=read(root+'/author_large_exponent_replay.json')['results'][0]
independent=read(root+'/final_addendum_controls.json')['large_exponent_independent_search_without_escape_guard']
if (large['b'],large['cycles'],large['primitive_triangle_candidates'])!=(independent['b'],independent['cycles'],independent['primitive_triangle_candidates']):raise ValueError('independent large-exponent mismatch')
receipts['large_exponent_control']={'b':large['b'],'signed_seed_count':large['signed_seed_count'],'primitive_triangle_candidates':large['primitive_triangle_candidates'],'escape_rejections':large['escape_rejections'],'cycles':large['cycles'],'normal_optimized_byte_equal':True,'independent_without_escape_guard_match':True}
receipts['critical_cycle_control_status']=read(root+'/final_addendum_controls.json')['status']
receipts['status']='PASS'
receipts['reviewed_sources_sha256']={p:hashlib.sha256(open(primary+'/'+p,'rb').read()).hexdigest() for p in ('rational_cycle_signature_probe.py','denominator_signature_theorem.md','denominators_through_400.json','denominators_through_1000.json','finite_lift_obstruction.py','rational_cycle_pair_probe.py','critical_cycle_controls.py')}
print(json.dumps(receipts,indent=2))
