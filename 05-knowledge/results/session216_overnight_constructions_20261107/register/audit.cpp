#define main prior_main
#include "three_core_audit.cpp"
#undef main
int main(){
 vector<S> seeds;for(int cut=1;cut<16;cut++)seeds.push_back(literal(cores[0],cut));
 const int LIMIT=100000;set<pair<I,I>>joint;vector<bool>hit(LIMIT+1);for(S z:seeds)joint.insert({z.a,z.e});
 for(int a=5;a<=LIMIT;a++){
  auto it=joint.lower_bound({a,-1});bool actual=it!=joint.end()&&it->first==a; // bounded support control
  hit[a]=actual;if(actual){if(2*a<=LIMIT)joint.insert({2*a,0});if(2*a+5<=LIMIT)joint.insert({2*a+5,5});}
  if(actual!=(a!=8))return 21;
 }
 for(int a=5;a<=LIMIT;a++)for(int e=0;e<=5;e++){bool expected=false;
 if(e==0)expected=(a>=10&&a%2==0&&a!=16);
 if(e==5)expected=(a>=15&&a%2==1&&a!=21);
 if(e==1)expected=set<I>{5,7,9}.count(a);
 if(e==2)expected=set<I>{6,12,14,16}.count(a);
 if(e==3)expected=set<I>{11,13,15,21}.count(a);
 if(e==4)expected=set<I>{18,20,22}.count(a);
 if(bool(joint.count({a,e}))!=expected)return 22;
 }
 S a=append(literal(cores[0],6),1),b=append(append(literal(cores[0],1),1),1);
 if(!(a.h==5&&a.a==35&&a.t==24&&a.e==5&&b.a==35&&b.t==26&&b.e==5))return 23;
 if(paths(close(padded(0,0,1),fullcut(0,6,0,1,0,1)))!=51||paths(close(padded(0,0,2),fullcut(0,1,0,2,0,3)))!=49)return 28;
 auto LA=padded(0,6,1),LB=padded(0,3,4);int ca=fullcut(0,10,6,1,59,1),cb=fullcut(0,10,3,4,7,3);a=literal(LA,ca);b=literal(LB,cb);
 if(!(a==S{5,4096,3807,3219,2984,5,5})||!(b==S{5,4096,3807,3219,2986,5,5}))return 24;
 auto ga=close(LA,ca),gb=close(LB,cb);if(paths(ga)!=4635||paths(gb)!=4633||!strong(ga)||!strong(gb))return 25;
 auto lA=padded(0,3,1),lB=padded(0,0,4);int cA=fullcut(0,10,3,1,3,1),cB=fullcut(0,10,0,4,0,3);auto x=literal(lA,cA),y=literal(lB,cB);
 if(!(x==S{5,512,223,398,168,0,5})||!(y==S{5,512,223,398,168,2,5}))return 26;
 int controls=0;for(int left=0;left<=3;left++)for(int right=0;left+right<=3;right++)for(int p=0;p<(1<<left);p++)for(int r=0;r<(1<<right);r++){
 S A=a,B=b;for(int i=0;i<left;i++){A=prepend(A,p>>i&1);B=prepend(B,p>>i&1);}for(int i=right-1;i>=0;i--){A=append(A,r>>i&1);B=append(B,r>>i&1);}
 S Bfixed=B;Bfixed.t-=2;if(!(A==Bfixed))return 27;controls++;
 }
 cout<<"{\"register_limit\":"<<LIMIT<<",\"register_support_pass\":true,\"joint_endpoint_support_pass\":true,\"two_sided_pair_literal_counts\":[4635,4633],\"local_endpoint_values\":[0,2],\"two_sided_continuation_controls\":"<<controls<<"}"<<endl;
}
