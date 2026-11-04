#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;using U=uint64_t;
struct Bits{vector<U>b;U count=0;void put(unsigned x){unsigned j=x/64;if(j>=b.size())b.resize(max((size_t)j+1,b.size()*2+1));U z=U(1)<<(x%64);if(!(b[j]&z)){b[j]|=z;count++;}}bool has(unsigned x)const{return x/64<b.size()&&(b[x/64]>>(x%64)&1);}};
vector<U> ends(vector<unsigned>r){int n=r.size(),N=1<<n;vector<U>d(N*n);for(int s=1;s<N;s++)for(int v=0;v<n;v++)if(s>>v&1){int prev=s^(1<<v);if(!prev){d[s*n+v]=1;continue;}U sum=0;for(int q=prev;q;q&=q-1){int u=__builtin_ctz((unsigned)q);if(r[u]>>v&1)sum+=d[prev*n+u];}d[s*n+v]=sum;}return d;}
bool strong(const vector<unsigned>&r){unsigned full=(1u<<r.size())-1;for(int start=0;start<(int)r.size();start++){unsigned seen=1u<<start,last=0;while(seen!=last){last=seen;for(int v=0;v<(int)r.size();v++)if(seen>>v&1)seen|=r[v];}if(seen!=full)return false;}return true;}
int main(int argc,char**argv){if(argc<3)return 2;ifstream input(argv[1]);ofstream rowsout(argv[2]);vector<Bits>sets(7);array<U,10>graphcounts{},rootcounts{},distinct{};array<unordered_set<U>,10>triples;vector<U>witness((702191-69)/2+1);int n,sc;U code;U A[11]={1,1,1},C[11]={0,1,0},D[11]={0,1,2};for(int k=3;k<=10;k++){A[k]=A[k-1]+A[k-2]+A[k-3];C[k]=C[k-1]+C[k-2]+C[k-3];D[k]=D[k-1]+D[k-2]+D[k-3];}
while(input>>n>>code>>sc){vector<unsigned>base(n);int bit=n*(n-1)/2-1;for(int j=1;j<n;j++)for(int i=0;i<j;i++,bit--){if(code>>bit&1)base[i]|=1u<<j;else base[j]|=1u<<i;}for(int op=0;op<(sc?1:2);op++){vector<unsigned>r=base;if(op)for(int i=0;i<n;i++)r[i]^=((1u<<n)-1)^(1u<<i);vector<unsigned>rev(n);for(int i=0;i<n;i++)rev[i]=r[i]^(((1u<<n)-1)^(1u<<i));auto E=ends(r),S=ends(rev);int N=1<<n,full=N-1;vector<U>H(N);H[0]=1;for(int s=1;s<N;s++)for(int i=0;i<n;i++)H[s]+=E[s*n+i];bool st=strong(r);graphcounts[n]++;
 for(int root=0;root<n;root++){rootcounts[n]++;U h=H[full],e=E[full*n+root],B=0;int rest=full^(1<<root);for(int sub=rest;;sub=(sub-1)&rest){B+=H[sub]*S[(full^sub)*n+root];if(!sub)break;}if(B<h||e>h||h>=1<<19||e>=1<<19||B>=1<<19)return 3;triples[n].insert(h|(e<<19)|(B<<38)|(U(st)<<57));
 auto mark=[&](unsigned val,int mode,int k){if(!(val&1))exit(4);sets[5].put(val);if(st)sets[6].put(val);if(val>=69&&val<=702191){unsigned i=(val-69)/2;if(!witness[i]||(st&&!(witness[i]>>51&1)))witness[i]=code|(U(root)<<36)|(U(n)<<40)|(U(op)<<44)|(U(mode)<<45)|(U(k)<<47)|(U(st)<<51);} };
 unsigned two=h+2*B;sets[0].put(two);if(st)sets[1].put(two);mark(two,0,2);sets[4].put(two);
 for(int k=1;k<=8;k++){unsigned v=h+(U(1)<<k)*B;sets[2].put(v);mark(v,1,k);}
 for(int k=3;k<=10;k++){unsigned v=h*A[k]-e*C[k]+B*D[k];sets[3].put(v);sets[4].put(v);mark(v,2,k);}
 sets[4].put(h-e+B);
 }
}}
U roots=0,distinctsum=0;cout<<"{\"orders\":[";for(int k=4;k<=9;k++){if(k>4)cout<<',';cout<<"{\"n\":"<<k<<",\"graphs\":"<<graphcounts[k]<<",\"roots\":"<<rootcounts[k]<<",\"distinct_triples\":"<<triples[k].size()<<'}';roots+=rootcounts[k];distinctsum+=triples[k].size();}cout<<"],\"rooted_total\":"<<roots<<",\"distinct_total\":"<<distinctsum<<",\"families\":[";
const char*names[]={"two","two_strong","nearjoin","tower3to10","tower1to10","source_union","source_union_strong"};for(int k=0;k<7;k++){if(k)cout<<',';unsigned v=69;if(!sets[k].has(69)){cerr<<"69 absent "<<k<<'\n';return 5;}while(sets[k].has(v+2))v+=2;cout<<"{\"family\":\""<<names[k]<<"\",\"maximum_complete_odd_interval\":"<<v<<",\"next_missing\":"<<v+2<<",\"distinct_values\":"<<sets[k].count<<",\"missed_below69\":[";bool comma=false;for(int x=1;x<69;x+=2)if(!sets[k].has(x)){if(comma)cout<<',';cout<<x;comma=true;}cout<<"]}";}cout<<"],\"height1_only_additions\":[";bool cc=false;for(unsigned x=1;x<sets[4].b.size()*64;x+=2)if(sets[4].has(x)&&!sets[0].has(x)&&!sets[3].has(x)){if(cc)cout<<',';cout<<x;cc=true;}cout<<"]}\n";
for(unsigned i=0;i<witness.size();i++){if(i+1==witness.size()){if(witness[i])return 6;break;}if(!witness[i]||!(witness[i]>>51&1))return 7;rowsout<<hex<<witness[i]<<'\n';}return 0;}
