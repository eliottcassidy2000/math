#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;using U=uint64_t;
vector<U> dp(const vector<unsigned>&r){int n=r.size(),N=1<<n;vector<U>out(N*n);for(int v=0;v<n;v++)out[(1<<v)*n+v]=1;for(int mask=1;mask<N;mask++)for(int v=0;v<n;v++){U x=out[mask*n+v];if(!x)continue;unsigned todo=r[v]&~mask;while(todo){unsigned b=todo&-todo;todo-=b;int u=__builtin_ctz(b);out[(mask|b)*n+u]+=x;}}return out;}
bool reaches(const vector<unsigned>&r){unsigned seen=1,last=0,full=(1u<<r.size())-1;while(last!=seen){last=seen;for(int i=0;i<(int)r.size();i++)if(seen>>i&1)seen|=r[i];}return seen==full;}

#include <fstream>
#include <map>
#include <array>
struct Request{int root,mode,k;U value;};
int main(int argc,char**argv){if(argc!=2)return 2;ifstream in(argv[1]);map<U,vector<Request>>groups;U packed;U v=69;while(in>>hex>>packed){U code=packed&((U(1)<<36)-1);int root=(packed>>36)&15,n=(packed>>40)&15,op=(packed>>44)&1,mode=(packed>>45)&3,k=(packed>>47)&15,st=(packed>>51)&1;if(!st||n<4||n>9||root>=n||mode>2)return 3;groups[code|(U(n)<<40)|(U(op)<<44)].push_back({root,mode,k,v});v+=2;}if(v!=702191)return 4;U A[11]={1,1,1},C[11]={0,1,0},D[11]={0,1,2};for(int k=3;k<=10;k++){A[k]=A[k-1]+A[k-2]+A[k-3];C[k]=C[k-1]+C[k-2]+C[k-3];D[k]=D[k-1]+D[k-2]+D[k-3];}U count=0;array<U,3>modes{};int maxorder=0;
for(auto&kv:groups){U key=kv.first,code=key&((U(1)<<36)-1);int n=(key>>40)&15,op=(key>>44)&1;vector<unsigned>r(n),rev(n);int bit=n*(n-1)/2-1;for(int j=1;j<n;j++)for(int i=0;i<j;i++,bit--){bool a=(code>>bit&1)^op;if(a)r[i]|=1u<<j;else r[j]|=1u<<i;}for(int i=0;i<n;i++)rev[i]=r[i]^(((1u<<n)-1)^(1u<<i));if(!reaches(r)||!reaches(rev))return 5;auto E=dp(r),S=dp(rev);int full=(1<<n)-1;vector<U>H(1<<n),B(n);H[0]=1;for(int mask=1;mask<=full;mask++)for(int i=0;i<n;i++)H[mask]+=E[mask*n+i];for(int root=0;root<n;root++){int rest=full^(1<<root);for(int a=rest;;a=(a-1)&rest){B[root]+=H[a]*S[(full^a)*n+root];if(!a)break;}}
for(auto q:kv.second){U want=0;int outorder=0;if(q.mode==0){if(q.k!=2)return 6;want=H[full]+2*B[q.root];outorder=n+2;}else if(q.mode==1){if(q.k<1||q.k>8)return 7;want=H[full]+(U(1)<<q.k)*B[q.root];outorder=n+q.k+1;}else{if(q.k<3||q.k>10)return 8;want=H[full]*A[q.k]-E[full*n+q.root]*C[q.k]+B[q.root]*D[q.k];outorder=n+q.k;}if(want!=q.value)return 9;count++;modes[q.mode]++;maxorder=max(maxorder,outorder);}}
cout<<"{\"status\":\"PASS\",\"verified_strong_constructions\":"<<count<<",\"distinct_physical_base_graphs\":"<<groups.size()<<",\"maximum_output_order\":"<<maxorder<<",\"mode_counts\":["<<modes[0]<<','<<modes[1]<<','<<modes[2]<<"]}\n";return 0;}
