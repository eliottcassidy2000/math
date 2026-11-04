#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;using Size=array<int,5>;
struct Form{vector<pair<int,int>> terms;long long den;int weight;};
int code(Size a){int c=0;for(int x:a)c=6*c+x;return c;}
const int nums[6][5]={{0},{1},{-1,1},{2,-3,1},{-6,11,-6,1},{24,-50,35,-10,1}};
const int fac[]={1,1,2,6,24,120};
array<uint64_t,3> count(const array<unsigned,5>&Q,Size sizes,unsigned a,unsigned c){int off[5],n=0;for(int i=0;i<5;i++){off[i]=n;n+=sizes[i];}int m=n+2;vector<unsigned>R(m);R[n]=1u<<(n+1);for(int i=0;i<5;i++)for(int u=0;u<sizes[i];u++){int v=off[i]+u;for(int j=u+1;j<sizes[i];j++)R[v]|=1u<<(off[i]+j);for(int j=0;j<5;j++)if(Q[i]>>j&1)R[v]|=((1u<<sizes[j])-1)<<off[j];if(a>>i&1)R[v]|=1u<<n;else R[n]|=1u<<v;if(c>>i&1)R[v]|=1u<<(n+1);else R[n+1]|=1u<<v;}unsigned N=1u<<m;vector<uint32_t>d(size_t(N)*m);for(int i=0;i<m;i++)d[size_t(1u<<i)*m+i]=1;for(unsigned S=1;S<N;S++){unsigned bs=S;while(bs){unsigned b=bs&-bs;bs-=b;int i=__builtin_ctz(b);uint32_t z=d[size_t(S)*m+i];if(!z)continue;unsigned next=R[i]&~S;while(next){unsigned t=next&-next;next-=t;d[size_t(S|t)*m+__builtin_ctz(t)]+=z;}}}vector<uint64_t>H(N),two(N);for(unsigned S=1;S<N;S++)for(int i=0;i<m;i++)H[S]+=d[size_t(S)*m+i];for(unsigned S=1;S<N;S++)for(unsigned A=(S-1)&S;A;A=(A-1)&S)two[S]+=H[A]*H[S^A];uint64_t three=0;unsigned S=N-1;for(unsigned A=(S-1)&S;A;A=(A-1)&S)three+=H[A]*two[S^A];return {H[S],two[S],three};}
void list_sizes(int i,int rest,Size a,vector<Size>&out){if(i==5){out.push_back(a);return;}for(int t=0;t<=rest;t++){a[i]=1+t;list_sizes(i+1,rest-t,a,out);}}
void terms(int i,Size k,Size a,int coef,const map<int,int>&idx,vector<pair<int,int>>&out){if(i==5){out.push_back({idx.at(code(a)),coef});return;}for(int j=1;j<=k[i];j++){a[i]=j;terms(i+1,k,a,coef*nums[k[i]][j-1],idx,out);}}
Form make_form(Size k,const map<int,int>&idx){Form f;f.den=1;f.weight=1;int q[]={0,1,14,36,24};for(int x:k){f.den*=fac[x];if(x<5)f.weight*=q[x];else f.weight=0;}Size a{};terms(0,k,a,1,idx,f.terms);return f;}
long long eval(const Form&f,const vector<uint64_t>&E){long long n=0;for(auto [i,c]:f.terms)n+=static_cast<long long>(c)*E[i];if(n%f.den||n<0)throw runtime_error("exact triangular interpolation");return n/f.den;}
int main(){Size one{1,1,1,1,1};vector<Size>sizes{one};for(int i=0;i<5;i++)for(int k:{2,3}){Size s=one;s[i]=k;sizes.push_back(s);}for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)for(int k:{2,3})for(int l:{2,3}){Size s=one;s[i]=k;s[j]=l;sizes.push_back(s);}map<int,int>idx;for(int i=0;i<(int)sizes.size();i++)idx[code(sizes[i])]=i;if(idx.size()!=51)throw runtime_error("51 models");vector<Form>forms;for(auto s:sizes)forms.push_back(make_form(s,idx));int cases;cin>>cases;long long calls=0;while(cases--){array<unsigned,5>Q;for(auto&v:Q)cin>>v;for(unsigned a=0;a<32;a++)if(__builtin_popcount(a)<=2)for(unsigned c=0;c<32;c++)if(__builtin_popcount(c)<=2)for(int rev=0;rev<2;rev++){auto R=Q;if(rev)for(int i=0;i<5;i++)R[i]=31u^(1u<<i)^Q[i];array<vector<uint64_t>,3>E;for(auto s:sizes){auto v=count(R,s,a,c);for(int h=0;h<3;h++)E[h].push_back(v[h]);calls++;}for(auto f:forms)for(int h=0;h<3;h++)cout<<(eval(f,E[h])%8)<<' ';cout<<'\n';}}cerr<<"physical_profiles "<<calls<<" max_names 11\n";return cin.fail()?1:0;}
