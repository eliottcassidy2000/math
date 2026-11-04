#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using I=__int128_t;
string dec(I v){if(v==0)return "0";bool neg=v<0;if(neg)v=-v;string s;while(v){s+=char('0'+v%10);v/=10;}if(neg)s+='-';reverse(s.begin(),s.end());return s;}
int main(){int count;cin>>count;cout<<"[";for(int g=0;g<count;g++){int k,z,target;cin>>k>>z>>target;vector<int>r(k),cap(k),stride(k);for(int&x:r)cin>>x;for(int&x:cap)cin>>x;int N=1,n=0;for(int i=0;i<k;i++){stride[i]=N;N*=cap[i]+1;n+=cap[i];}vector<I>H((size_t)N*k*n),D(H.size());for(int i=0;i<k;i++)H[((size_t)stride[i]*k+i)*n]=1;vector<I>sumH(n),sumD(n);for(int state=1;state<N;state++){vector<int>used(k);int len=0,ones=0,twos=0;bool complete=true;for(int i=0;i<k;i++){used[i]=(state/stride[i])%(cap[i]+1);len+=used[i];if(!used[i])complete=false;if(cap[i]==2){ones+=used[i]==1;twos+=used[i]==2;}}
for(int last=0;last<k;last++)if(used[last]){I*h=&H[((size_t)state*k+last)*n],*d=&D[((size_t)state*k+last)*n];if(complete){vector<I>w(ones+1);w[0]=I(1)<<twos;for(int a=1;a<=ones;a++)w[a]=w[a-1]*2*(ones-a+1)/a;for(int a=0;a<len;a++)for(int b=0;b<=ones;b++){sumH[a+b]+=h[a]*w[b];sumD[a+b]+=d[a]*w[b];}}
for(int v=0;v<k;v++)if(v!=last&&used[v]<cap[v]){int next=state+stride[v],edge=(r[last]>>v)&1,sgn=(last==z&&v==target)-(last==target&&v==z);I*hh=&H[((size_t)next*k+v)*n],*dd=&D[((size_t)next*k+v)*n];for(int a=0;a<len;a++){hh[a]+=edge*h[a];hh[a+1]+=h[a];dd[a]+=edge*d[a]+sgn*h[a];dd[a+1]+=d[a];}}}}
if(g)cout<<",";cout<<"{\"id\":"<<g<<",\"target\":"<<target<<",\"H\":[";for(int a=0;a<n;a++){if(a)cout<<",";cout<<dec(sumH[a]);}cout<<"],\"skew\":[";for(int a=0;a<n;a++){if(a)cout<<",";cout<<dec(sumD[a]);}cout<<"]}";cerr<<"case "<<g<<" done\n";}cout<<"]\n";}
