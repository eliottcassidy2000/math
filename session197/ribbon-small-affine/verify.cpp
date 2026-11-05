// Exact modular certificates and exhaustive low-output word enumeration.
#include <array>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
using namespace std;
constexpr int N=62,K=8;using Mat=array<int,N*N>;
vector<pair<int,int>> sp[K][N];vector<array<int,K>> dfa;
constexpr array<int,5> primes={3,5,7,11,13};
long long short_classes[7]={},short_pass[7][5]={};
long long low_words[12]={},active_words[12]={},long_classes[12]={},long_pass[12][5]={};
set<string> candidates; int multiplier=4;
int mod(int x,int p){x%=p;return x<0?x+p:x;}
Mat identity(){Mat A{};for(int i=0;i<N;i++)A[i*N+i]=1;return A;}
Mat multiply(const Mat&A,int a,int p){Mat C{};for(int i=0;i<N;i++)for(int k=0;k<N;k++)if(A[i*N+k])for(auto [j,c]:sp[a][k])C[i*N+j]+=A[i*N+k]*c;for(int&x:C)x=mod(x,p);return C;}
Mat product(const string&w,int p){Mat A=identity();for(char a:w)A=multiply(A,a-'0',p);return A;}
int determinant_minus(Mat A,int target,int p){
 for(int i=0;i<N;i++)A[i*N+i]=mod(A[i*N+i]-target,p);int d=1;
 for(int k=0;k<N;k++){
  int r=k;while(r<N&&!A[r*N+k])r++;if(r==N)return 0;
  if(r!=k){for(int j=k;j<N;j++)swap(A[k*N+j],A[r*N+j]);d=mod(-d,p);}
  int x=A[k*N+k],iv=1;while(iv*x%p!=1)iv++;d=d*x%p;
  for(int i=k+1;i<N;i++)if(A[i*N+k]){int z=A[i*N+k]*iv%p;A[i*N+k]=0;for(int j=k+1;j<N;j++)A[i*N+j]=mod(A[i*N+j]-z*A[k*N+j],p);}
 }
 return d;
}
string rotation(const string&w,int k){return w.substr(k)+w.substr(0,k);}
bool canonical(const string&w){for(int k=1;k<(int)w.size();k++)if(rotation(w,k)<w)return false;return true;}
int certify(const string&w,const Mat*mod3=nullptr){for(int z=0;z<5;z++){Mat A=(z==0&&mod3)?*mod3:product(w,primes[z]);if(determinant_minus(A,1,primes[z])||determinant_minus(A,multiplier,primes[z]))return z;}cerr<<"UNRESOLVED "<<w<<"\n";exit(1);}
void short_dfs(string&w,const Mat&A){
 if(!w.empty()&&canonical(w)){int l=w.size();short_classes[l]++;short_pass[l][certify(w,&A)]++;}
 if(w.size()==6)return;
 for(int a=0;a<8;a++){w.push_back('0'+a);auto C=multiply(A,a,3);short_dfs(w,C);w.pop_back();}
}
bool cyclically_low(const string&w,string&best){best=w;for(int k=0;k<(int)w.size();k++){string z=rotation(w,k);best=min(best,z);int state=0;for(char a:z){state=dfa[state][a-'0'];if(state<0)return false;}}return true;}
void low_dfs(string&w,int state,array<unsigned,11> reach){
 int n=w.size();if(n){low_words[n]++;if(n>=7){bool active=false;for(int i=n-3;i<n;i++)if(reach[i]&7u)active=true;if(active){active_words[n]++;string best;if(cyclically_low(w,best))candidates.insert(best);}}}
 if(n==11)return;
 for(int a=0;a<8;a++){
  int next=dfa[state][a];if(next<0)continue;auto rr=reach;unsigned out=0;
  for(int i=max(0,n-3);i<n;i++)if(a>>(i-(n-3))&1)out|=1u<<i;
  unsigned in=((1u<<n)-1)^out,newreach=1u<<n;
  for(int i=0;i<n;i++)if(out>>i&1)newreach|=reach[i];
  for(int i=0;i<n;i++)if(reach[i]&in)rr[i]|=newreach;
  rr[n]=newreach;w.push_back('0'+a);low_dfs(w,next,rr);w.pop_back();
 }
}
int main(int argc,char**argv){
 if(argc!=2)return 4;multiplier=stoi(argv[1]);if(multiplier!=3&&multiplier!=4)return 5;
 int n,k;cin>>n>>k;if(n!=N||k!=K)return 2;
 for(int a=0;a<K;a++)for(int i=0;i<N;i++)for(int j=0;j<N;j++){int x;cin>>x;if(x)sp[a][i].push_back({j,x});}
 int states;cin>>states;dfa.resize(states);for(auto&row:dfa)for(int&x:row)cin>>x;if(!cin)return 3;
 string word;short_dfs(word,identity());array<unsigned,11>reach{};low_dfs(word,0,reach);
 for(const string&w:candidates){int l=w.size();long_classes[l]++;long_pass[l][certify(w)]++;}
 cout<<"{\n  \"multiplier\": "<<multiplier<<",\n  \"dfa_states\": "<<states<<",\n  \"short\": [\n";
 for(int l=1;l<=6;l++){cout<<"    {\"length\": "<<l<<", \"classes\": "<<short_classes[l]<<", \"first_success_by_prime\": [";for(int z=0;z<5;z++)cout<<(z?", ":"")<<short_pass[l][z];cout<<"]}"<<(l<6?",":"")<<"\n";}
 cout<<"  ],\n  \"long\": [\n";
 for(int l=7;l<=11;l++){cout<<"    {\"length\": "<<l<<", \"low_words\": "<<low_words[l]<<", \"active_words\": "<<active_words[l]<<", \"filtered_classes\": "<<long_classes[l]<<", \"first_success_by_prime\": [";for(int z=0;z<5;z++)cout<<(z?", ":"")<<long_pass[l][z];cout<<"]}"<<(l<11?",":"")<<"\n";}
 cout<<"  ],\n  \"unresolved\": 0\n}\n";
}
