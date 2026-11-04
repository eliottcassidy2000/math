#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
struct State{int raw,avail;int next[11];unsigned char mult[11];bool done;};
vector<State>states;vector<int>indexer;int MOD=8;
void init(int bits){MOD=1<<bits;indexer.assign(1<<20,-1);for(int raw=0;raw<(1<<20);raw++){int repeats=0;bool done=(raw>>18)==3;for(int i=0;i<9;i++){int d=(raw>>(2*i))&3;repeats+=d>=2;if(!d)done=false;}if(repeats>=bits)continue;State s{};s.raw=raw;s.done=done;for(int i=0;i<11;i++){int d=i<9?(raw>>(2*i))&3:(raw>>(i+9))&1,cap=i<9?3:1;if(d<cap){s.avail|=1<<i;s.next[i]=raw+(1<<(i<9?2*i:i+9));s.mult[i]=d==1&&i<9?2:1;}}indexer[raw]=states.size();states.push_back(s);}for(auto&s:states)for(int i=0;i<11;i++)if(s.avail>>i&1){s.next[i]=indexer[s.next[i]];if(s.next[i]<0)s.avail^=1<<i;}}
vector<int> count(const vector<int>&r){vector<unsigned char>d(states.size()*11);for(int i=0;i<11;i++)d[(size_t)indexer[1<<(i<9?2*i:i+9)]*11+i]=1;vector<int> ans(46);for(size_t st=1;st<states.size();st++){auto&s=states[st];unsigned sum=0;for(int last=0;last<11;last++){unsigned v=d[st*11+last];if(!v)continue;if(s.done)sum=(sum+v)&7;int todo=r[last]&s.avail;while(todo){int j=__builtin_ctz((unsigned)todo);todo&=todo-1;auto&cell=d[(size_t)s.next[j]*11+j];cell=(cell+v*s.mult[j])&7;}}if(s.done&&sum){ans[0]=(ans[0]+sum)&7;int t[9];for(int j=0;j<9;j++){int count=(s.raw>>(2*j))&3;t[j]=count==2?3:count==3?2:1;ans[1+j]=(ans[1+j]+sum*(count>=2?3:1))&7;}int idx=10;for(int j=0;j<9;j++)for(int k=j+1;k<9;k++){ans[idx]=(ans[idx]+sum*t[j]*t[k])&7;idx++;}}}return ans;}
vector<int>aug(const vector<int>&q,int a,int c,bool op){vector<int>r(11);for(int i=0;i<9;i++)r[i]=(op?(511^q[i]^(1<<i)):q[i])|(((a>>i)&1)<<9)|(((c>>i)&1)<<10);r[9]=(511^a)|1024;r[10]=511^c;return r;}
int main(){init(3);int cases;cin>>cases;cout<<"node\ta\tc";for(int i=0;i<46;i++)cout<<"\tS"<<i<<"\tT"<<i;cout<<"\n";for(int h=0;h<cases;h++){int id,a,c;cin>>id>>a>>c;vector<int>Q(9);for(int&x:Q)cin>>x;auto S=count(aug(Q,a,c,false)),T=count(aug(Q,a,c,true));cout<<id<<"\t"<<a<<"\t"<<c;for(int j=0;j<46;j++)cout<<"\t"<<S[j]<<"\t"<<T[j];cout<<"\n";if(h%100==0)cerr<<"completed "<<h+1<<" / "<<cases<<"\n";}}
