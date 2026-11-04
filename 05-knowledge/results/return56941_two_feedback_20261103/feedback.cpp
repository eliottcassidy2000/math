#define main inherited_search_main
#include "../adaptive_four_returns_20261103/lift_search.cpp"
#undef main
T family(int A,int m,int t=11){int n=m+t+1,p=n-1;T Tn(n);
 for(int i=0;i<p;i++){for(int j=i+1;j<p;j++)Tn[i]|=1<<j;bool b=i<m||!(A>>(i-m)&1);if(b)Tn[p]|=1<<i;else Tn[i]|=1<<p;}return Tn;}
int mask(int s,int k,int c){int m=s+k,z=c&((1<<s)-1);if(c>>s&1)z|=((1<<k)-1)<<s;z|=(c>>(s+1))<<m;return z;}
#ifndef FEEDBACK_NO_MAIN
int main(int argc,char**argv){int A=1243,maxs=argc>1?stoi(argv[1]):3;
 for(int s=0;s<=maxs;s++){vector<vector<I>>H;for(int k=1;k<=3;k++){T t=family(A,s+k);H.push_back(ears(t));}
  int count=1<<(bit_width((unsigned)A)+s+2);I target2=3729LL<<s,target1=2*target2-3181,best=INT64_MAX;vector<array<I,4>>close;int zeros=0;
  for(int c=0;c<count;c++){I h1=H[0][mask(s,1,c)],h2=H[1][mask(s,2,c)],h3=H[2][mask(s,3,c)];I n2=h2-h1,n3=h3-3*h2+2*h1;
   if(n2%2||n3%6||n2<0||n3<0)throw runtime_error("run coefficients");I c2=n2/2,c3=n3/6;if(!c3)zeros++;
   if(!c3&&c2==target2){I d=abs(h1-target1);if(d<best){best=d;close.clear();}if(d==best)close.push_back({c,h1,c2,c3});}
   if(!c3&&c2==target2&&h1==target1)cout<<"{\"kind\":\"hit\",\"s\":"<<s<<",\"cut_code\":"<<c<<"}\n";
  }
  cout<<"{\"kind\":\"layer\",\"s\":"<<s<<",\"cuts\":"<<count<<",\"zero_three_run\":"<<zeros<<",\"minimum_constant_error\":"<<best<<",\"closest\":[";
  for(int j=0;j<(int)close.size();j++){if(j)cout<<',';cout<<'['<<close[j][0]<<','<<close[j][1]<<','<<close[j][2]<<','<<close[j][3]<<']';}cout<<"]}\n";cout.flush();
 }
}
#endif
