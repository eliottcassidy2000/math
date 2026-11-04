#define FEEDBACK_NO_MAIN
#include "feedback.cpp"
int main(int argc,char**argv){int limit=argc>1?stoi(argv[1]):2;int original=1243;for(int flip=0;flip<(1<<11);flip++){if(popcount((unsigned)flip)>limit)continue;int A=original^flip;vector<vector<I>>H;
 for(int k=1;k<=3;k++)H.push_back(ears(family(A,k)));
 I target2=3729;int count=1<<13;for(int c=0;c<count;c++){
  I h1=H[0][mask(0,1,c)],h2=H[1][mask(0,2,c)],h3=H[2][mask(0,3,c)];
  I c2=(h2-h1)/2,c3=(h3-3*h2+2*h1)/6;
  for(int offset:{2,4})if(!c3&&c2==target2&&h1==2*target2-3183+offset)
   cout<<"{\"flip\":"<<flip<<",\"edited_A\":"<<A<<",\"cut_code\":"<<c<<",\"offset\":"<<offset<<",\"runs\":["<<h1<<','<<c2<<','<<c3<<"]}\n";
 }cerr<<"done flip "<<flip<<'\n';cout.flush();}}
