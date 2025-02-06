#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define inf 0x3f3f3f3f
#define mod 998244353
int ksm(int a,int b)
{
    int ans=1;
    while(b)
    {
        if(b&1)ans=(ll)ans*a%mod;
        a=(ll)a*a%mod;
        b>>=1;
    }
    return ans;
}
//多项式开根号，需要NTT+INV+二次剩余(cipolla)
const int gen=3,gi=ksm(gen,mod-2),N=(1<<22)|10,ny2=ksm(2,mod-2);
int main()
{
    void solve();
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    solve();
    return 0;
}
struct comp{int r,i;};
comp image_mul(comp a,comp b,int w)
{
    comp ans;
    ans.r=((ll)a.r*b.r%mod+(ll)a.i*b.i%mod*w%mod)%mod;
    ans.i=((ll)a.r*b.i%mod+(ll)a.i*b.r%mod)%mod;
    return ans;
}
comp ksm_image(comp a,int b,int w)
{
    comp res;
    res.r=1;res.i=0;
    while(b)
    {
        if(b&1)res=image_mul(res,a,w);
        a=image_mul(a,a,w);
        b>>=1;
    }
    return res;
}
int cipolla(int n)  //求二次剩余，仅用于求常数mod意义下的二次方根
{
    int w,a=rand()%mod;
    w=((ll)a*a-n+mod)%mod;
    if(ksm(n,(mod-1)>>1)==0)return 0;
    while(ksm(w,(mod-1)>>1)!=mod-1)a=rand()%mod,w=((ll)a*a-n+mod)%mod;
    comp ans;
    ans.r=a;
    ans.i=1;
    ans=ksm_image(ans,((ll)mod+1)>>1,w);
    int a1=ans.r,a2=(mod-a1)%mod;
    if(a1>a2)swap(a1,a2);
    return a1;
}

int a[N],b[N],c[N],d[N],e[N],rev[N]={0};
void chang(int *A,int n)  
{
    for(int i=0;i<n;i++)rev[i]=(rev[i>>1]>>1)|((i&1)?(n>>1):0);
    for(int i=0;i<n;i++)if(i<rev[i])swap(A[i],A[rev[i]]);
}
void NTT(int *A,int n,int f)  //多项式相乘
{
    chang(A,n);
    for(int mid=1;mid<n;mid<<=1)
    {
        int g1=ksm((f==1)?gen:gi,(mod-1)/(mid<<1));
        for(int j=0;j<n;j+=(mid<<1))
        {
            int gk=1;
            for(int k=0;k<mid;k++,gk=(ll)gk*g1%mod)
            {
                int x=A[j+k],y=(ll)A[j+k+mid]*gk%mod;
                A[j+k]=((ll)x+(ll)y)%mod;
                A[j+k+mid]=((ll)x-(ll)y+mod)%mod;
            }
        }
    }
    if(f==1)return ;
    int inv=ksm(n,mod-2);
    for(int i=0;i<n;i++)A[i]=(ll)A[i]*inv%mod;
}
void INV(int n,int *A,int *B)//多项式求逆
{
    stack<int> stk;
    while(n!=1){stk.push(n);n=(n+1)>>1;}
    stk.push(1);
    int mx=1;
    while(!stk.empty())
    {
        n=stk.top();
        stk.pop();
        if(n==1){B[0]=ksm(A[0],mod-2);continue;}
        while(mx<(n<<1))mx<<=1;
        for(int i=0;i<n;i++)c[i]=A[i];
        for(int i=n;i<mx;i++)c[i]=0;
        NTT(c,mx,1);NTT(B,mx,1);
        for(int i=0;i<mx;i++)B[i]=(2ll-(ll)c[i]*B[i]%mod+mod)%mod*B[i]%mod;
        NTT(B,mx,-1);
        for(int i=n;i<mx;i++)B[i]=0;
    }
}
void SQRT(int n,int *A,int *B)  //多项式开根
{
    //思想和INV相同：分治、套公式
    stack<int> stk;
    while(n!=1){stk.push(n);n=(n+1)>>1;}
    stk.push(1);
    int mx=1;
    while(!stk.empty())
    {
        n=stk.top();
        stk.pop();
        if(n==1){B[0]=cipolla(A[0]);continue;}  //常数求二次剩余
        while(mx<(n<<1))mx<<=1;
        for(int i=0;i<n;i++)d[i]=A[i],e[i]=0;
        for(int i=n;i<mx;i++)d[i]=0,e[i]=0;
        //e是B的逆元

        INV(n,B,e);NTT(d,mx,1);NTT(B,mx,1);NTT(e,mx,1); 

        for(int i=0;i<mx;i++)B[i]=((ll)d[i]*e[i]%mod+(ll)B[i])%mod*ny2%mod;//套公式
        
        NTT(B,mx,-1);
        
        for(int i=n;i<mx;i++)B[i]=0;
    }
}
void solve()
{
    int n;
    cin>>n;
    for(int i=0;i<n;i++)cin>>a[i];
    SQRT(n,a,b);
    for(int i=0;i<n;i++)cout<<b[i]<<" ";
    cout<<"\n";
}