#include<bits/stdc++.h>
using namespace std;
#define	ll long long
#define N 100005
const ll inf=1e17+7;
signed main()
{
	ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
	void solve();
	int t=1;
    while(t--)solve();
    return 0;
}
vector<int> a;
int fa[N];
int find(int u)
{
    return (fa[u]==u)?u:fa[u]=find(fa[u]);
}
struct left_tree
{
    int rs,ls,fa;
    ll val;
    int d=-1;
}tree[N];
void init(int x,ll v)
{
    tree[x].val=v;
    tree[x].d=tree[0].d+1;
    tree[x].fa=tree[x].ls=tree[x].rs=0;
}
int merge(int x,int y)
{
    if(!x||!y)return x|y;//存在空堆
    if((tree[x].val==tree[y].val)?x>y:tree[x].val>tree[y].val)swap(x,y);  //维护小根堆，选择权小的作为根x
    tree[x].rs=merge(tree[x].rs,y);  //合并根的右儿子和y
    if(tree[tree[x].rs].d>tree[tree[x].ls].d)swap(tree[x].rs,tree[x].ls);//维护左偏性
    tree[x].d=tree[tree[x].rs].d+1;  //距离增加
    tree[tree[x].rs].fa=x;  //更新父子关系
    return x;
}
void pushup(int x)
{
    if(!x)return;  //无可更新
    if(tree[tree[x].rs].d>tree[tree[x].ls].d)swap(tree[x].rs,tree[x].ls);//维护左偏性
    if(tree[x].d!=tree[tree[x].rs].d+1)  //维护距离
    {
        tree[x].d=tree[tree[x].rs].d+1;
        pushup(tree[x].fa);  //向上更新
    }
}
void erase(int x)
{
    int y=merge(tree[x].rs,tree[x].ls);  //合并x结点的左右儿子
    tree[y].fa=tree[x].fa;  //更新新子树的根的父亲

    if(tree[tree[x].fa].ls==x)tree[tree[x].fa].ls=y;   //维护父子关系
    else if(tree[tree[x].fa].rs==x)tree[tree[x].fa].rs=y;
    
    pushup(tree[y].fa);  //向上维护左偏性
}
void uni(int x,int y)//合并x和y所在的堆
{
    if(a[x]==-1||a[y]==-1)return ;
    int fax=find(x),fay=find(y);
    if(fax==fay)return;
    fa[fax]=fa[fay]=merge(fax,fay);
}
int pop(int x)
{
    if(a[x]==-1)return -1;
    int fa1=find(x);
    int ans=a[fa1];
    a[fa1]=-1;
    fa[tree[fa1].ls]=fa[tree[fa1].rs]=fa[fa1]=merge(tree[fa1].ls,tree[fa1].rs);
    tree[fa[tree[fa1].ls]].fa=0;
    return ans;
}
void solve()
{
    int n,m;
    cin>>n>>m;
    a.resize(n+10);
    tree[0].d=-1;
    for(int i=1;i<=n;i++)
    {
        cin>>a[i];
        init(i,a[i]);
    }
    for(int i=1;i<=n;i++)
    {
        fa[i]=i;
    }
    while(m--)
    {
        int op;
        cin>>op;
        if(op==1)
        {
            int x,y;
            cin>>x>>y;
            uni(x,y);
        }
        else if(op==2)
        {
            int x;
            cin>>x;
            int ans=pop(x);
            cout<<ans<<"\n";
        }
    }
}