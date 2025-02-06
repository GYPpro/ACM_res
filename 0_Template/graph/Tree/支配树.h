#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define pb push_back
#define emp empty
#define fi first
#define se second
const int N=3e5+7;
const ll inf=1e16+7;
const ll mod=998244353;
signed main()
{
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    void solve();
	int t=1;
//    cin>>t;
	while(t--)solve();
    return 0;
}
struct node
{
	int v;
	int nex;
}a[N<<2];
int n,m,u,v,cnt=0,dfc=0;  
int siz[N],dfn[N],pos[N],sdm[N],idm[N],fa[N],fth[N],mn[N];
//siz儿子数、dfn深搜序、pos深搜序对应结点，sdm半支配，idm最近支配，fa并查集，fth深搜树前驱，mn半支配点深搜序最小的祖宗
int h[3][N<<1];//三幅图的表头，分别是原图、反图、支配树图
void clear()
{
	memset(h,0,sizeof(h));
	memset(dfn,0,sizeof(dfn));
	memset(siz,0,sizeof(siz));
	cnt=0;
	dfc=0;
}
void add(int u,int v,int x)
{
	a[++cnt]=(node){v,h[x][u]};
	h[x][u]=cnt;
}
void dfs(int u)
{
	dfn[u]=++dfc;
	pos[dfc]=u;
	for(int i=h[0][u];i;i=a[i].nex)
	{
		int v=a[i].v;
		if(dfn[v])continue;
		dfs(v);
		fth[v]=u;
	}
}
int find(int u)//这是一个带权并查集？
{
	if(fa[u]==u)return u;
	int temp=fa[u];
	fa[u]=find(fa[u]);
	if(dfn[sdm[mn[temp]]]<dfn[sdm[mn[u]]])
	{
		mn[u]=mn[temp];
	}
	return fa[u];
}
void tar(int st)
{
	dfs(st);
	for(int i=1;i<=n;i++)
	{
		sdm[i]=fa[i]=mn[i]=i;
	}
	for(int i=dfc;i>=2;i--)
	{
		int u=pos[i];
		int res=mod;
		for(int j=h[1][u];j;j=a[j].nex)
		{
			int v=a[j].v;
			if(!dfn[v])continue;
			find(v);
			if(dfn[v]<dfn[u])res=min(res,dfn[v]);
			else res=min(res,dfn[sdm[mn[v]]]);
		}
		sdm[u]=pos[res];
		fa[u]=fth[u];
		add(sdm[u],u,2);
		u=fth[u];
		for(int j=h[2][u];j;j=a[j].nex)
		{
			int v=a[j].v;
			find(v);
			if(u==sdm[mn[v]])
			{
				idm[v]=u;
			}
			else 
			{
				idm[v]=mn[v];
			}
		}
		h[2][u]=0;
	}
	for(int i=2;i<=dfc;i++)
	{
		int u=pos[i];
		if(idm[u]!=sdm[u])idm[u]=idm[idm[u]];
	}
	for(int i=dfc;i>=2;i--)
	{
		int u=pos[i];
		++siz[u];
		siz[idm[u]]+=siz[u];
	}
	++siz[st];
}//这里支配树消失了，但可以通过idm重建
void solve()
{
	cin>>n>>m;
	clear();
	for(int i=1;i<=m;i++)
	{
		int u,v;
		cin>>u>>v;
		add(u,v,0);
		add(v,u,1);
	}
	tar(1);
	for(int i=1;i<=n;i++)cout<<siz[i]<<" ";
	cout<<"\n";
}
