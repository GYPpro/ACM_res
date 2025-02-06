//题目：有一张简单无向图，问加至少加多少条新边，才能使得图上既有奇环又有偶环，或不可能

//涉及：树的直径、二分图染色
//     貌似只能判定、不能求数量

struct node  //存边
{
	int u;
	int v;
	int id;
};
struct D  //用于求树的直径的结构体
{
	int d;   //直径   d = maxn + submaxn
	int maxn=0;    //最大子树结点数
	int submaxn=0;  //次大子树结点数
}d[N];
vector<node> edge;  //邻接表

int odd=0,even=0;  //是否存在奇/偶环

vector<int> e[N];  //存边
int col[N]={0};  //二分图染色法
int vis[N]={0};  //
int dfn[N]={0},cnt=0;  //深搜序

int fa[N]={0};   //并查集
int siz[N]={0};  //存连通块的直径

int pre[N]={0},cov[N]={0};  //前驱点、染色体
int pre_edge[N]={0};   //前驱边

int find(int u) //并查集
{
	if(fa[u]==u)return u;
	else return fa[u]=find(fa[u]);
}
void uni(int u,int v)
{
	int fa1=find(u);
	int fa2=find(v);
	fa[fa1]=min(fa1,fa2);
	fa[fa2]=min(fa1,fa2);
}

void dfs(int u,int fa)
{
	if(fa!=0)uni(u,fa);  //深搜时顺便建立父子关系

	dfn[u]=++cnt;   //记录深搜序
	col[u]=col[fa]^1;   //染色
	vis[u]++;     //搜过的点

	d[u].maxn=0;  //系列初始化
	d[u].submaxn=0;
	d[u].d=0;

	for(auto re:e[u])
	{
		int v=u^edge[re].u^edge[re].v;   //取儿子结点

		if(v==fa)continue;  //跳过父亲

		if(col[v]<0)   //儿子未染色
		{
			pre[re]=u;      //记录 边 的前驱点
			pre_edge[v]=re;  //记录 点 的前驱边
			dfs(v,u);

			if(d[u].maxn<d[v].maxn+1)  //回溯时更新最大和次大子树，用于求直径
			{
				d[u].submaxn=d[u].maxn;
				d[u].maxn=d[v].maxn+1;
			}
			else if(d[u].maxn==d[v].maxn+1)
			{
				d[u].submaxn=d[u].maxn;
			}
			else if(d[u].submaxn<d[v].maxn+1)
			{
				d[u].submaxn=d[v].maxn+1;
			}

		}
		else if(col[v]==col[u]^1) //找到偶环
		{
		    even++;//记录
		}
		else if(col[v]==col[u])//找到奇环
		{
			odd++;//记录
			if(dfn[v]>dfn[u])continue;  //儿子的深搜序更大，说明是个孙子，跳过，返祖边只走一次
			//否则遇到了祖宗
			cov[re]++;    //标记这条边

			int pnt=u;//记录当前点

			int pe=pre_edge[u];//取出当前点的前驱边，准备走返祖边

			while(!even)    //若没有发现偶环，则尝试用两个奇环制造偶环
			{
				if(cov[pe])even++;//找到了一条被cov标记的前驱边，说明本奇环与其他奇环相交，可以制造偶环

				cov[pe]++;  // 标记沿途的边

				pnt=pre[pe];// 取出 边 的前驱点
				pe=pre_edge[pnt];  //再取 前驱点 的 前驱边

				if(pnt==v||pnt==-1||pe==0)break;  //若回到了祖宗v、找到不存在的前驱点或边，直接退出 
			}
		}
	}
	d[u].d=d[u].maxn+d[u].submaxn;  //求直径
}
void solve()
{
	//这里要补初始化//
	edge.clear();
    int n,m;
    cin>>n>>m;

    edge.push_back((node){-1,-1,0});//这个是凑数用的
    pre[0]=-1;

    for(int i=1;i<=m;i++)//存边
	{
    	int u,v;
    	cin>>u>>v;
    	edge.push_back((node){u,v,i});

    	e[u].push_back(i);  //注意:邻接表存边号
    	e[v].push_back(i);
	}

	for(int i=1;i<=n;i++)col[i]=-1,fa[i]=i;  //初始化
	col[0]=0;//0点不存在，但初始化为0备用

	for(int i=1;i<=n;i++)//深搜，找环，顺便求连通块的直径
	{
		if(!vis[i])dfs(i,0);
	}

	for(int i=1;i<=n;i++)
	{
		int f=find(i);
		siz[f]=max(siz[f],d[i].d);  //更新连通块的直径
	}

	if(n<4)cout<<-1<<"\n"; //n<4 is -1

	else if(odd&&even)  //若有奇又有偶，直接输出
	{
		cout<<0<<"\n";
	}
	else if(even)   //仅有偶环，只需再加一条边
	{
		cout<<1<<"\n";
	}
	else // 若无偶环，则需要尝试用若干直径制造偶环
	{
		map<int,int> is;
		priority_queue<int> q;
		for(int i=1;i<=n;i++)  //找连通块
		{
			int f=find(i);
			if(!is[f])
			{
				q.push(siz[f]);  //记录连通块的直径(最长路)
			}
			is[f]++;  //标记连通块
		}
		int temp=0;
		int ans=0;
		while(!q.empty())
		{
			temp+=q.top();   //取出一个直径，尝试连出一条长度不小于4的路
			q.pop();
			ans++;     //新边数量增加
			temp++;    //路长增加
			if(temp>=4)break;  //路长不小于4则退出，已经有偶环了
		}
		if(!odd)ans++;  //如果没有奇环，还需多加一条边
		cout<<ans<<"\n"; 
	}
}
/*
7 9
1 2
2 3
3 1
3 5
5 4
4 7
6 7
6 4
4 3
*/