//区间修改，区间查询 如果N=4000，内存大概在300多MB 注意空间 
#define int long long
const int N=2050;
int t1[N][N],t2[N][N],t3[N][N],t4[N][N];
int n,m;
int lowbit(int x){
	return x&-x;
}
void add(int x,int y,int val){
	for(int i=x;i<=n;i+=lowbit(i)){
		for(int j=y;j<=m;j+=lowbit(j)){
			t1[i][j]+=val;
			t2[i][j]+=val*x;
			t3[i][j]+=val*y;
			t4[i][j]+=val*x*y;
		}
	}
}
void range_add(int xa,int ya,int xb,int yb,int val){
	add(xa,ya,val); add(xa,yb+1,-val);
	add(xb+1,ya,-val); add(xb+1,yb+1,val);
}
int ask(int x,int y){
	int ans=0;
	for(int i=x;i;i-=lowbit(i)){
		for(int j=y;j;j-=lowbit(j)){
			ans+=(x+1)*(y+1)*t1[i][j]-(y+1)*t2[i][j]-(x+1)*t3[i][j]+t4[i][j];
		}
	}
	return ans;
}
int range_ask(int xa,int ya,int xb,int yb){
	return ask(xb,yb)-ask(xb,ya-1)-ask(xa-1,yb)+ask(xa-1,ya-1);
}
void solve(){
	cin>>n>>m;
	int op;
	while(cin>>op){
		if(op==1){
			int a,b,c,d,val;
			cin>>a>>b>>c>>d>>val;
			range_add(a,b,c,d,val);
		}
		else{
			int a,b,c,d;
			cin>>a>>b>>c>>d;
			cout<<range_ask(a,b,c,d)<<"\n";
		}
	}
}


// 点修 区间查询 可以建一个数组，可以开更大的N，防止MLE 
#define int long long
const int N=4100;
int t[N][N];
int n,m;
int lowbit(int x){
	return x&-x;
}
void add(int x,int y,int val){
	for(int i=x;i<=n;i+=lowbit(i)){
		for(int j=y;j<=m;j+=lowbit(j)){
			t[i][j]+=val;
		}
	}
} 
int sum(int x,int y){
	int ans=0;
	for(int i=x;i;i-=lowbit(i)){
		for(int j=y;j;j-=lowbit(j)){
			ans+=t[i][j];
		}
	}
	return ans;
}
int ask(int x1,int y1,int x2,int y2){
	return sum(x2,y2)-sum(x2,y1-1)-sum(x1-1,y2)+sum(x1-1,y1-1);
}
void solve(){
	cin>>n>>m;
	int op;
	while(cin>>op){
		if(op==1){
			int x,y,val;
			cin>>x>>y>>val;
			add(x,y,val);
		}
		else{
			int a,b,c,d;
			cin>>a>>b>>c>>d;
			cout<<ask(a,b,c,d)<<"\n";
		}
	}
}

// 第k大
int t[N]; 
int kth(int k){
	int ans=0,x=0;
	for(int i=log2(n);~i;i--){
		x+=(1<<i);
		if(x>=n||sum+t[x]>=k) x-=(1<<i);
		else sum+=t[x];
	}
	return x+1;
}

// 区间修改 单点查询
const int N=4100;
int t[N][N];
int n,m;
int lowbit(int x){
	return x&-x;
}
void add(int x,int y,int val){
	for(int i=x;i<=n;i+=lowbit(i)){
		for(int j=y;j<=m;j+=lowbit(j)){
			t[i][j]+=val; 
		}
	}
} 
void range_add(int x1,int y1,int x2,int y2,int val){
	add(x1,y1,val); add(x1,y2+1,-val);
	add(x2+1,y1,-val); add(x2+1,y2+1,val);
}
int sum(int x,int y){
	int ans=0;
	for(int i=x;i;i-=lowbit(i)){
		for(int j=y;j;j-=lowbit(j)){
			ans+=t[i][j];
		}
	}
	return ans;
}
