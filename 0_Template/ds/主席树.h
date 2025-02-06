#define int long long

const int N=2e5+10;
#define lc(x) tr[x].l
#define rc(x) tr[x].r
struct node{
	int l,r,s; //l,r->son  s->cnt 
}tr[N*20];
int root[N],idx;
int n,m,a[N];
vector<int> b;
void build(int &x,int l,int r){
	x=++idx;
	if(l==r) return;
	int m=l+r>>1;
	build(lc(x),l,m);
	build(rc(x),m+1,r);
}
void insert(int x,int &y,int l,int r,int k){
	y=++idx; tr[y]=tr[x]; tr[y].s++;
	if(l==r) return;
	int m=l+r>>1;
	if(k<=m) insert(lc(x),lc(y),l,m,k);
	else insert(rc(x),rc(y),m+1,r,k);
}
int query(int x,int y,int l,int r,int k){
	if(l==r) return l;
	int m=l+r>>1;
	int s=tr[lc(y)].s-tr[lc(x)].s;
	if(k<=s) return query(lc(x),lc(y),l,m,k);
	else return query(rc(x),rc(y),m+1,r,k-s);
}

void solve(){
	int n,m;
	cin>>n>>m;
	for(int i=1;i<=n;i++){
		cin>>a[i];
		b.push_back(a[i]);
	}
	sort(b.begin(),b.end());
	b.erase(unique(b.begin(),b.end()),b.end());
	int bn=b.size();
	build(root[0],1,bn);
	for(int i=1;i<=n;i++){
		int id=lower_bound(b.begin(),b.end(),a[i])-b.begin()+1;
		insert(root[i-1],root[i],1,bn,id);
	}
	while(m--){
		int l,r,k;
		cin>>l>>r>>k;
		int id=query(root[l-1],root[r],1,bn,k)-1;
		cout<<b[id]<<"\n";
	}
}

