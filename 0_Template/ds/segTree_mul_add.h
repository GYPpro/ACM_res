/*
三种操作：（模数为mod） 
1.区间乘val 
2.区间加val
3.区间和 
*/
#define int long long 
#define lc p<<1
#define rc p<<1|1
const int N=1e5+10;
int mod=1e9+7;
struct node{
	int l,r,sum,mul,add;
}tr[4*N];
vector<int> a(N);
void pushup(int p){
	tr[p].sum=tr[lc].sum+tr[rc].sum;
	tr[p].sum%=mod;
}
void build(int p,int l,int r){
	tr[p].l=l;tr[p].r=r; 
	tr[p].mul=1;tr[p].add=0;
	if(l==r){
		tr[p].sum=a[l];
		return;
	}
	int m=l+r>>1;
	build(lc,l,m);
	build(rc,m+1,r);
	pushup(p);
}
void pushdown(int p){
	tr[lc].sum=(tr[lc].sum*tr[p].mul+tr[p].add*(tr[lc].r-tr[lc].l+1))%mod;
	tr[rc].sum=(tr[rc].sum*tr[p].mul+tr[p].add*(tr[rc].r-tr[rc].l+1))%mod;
	tr[lc].mul=tr[lc].mul*tr[p].mul%mod;
	tr[rc].mul=tr[rc].mul*tr[p].mul%mod;
	tr[lc].add=(tr[lc].add*tr[p].mul+tr[p].add)%mod;
	tr[rc].add=(tr[rc].add*tr[p].mul+tr[p].add)%mod;
	tr[p].mul=1;tr[p].add=0;
}
void update_mul(int p,int l,int r,int val){
	if(l<=tr[p].l&&tr[p].r<=r){
		tr[p].sum=tr[p].sum*val%mod;
		tr[p].mul=tr[p].mul*val%mod;
		tr[p].add=tr[p].add*val%mod;
		return;
	}
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	if(l<=m) update_mul(lc,l,r,val);
	if(m<r) update_mul(rc,l,r,val);
	pushup(p);
	return;
}
void update_add(int p,int l,int r,int val){
	if(l<=tr[p].l&&tr[p].r<=r){
		tr[p].add=(tr[p].add+val)%mod;
		tr[p].sum=(tr[p].sum+val*(tr[p].r-tr[p].l+1))%mod;
		return;
	}
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	if(l<=m) update_add(lc,l,r,val);
	if(m<r) update_add(rc,l,r,val);
	pushup(p);
	return;
}
int query(int p,int l,int r){
	if(l<=tr[p].l&&tr[p].r<=r) return tr[p].sum;
	pushdown(p);
	int m=tr[p].l+tr[p].r>>1;
	int sum=0;
	if(l<=m) sum+=query(lc,l,r);
	if(m<r) sum+=query(rc,l,r);
	return sum%mod;
}
void solve(){
	int n,m;
	cin>>n>>m>>mod;
	for(int i=1;i<=n;i++) cin>>a[i];
	build(1,1,n);
	while(m--){
		int op;
		cin>>op;
		if(op==1){
			int x,y,val;
			cin>>x>>y>>val;
			update_mul(1,x,y,val);
		}
		else if(op==2){
			int x,y,val;
			cin>>x>>y>>val;
			update_add(1,x,y,val);
		}
		else{
			int x,y;
			cin>>x>>y;
			cout<<query(1,x,y)<<"\n";
		}
	}
	return;
}
