// ##IGNORE##
template<class TYPE_NAME>
class profixSum //HACK
{// 下标从0、1开始均可
private:
    vector<TYPE_NAME> initdata;
    vector<TYPE_NAME> sum;
public:
    profixSum(vector<TYPE_NAME> _inidata)
    {
        initdata = _inidata;
        sum.resize(initdata.size());
        sum[0] = initdata[0];
        for(int i = 1;i < initdata.size();i ++)
            sum[i] = sum[i-1] + initdata[i];
    }
    profixSum()
    {}
    profixSum(int n)
    {
        initdata.resize(n);
        sum.resize(n);
    }

    void push_back(TYPE_NAME _x)
    {
        initdata.push_back(_x);
        sum.push_back(_x + sum[sum.size()-1]);
    }

    void reset(vector<TYPE_NAME> _inidata)
    {
        initdata = _inidata;
        sum.resize(initdata.size());
        sum[0] = initdata[0];
        for(int i = 1;i < initdata.size();i ++)
            sum[i] = sum[i-1] + initdata[i];
    }

    TYPE_NAME getSum(int l,int r)
    {
        if(l == 0) return sum[r];
        else
        return sum[r]-sum[l-1]
    }

};
