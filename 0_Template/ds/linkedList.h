#include <template_overAll.h>

#define NEL 1145141919810LL

//HACK
template<class T>
class linkedList{
    class LLnode{
        public:
        T data;
        int r;
        int l;

        LLnode()
        {
            r = NEL;
            l = NEL;
        }
        LLnode(T _data)
        {
            r = NEL;
            l = NEL;
            data = _data;
        }
    };

    public:
    vector<LLnode> dt;
    void init(const vector<T> &indata)
    {
        dt.clear;
        dt.resize(indata.size());
        for(int i = 1;i < indata.size()-1;i ++)
        {
            dt[i-1].r = i;
            dt[i+1].l = i;
            dt[i].l = i-1;
            dt[i].r = i+1;
        }
        for(int i = 0;i < indata.size()-1;i ++)
            dt[i].data = indata[i];
    }
    linkedList(){;}
    linkedList(const vector<T> &indata)
    {
        dt.clear;
        dt.resize(indata.size());
        for(int i = 1;i < indata.size()-1;i ++)
        {
            dt[i-1].r = i;
            dt[i+1].l = i;
            dt[i].l = i-1;
            dt[i].r = i+1;
        }
        for(int i = 0;i < indata.size()-1;i ++)
            dt[i].data = indata[i];
    }

    T& operator[] (const int i)
    {
        return dt[i];
    }

    void rm(int i)
    {
        if(dt[i].r != NEL) dt[dt[i].r].l = dt[i].l;
        if(dt[i].l != NEL) dt[dt[i].l].r = dt[i].r;
    }
};