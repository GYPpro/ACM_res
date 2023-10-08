#include <template_overAll.h>


template <class T>
class stackElememt{
public:
    stackElememt *top;
    stackElememt *base;
    T data;
};

template <class T>
class mystack{
private:
    // stackElememt<T> *bedrock;
    stackElememt<T> *upper;
public:
    mystack()
    {
        stackElememt<T> *st = new stackElememt<T>;
        // bedrock = st;
        // bedrock->base = bedrock;
        // bedrock->top = bedrock;
        upper = st;
        upper->base = nullptr;
        upper->top = upper;
    }
    void push(T _t)
    {
        stackElememt<T> *st = new stackElememt<T>;
        st->data = _t;
        upper->top = st;
        st->top = st;
        st->base = upper;
        upper = upper->top;
    }
    T pop()
    {
        T rt = upper->data;
        stackElememt<T> *st = upper;
        upper->base->top = upper->base;
        upper = upper->base;
        free(st);
        return rt;
    }
    ~mystack()
    {
        for(;;)
        {
            if(upper == nullptr) break;
            stackElememt<T> *st = upper;
            upper = upper->base;
            free(st);
        }
    }
};
