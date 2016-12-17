// Powell.cpp: определяет точку входа для консольного приложения.
//
/**/
#include <math.h>
#include "nrutil.h"
#define TINY 1.0e-25
#define ITMAX 200
#define TOL 2.0e-4
#define GOLD 1.618034
#define GLIMIT 100.0
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define FMAX(a,b) (a>b)?a:b
#define CGOLD 0.3819660
#define SQR(a) a*a
#define SIGN(a,b) (b>0)?abs(a):(-abs(a))
#define ZEPS 1.0e-10

float brent(float ax, float bx, float cx, float (*f)(float), float tol,
            float *xmin)
{
    int iter;
    float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    float e=0.0;
    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(x);
    for (iter=1; iter<=ITMAX; iter++)
    {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a)))
        {
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1)
        {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=CGOLD*(e=(x >= xm ? a-x : b-x));

            else
            {
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
            }
        }
        else
        {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=(*f)(u);

        if (fu <= fx)
        {
            if(u >= x) a=x;
            else b=x;
            SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
        }
        else
        {
            if (u < x) a=u;
            else b=u;
            if (fu <= fw || w == x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v=u;
                fv=fu;
            }
        }
    }
    nrerror("Too many iterations in brent");
    *xmin=x;
    return fx;
}

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
            float (*func)(float))
{
    float ulim,u,r,q,fu,dum;
    *fa=(*func)(*ax);
    *fb=(*func)(*bx);
    if (*fb > *fa)
    {
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(*func)(*cx);
    while (*fb > *fc)
    {
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        if ((*bx-u)*(u-*cx) > 0.0)
        {
            fu=(*func)(u);
            if (fu < *fc)
            {
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            }
            else if (fu > *fb)
            {
                *cx=u;
                *fc=fu;
                return;
            }
            u=(*cx)+GOLD*(*cx-*bx);
        }
        else if ((*cx-u)*(u-ulim) > 0.0)
        {
            fu=(*func)(u);
            if (fu < *fc)
            {
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                SHFT(*fb,*fc,fu,(*func)(u))
            }
        }
        else if ((u-ulim)*(ulim-*cx) >= 0.0)
        {
            u=ulim;
            fu=(*func)(u);
        }
        else
        {
            fu=(*func)(u);
        }
        SHFT(*ax,*bx,*cx,u)
        SHFT(*fa,*fb,*fc,fu)
    }
}

void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret, float(*func)(float[]))
{
    void linmin(float p[], float xi[], int n, float *fret, float(*func)(float[]));
    int i, ibig, j;
    float del, fp, fptt, t, *pt, *ptt, *xit;

    pt = vector(1, n);
    ptt = vector(1, n);
    xit = vector(1, n);
    *fret = (*func)(p);
    for (j = 1; j <= n; j++) pt[j] = p[j];
    for (*iter = 1;; ++(*iter))
    {
        fp = (*fret);
        ibig = 0;
        del = 0.0;
        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n; j++) xit[j] = xi[j][i];
            fptt = (*fret);
            linmin(p, xit, n, fret, func);
            if (fptt - (*fret) > del)
            {
                del = fptt - (*fret);
                ibig = i;

            }
        }
        if (2.0*(fp - (*fret)) <= ftol*(fabs(fp) + fabs(*fret))+TINY)
        {
            free_vector(xit, 1, n);
            free_vector(ptt, 1, n);
            free_vector(pt, 1, n);
            return;

        }
        if (*iter == ITMAX) nrerror("powell exceeding maximum iterations");
        for (j = 1; j <= n; j++)
        {
            ptt[j] = 2.0*p[j] - pt[j];
            xit[j] = p[j] - pt[j];
            pt[j] = p[j];
        }
        fptt = (*func)(ptt);
        if (fptt < fp)
        {
            t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
            if (t < 0.0)
            {
                linmin(p, xit, n, fret, func);
                for (j = 1; j <= n; j++)
                {
                    xi[j][ibig] = xi[j][n];
                    xi[j][n] = xit[j];
                }
            }
        }
    }
}




int ncom;
float *pcom, *xicom, (*nrfunc)(float[]);
void linmin(float p[], float xi[], int n, float *fret, float(*func)(float[]))
{
    float brent(float ax, float bx, float cx,float(*f)(float), float tol, float *xmin);
    float f1dim(float x);
    void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,float *fc, float(*func)(float));
    int j;
    float xx, xmin, fx, fb, fa, bx, ax;
    ncom = n;
    pcom = vector(1, n);
    xicom = vector(1, n);
    nrfunc = func;
    for (j = 1; j <= n; j++)
    {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }
    ax = 0.0;
    xx = 1.0;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
    *fret = brent(ax, xx, bx, f1dim, TOL, &xmin);
    for (j = 1; j <= n; j++)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_vector(xicom, 1, n);
    free_vector(pcom, 1, n);
}

extern int ncom;
extern float *pcom, *xicom, (*nrfunc)(float[]);
float f1dim(float x)
{
    int j;
    float f, *xt;
    xt = vector(1, ncom);
    for (j = 1; j <= ncom; j++) xt[j] = pcom[j] + x*xicom[j];
    f = (*nrfunc)(xt);
    free_vector(xt, 1, ncom);
    return f;
}

int main()
{
    return 0;
}

