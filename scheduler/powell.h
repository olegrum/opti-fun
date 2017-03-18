#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <functional>
#include <utility>
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


using namespace std;
typedef std::vector<tuple<double*, int, double, function<double(double*)>>> funcTupleList;

double brent(double ax, double bx, double cx, double(*f)(double), double tol,
	double *xmin)
{
	int iter;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;
	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = (*f)(x);
	for (iter = 1; iter <= ITMAX; iter++) //ITMAX
	{
		xm = 0.5*(a + b);
		tol2 = 2.0*(tol1 = tol*fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5*(b - a)))
		{
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1)
		{
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2.0*(q - r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a - x) || p >= q*(b - x))
				d = CGOLD*(e = (x >= xm ? a - x : b - x));

			else
			{
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		}
		else
		{
			d = CGOLD*(e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = (*f)(u);

		if (fu <= fx)
		{
			if (u >= x) a = x;
			else b = x;
			SHFT(v, w, x, u)
				SHFT(fv, fw, fx, fu)
		}
		else
		{
			if (u < x) a = u;
			else b = u;
			if (fu <= fw || w == x)
			{
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w)
			{
				v = u;
				fv = fu;
			}
		}
	}
	//nrerror("Too many iterations in brent");
	*xmin = x;
	return fx;
}

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	double(*f1dim)(double))
{
	double ulim, u, r, q, fu, dum;
	*fa = (*f1dim)(*ax);
	*fb = (*f1dim)(*bx);
	if (*fb > *fa)
	{
		SHFT(dum, *ax, *bx, dum)
			SHFT(dum, *fb, *fa, dum)
	}
	*cx = (*bx) + GOLD*(*bx - *ax);
	*fc = (*f1dim)(*cx);
	while (*fb > *fc)
	{
		r = (*bx - *ax)*(*fb - *fc);
		q = (*bx - *cx)*(*fb - *fa);
		u = (*bx) - ((*bx - *cx)*q - (*bx - *ax)*r) / (2.0*SIGN(FMAX(fabs(q - r), TINY), q - r));
		ulim = (*bx) + GLIMIT*(*cx - *bx);
		if ((*bx - u)*(u - *cx) > 0.0)
		{
			fu = (*f1dim)(u);
			if (fu < *fc)
			{
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;
			}
			else if (fu > *fb)
			{
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD*(*cx - *bx);
		}
		else if ((*cx - u)*(u - ulim) > 0.0)
		{
			fu = (*f1dim)(u);
			if (fu < *fc)
			{
				SHFT(*bx, *cx, u, *cx + GOLD*(*cx - *bx))
					SHFT(*fb, *fc, fu, (*f1dim)(u))
			}
		}
		else if ((u - ulim)*(ulim - *cx) >= 0.0)
		{
			u = ulim;
			fu = (*f1dim)(u);
		}
		else
		{
			fu = (*f1dim)(u);
		}
		SHFT(*ax, *bx, *cx, u)
			SHFT(*fa, *fb, *fc, fu)
	}
}

int ncom;
double *pcom, *xicom;

function<double(double*)> nrfunc;
void linmin(double p[], double xi[], int n, double *fret, function<double(double*)> func)
{
	double brent(double ax, double bx, double cx, double(*f)(double), double tol, double *xmin);
	double f1dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double(*func)(double));
	int j;
	double xx, xmin, fx, fb, fa, bx, ax;
	ncom = n;
	pcom = fvector(1, n);
	xicom = fvector(1, n);
	nrfunc = func;
	for (j = 0; j < n; j++)
	{
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	ax = 0.0;
	bx = 1.0;
	mnbrak(&ax, &bx, &xx, &fa, &fx, &fb, f1dim);//золотое сечение
	*fret = brent(ax, bx, xx, f1dim, TOL, &xmin);
	for (j = 0; j <n; j++)
	{
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom, 1, n);
	free_vector(pcom, 1, n);
}
extern int ncom;
extern double *pcom, *xicom;
extern function<double(double*)> nrfunc;
double f1dim(double x)
{
	int j;
	double f, *xt;
	xt = fvector(1, ncom);
	for (j = 0; j < ncom; j++) xt[j] = pcom[j] + x*xicom[j];
	f = nrfunc(xt);
	free_vector(xt, 1, ncom);
	return f;
}

/*void powell_print(double *p, int n, double precision, double e, function<double(double*)> func)
{
	for (int i = 0; i < n; ++i) {
		//if (abs(p[i] - precision) < e) {
			cout << "p[" << i << "] : " << p[i] << endl;
		}
		//else {
		//cout << "p[" << i << "] : " << p[i] << " (bad func)" << endl;
		//}
	//}
	if (abs(func(p) - precision) < e) {
		cout << "f(p) = " << func(p) << endl;
	}
	else {
		cout << "f(p) = " << func(p) << " (bad func)" << endl;
	}
}*/


/*
Входные параметры :
p[1..n]-начальная точка
n - кол-во переменных
ftol - точность
func - сама функция
*/
double powell(double *p, const int n, double ftol, function<double(double*)> func)
{

	double e = 0.00001;
	int i, ibig, j; //переменные для циклов :)
				//double xi[2][2];
	double **xi = new double*[n];// двумерный динамический массив  массив содержащий начальный набор направлений , в теории единичная матрица
	for (i = 0; i < n; ++i)
	{
		xi[i] = new double[n];
	}
	// заполняем массив
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j)
		{
			if (i == j) { xi[i][j] = 1; }
			else xi[i][j] = 0;
		}
	}

	double del, fp, fptt, t, *pt, *ptt, *xit, krit; // del-дельта , fp и fptt - ф-ция в р и рtt

	int iter = 0;
	pt = fvector(1, n);
	ptt = fvector(1, n);
	xit = fvector(1, n);
	double fret;
	fret = func(p);


	for (j = 0; j < n; j++) pt[j] = p[j]; //сохраняем начальную точку


	for (iter;; ++iter)
	{
		fp = fret;
		ibig = 0; // номер направления
		del = 0.0;   // разница в предыдущем иксе и нынешнем
		for (i = 0; i < n; i++) // в каждой итерации , цикл по всем направлениям
		{
			for (j = 0; j < n; j++) xit[j] = xi[j][i]; // копируем направления
			fptt = fret;
			linmin(p, xit, n, &fret, func);  //минимизируем
			if (fptt - fret > del) // записываем del, если это найбольшее уменьшение
			{
				del = fptt - fret;
				ibig = i;

			}
		}
		if (fret == fptt) {


			free_vector(xit, 1, n);
			free_vector(ptt, 1, n);
			free_vector(pt, 1, n);
			//время освободить память от двумерного массива
			for (i = 0; i < n; ++i)
			{
				delete[]xi[i];
			}
			//красивый вывод результатов ))))))
			// 			for (i = 0; i < n; ++i) {
			// 				cout << "p[" << i << "] : " << p[i] << endl;
			// 			}
			// 			cout << "f(p) = " << (func)(p) << endl;
			//powell_print(p, n, precision, e, func);
			return func(p);
		}
		if (iter == ITMAX) nrerror("powell exceeding maximum iterations");
		for (j = 0; j < n; j++) // меняем направление и  сохраняем старую точку
		{
			ptt[j] = 2.0*p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
		fptt = func(ptt); // значение функции в новой точке

		if ((fptt<fp) || (2.0*(fp - 2.0*fret + fptt)*SQR(fp - fret - del) - del*SQR(fp - fptt)<0))

		{
			linmin(p, xit, n, &fret, func); // переходим к минимуму нового направления и сохраняем его
			for (j = 0; j < n; j++)
			{
				xi[j][ibig] = xi[j][n];
				xi[j][n] = xit[j];
			}

		}
		if (abs(fret)<0.000001) //критерий завершения программы
		{
			krit = 0;

			for (i = 0; i < n;++i)
			{
				if (((p[i] - pt[i]))<ftol*pt[i]) ++krit;
			}
			if (krit == n)
			{
				free_vector(xit, 1, n);
				free_vector(ptt, 1, n);
				free_vector(pt, 1, n);
				//время освободить память от двумерного массива
				for (i = 0; i < n; ++i)
				{
					delete[]xi[i];
				}
				//красивый вывод результатов ))))))
				// 				for (i = 0; i < n; ++i) {
				// 					if ((p[i] - precision) < e) {
				// 						cout << "p[" << i << "] : " << p[i] << endl;
				// 					} else {
				// 						cout << "p[" << i << "] : " << p[i] << " (bad func)" << endl;
				// 					}

				// 				}
				// 				cout << "f(p) = " << func(p) << endl;
				//powell_print(p, n, precision, e, func);
				return func(p);
			}

		}
	}


}
