#include "stdafx.h"
#include "iostream"
#include <math.h>
#include <string>
//#include <map>
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
typedef std::vector<tuple<float*, int, float, function<float(float*)>>> funcTupleList;

float brent(float ax, float bx, float cx, float(*f)(float), float tol,
	float *xmin)
{
	int iter;
	float a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	float e = 0.0;
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

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
	float(*f1dim)(float))
{
	float ulim, u, r, q, fu, dum;
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
float *pcom, *xicom;

function<float(float*)> nrfunc;
void linmin(float p[], float xi[], int n, float *fret, function<float(float*)> func)
{
	float brent(float ax, float bx, float cx, float(*f)(float), float tol, float *xmin);
	float f1dim(float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float(*func)(float));
	int j;
	float xx, xmin, fx, fb, fa, bx, ax;
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
extern float *pcom, *xicom;
extern function<float(float*)> nrfunc;
float f1dim(float x)
{
	int j;
	float f, *xt;
	xt = fvector(1, ncom);
	for (j = 0; j < ncom; j++) xt[j] = pcom[j] + x*xicom[j];
	f = nrfunc(xt);
	free_vector(xt, 1, ncom);
	return f;
}

void powell_print(float *p, int n, float precision, float e, function<float(float*)> func)
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
		cout << func(p) - precision << endl;
	}
	else {
		cout << "f(p) = " << func(p) /*<< " (bad func)"*/ << endl;
		cout << abs(func(p) - precision) << endl;
	}
}


/*
Входные параметры :
p[1..n]-начальная точка
n - кол-во переменных
ftol - точность
func - сама функция
*/
void powell(float *p, const int n, float ftol, const float precision, function<float(float*)> func)
{
	float e = 0.000000001;
	int i, ibig, j; //переменные для циклов :)
					//float xi[2][2];
	float **xi = new float*[n];// двумерный динамический массив  массив содержащий начальный набор направлений , в теории единичная матрица  
	for (i = 0; i < n; ++i)
	{
		xi[i] = new float[n];
	}
	// заполняем массив
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j)
		{
			if (i == j) { xi[i][j] = 1; }
			else xi[i][j] = 0;
		}
	}

	float del, fp, fptt, t, *pt, *ptt, *xit, krit; // del-дельта , fp и fptt - ф-ция в р и рtt

	int iter = 0;
	pt = fvector(1, n);
	ptt = fvector(1, n);
	xit = fvector(1, n);
	float fret;
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
			powell_print(p, n, precision, e, func);
			return;
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
				if (((p[i] - pt[i]))<ftol*pt[i]+ftol) ++krit;
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
				powell_print(p, n, precision, e, func);
				return;
			}

		}
	}


}

int main()
{
	int iter = 0;
	float ftol = 0.001;

	funcTupleList funcVec = { 
		//std::make_tuple(new float[2]{ 8, 9 }, 2, 0, [](float *p) -> float { return 4 * (p[0] - 1)*(p[0] - 1) + (p[1] - 6)*(p[1] - 6);}),
		 		//std::make_tuple(new float[2]{ 1,2 }, 2, 0, [](float *p) -> float { return 1 / ((p[0] + 1)*(p[0] + 1) + p[1] * p[1]);}),
				//std::make_tuple(new float[3]{ -1, 2, -2}, 3, 0, [](float *p) -> float { return p[0] * p[0] + 3 * p[1] * p[1] + 5 * p[2] * p[2];}),
				//std::make_tuple(new float[4]{ -1, 2, -2, 1 }, 4, 0, [](float *p) -> float { return p[0] * p[0] + 3 * p[1] * p[1] + 5 * p[2] * p[2]+4*(p[3]*p[3]);}),
		 		//std::make_tuple(new float[4]{ -1, 2, -1, 1 }, 4, 0, [](float *p) -> float { return p[0] * p[0] + 3 * p[1] * p[1] + 5 * p[2] * p[2] + (p[3] - 1)*(p[3] - 1);}),

				//std::make_tuple(new float[2]{ 1, -2}, 2, 0, [](float *p) -> float { return 2*p[0] * p[0] +  p[1] * p[1] +  p[1] * p[2];}),
				//std::make_tuple(new float[2]{ -1, 1 }, 2, 0, [](float *p) -> float { return (p[1] - p[0] * p[0])*(p[1] - p[0] * p[0]) + 100 * (1 - p[0])*(1 - p[0]);}),
				//std::make_tuple(new float[2]{ -1, 1 }, 2, 0, [](float *p) -> float { return 100 * (p[1] - p[0] * p[0] * p[0])*(p[1] - p[0] * p[0] * p[0]) + (1 - p[0])*(1 - p[0]);}),
			//std::make_tuple(new float[2]{ -1, 1 }, 2, 0, [](float *p) -> float { return exp(p[1] * p[1] - p[0] * p[0])*(2 * p[0] * p[0] + 3 * p[1] * p[1]);}),

		 		//std::make_tuple(new float[2]{ 1, 2 }, 2, 0, [](float *p) -> float { return (tan(p[1] - p[0]))*(tan(p[1] - p[0]));}),
			//	std::make_tuple(new float[2]{ 1, 2 }, 2, 0, [](float *p) -> float { return (sin(p[0] - p[1]))*(cos(p[0] - p[1]));}),
		 		//std::make_tuple(new float[4]{1, 2, 2, 2 }, 4, 0, [](float *p) -> float { return p[0]*p[0]+100 * (p[1] - p[2])*(p[1] - p[2]) + (tan(p[2] - p[3]))*(tan(p[2] - p[3]));}),
		 		//std::make_tuple(new float[2]{ 3, -2 }, 2,0,  [](float *p) -> float { return (p[0] * p[0] + p[1] - 11)*(p[0] * p[0] + p[1] - 11) + (p[0] + p[1] * p[1] - 7)*(p[0] + p[1] * p[1] - 7);}),
		 		//std::make_tuple(new float[2]{ 1, 1 }, 2, 0, [](float *p) -> float { return (p[0] * p[0] + p[1] - 11)*(p[0] * p[0] + p[1] - 11) + (p[0] + p[1] * p[1] - 7)*(p[0] + p[1] * p[1] - 7);}),

		 		//std::make_tuple(new float[3]{ 2, 1, 1 }, 3, 0.5, [](float *p) -> float { return p[2]*p[2]+p[0]*p[0]+(p[0]-1)*(p[0]-1)+(1-p[1])*(1 - p[1]);}),
				

		 		
		std::make_tuple(new float[2]{ 8, 9 }, 2, 0, [](float *p) -> float { return 4 * (p[0] - 1)*(p[0] - 1) + (p[1] - 6)*(p[1] - 6)+10;})
	
	};

	for (auto const& item : funcVec) {
		powell(get<0>(item), get<1>(item), ftol, get<2>(item), get<3>(item));
	}

	system("pause");

	return 0;
};