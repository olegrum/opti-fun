#include <cmath>
#include <vector>
#include <algorithm>


#define DIM 30
#define MAXDIM 1111
#define MAXDEG 13


using namespace std;

int POLY[MAXDIM];
int VINIT[MAXDIM][MAXDEG];
double RECIPD;
int LASTQ[MAXDIM], V[MAXDIM][DIM];
int I;
int MAXCOL;
int COUNT, TAUS;
double QUASI[MAXDIM]; // координаты одной точки
double ANSWER[1048576][DIM]; // матрица, содержащия кординаты точек сетки

/*
Массив POLY выдает последовательные простые полиномы в двоичной записи,
Например: 45 = 100101; имеет единицы на 5, 2, и 0 местах (Подсчитывая справа налево)
и поэтому представляет: X**5 + X**2 + X**0
Эти полиномы испульзуются в работах Соболя в СССР вычислительной математики и физике 16 (1977),
236-242. Более завершенная таблица была дана Соболем и Левитаном, поиском точек равномерно
заполняющих н-мерный гиперкуб (В России), Препринт Академия наук СССР, 40, Москва 1976.

Инициализация массива VINIT из тестируемого списка. Для многочлена степени M, M исходных
значений необходимо: эти значения приводятся здесь. Последующие значения расчитываются в "INSOBL".
*/

struct point
{
    int point;
    double value;
};


void f(){
    int w[13] = {0, 3, 4, 6, 8, 14, 20, 38, 54, 102, 162, 338, 482};
    ifstream fin("data.txt");
    for (int i = 0; i < MAXDIM - 1; ++i)
        fin>>POLY[i];
    for (int j = 0; j < MAXDIM; ++j)
            VINIT[j][0] = 1;
    for (int i = 1; i < MAXDEG; ++i)
        for (int j = w[i] - 1; j < MAXDIM; ++j)
            fin>>VINIT[j][i];
    fin.close();
}


/*
Сначала проверяются заданные пользователем значение "DIMEN" - размерность квазислучайных
веторов строго между 0 и 30.
Затем проверяется "ATMOST", верхняя числовая грань вызовов "GOSOBL", сделанных пользователем.
Они должны быть положительны и меньше 2**30.
(Мы предполагаем, что работаем на компьютере с минимальной разрядностью 31 битов)
Число столбцов матрицы V, которые инициализируются как MAXCOL = число разрядов в ATMOST.
В "GOSOBL" мы проверяем, что их число не превышается.
Ведущие элементы каждого ряда V инициализируются, используя "VINIT" из "BDSOBL".
Каждая строчка соответствует простому многочлену(см "BDSOBL").
Если многочлен степени M, элементы после первого M рассчитываются.
Числа в  V на самом деле являются двоичными фукциями.
"RECIPD" содержит 1/(общий знаменатель их всех).
"INSOBL" неявно вычисляет все первые нулевые векторы, но не возвращает их на вызов программы.
Следующие векторы берутся из "GOSOBL".
"LASTQ" хранит числители последнего сгенерированного вектора.
"TAUS" для определения "благоприятных" значений.
Как и говорилось в BRATLEY/FOX, они имеют вид N = 2**K ,где K по определению (TAUS+S-1)
для интеграциии K и большие TAUS для глобальной оптимизации.
Входные данные :
Из программы пользователя : DIMEN, ATMOST
для блока DATA "BDSOBL" : POLY, VINIT
Выходные данные :
Для программы пользователя : FLAG, TAUS
Для "GOSOBL" VIA /SOBOL/ : V, S, MAXCOL, COUNT, LASTQ, RECIPD
*/
int INSOBL(int DIMEN, int ATMOST)
{
    int const MAXBIT = 30;
    int L, M, NEWV;
    int TAU[12] = {0, 1, 3, 5, 8, 11, 15, 19, 23, 27, 31, 35};
    char INCLUD[MAXDEG];
    //Проверка параметров:
    if (!((DIMEN >= 1) && (DIMEN <= MAXDIM) && (ATMOST > 0) && (ATMOST < pow(2,MAXBIT))))
    {
        cout << "error"<<endl;
        return 1;
    }
    (DIMEN < MAXDEG)? TAUS = TAU[DIMEN-1] : TAUS = -1;

    int I = ATMOST;
    MAXCOL = int(round(log2(I))); //Найти число разрядов в ATMOST
   /* while (I)
    {
        MAXCOL++;
        I >>= 1;
    }
*/
    //Инициализация первой строки матрицы V:
    for (int I = 0; I < MAXCOL; ++I)
        V[0][I] = 1;

    //Инициализация остальных строк матрицы V:
    for (int I = 1; I <= DIMEN; ++I)
    {
        //Комбинация разрядов полинома I(смотрите комментарии "BDSOBL")
        //найти степень многочлена I из бинарного кодирования.

        int J = POLY[I-1];
        M = int(round(log2(J)))-1;
        //J >>= 1;
      /*  while (J)
        {
            J >>= 1;
            ++M;
        }*/
        //Мы раскладываем этот бит на компоненты
        J = POLY[I-1];
        for (int K = M - 1; K >= 0; --K)
        {
            INCLUD[K] = J % 2;
            J >>= 1;
        }

        //Ведущие элементы ряда I берутся из VINIT
        for (int J = 0; J < M; ++J)
            V[I][J] = VINIT[I][J];

        //Вычисляем остальные элементы ряда I
        for (int J = M; J < MAXCOL; ++J)
        {
            NEWV = V[I][J - M];
            L = 1;

            for (int K = 0; K < M; ++K)
            {
                L <<= 1;
                if (INCLUD[K])
                    NEWV ^= (L*V[I][J - K-1]);
            }
            V[I][J] = NEWV;
        }
    }

    //Перемножаем строки V на соответствующие степени двойки
    L = 1;
    for (int J = MAXCOL - 2; J >= 0; --J)
    {
        L <<= 1;
        for( I = 0; I < DIMEN; ++I)
            V[I][J] = V[I][J]*L;
    }
    RECIPD = 1.0/(L<<1); //RECIPD IS 1/(общий знаменатель элементов в V)
    COUNT = 1;
    //Находим первый вектор и его значение для "GOSOBL"
    for(I = 0; I < DIMEN; ++I)
        LASTQ[I] = 0;

    return 0;
}


/*
Эта подпрограмма "GOSOBL" создает квазислучайные векторы с каждым вызовом
Пользователь должен вызывать "INSOBL" до вызова "GOSOBL".  После вызова "INSOBL",
Входные данные:
V        Таблица направляющих чисел
S        Размерность
MAXCOL   последня строчка V, которая будет использована
COUNT    Порядковый номер вызова
LASTQ    Числители для последнего сгенерированного вектора
RECIPD   (1/знаменатель) для этих числителей находим позицию справа от нуля в COUNT
*/
int GOSOBL(int c, int DIMEN)
{
    int I, L;
    L = 1;
    I = c;
    while(I % 2)
    {
       ++L;
       I >>= 1;
    }
    --L;
    if (L > MAXCOL)
        {
            cout<<" TOO MANY CALLS ON GOSOBL"<< endl;
            return 1;
        }

    //Расчет новых компонентов QUASI
    for(I = 0; I < DIMEN; I++)
        QUASI[I] = (LASTQ[I]^=V[I][L])*RECIPD;
    return 0;
}





int cmp(const void *a, const void *b){
    struct point *x = (struct point *) a;
    struct point *y = (struct point *) b;
    return 2*((x->value - y->value)>0)-1;
}


// входные параметры:
// n - размерность пространства
// a,b - начало и конец ребра n мерного куба
// m - количество точек, которые хотим получить на выходе
double** tms(int n, double a, double b, int m, double(*func)(double[]))
{
    if ((2 > n)&&(n > DIM))
    {
        cout<<"Incorrect input space dimension";
        return NULL;
    }
    if ((m > n)&&(m < 1))
    {
        cout<<"Incorrect input the number of output points";
        return NULL;
    }
        if (b <= a)
    {
        cout<<"Incorrect input a, b";
        return NULL;
    }

    int r[DIM-2] = {32, 64, 32, 16, 8, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    int k = pow(r[n-2],n); // количество точек в сетке

    f();

   // ofstream fout("out_all.txt"); // файл, в котором хранятся k точек сетки
   // ofstream fout1("out_min.txt"); // файл, в котором хранятся все m точек, где функция минимальна
    //values = (struct point *) malloc(k*sizeof(struct point));

    struct point *values = new struct point[k+1];

    INSOBL(n,k);

    for (int i = 0; i < n; ++i)
    {
        QUASI[i] = a;
        //fout << QUASI[i] << " ";
        ANSWER[0][i] = QUASI[i];
    }
    //fout << endl;
    values[0].value = (*func)(QUASI);
    values[0].point = 0;
    for (int j = 0; j < k - 1; ++j)
    {
        GOSOBL(j,n); // формируем одну точку сетки, ее координаты записываются в массив QUASI
        for (int i = 0; i < n; ++i)
        {
            ANSWER[j+1][i] = QUASI[i] = a + (b - a) * QUASI[i]; // записываем в матрицу координаты точки
            //fout << QUASI[i] << " "; // записываем в файл координаты точки
        }
       // fout<<endl;
        values[j+1].value = (*func)(QUASI); //записываем в вектор значение функции в данной точке
        values[j+1].point = j+1; //записываем в вектор номер данной точки
    }
    //fout.close();

    double** ANSWER1;
    try{
        ANSWER1 = new double * [m];
        for (int i=0; i<m; i++)
            ANSWER1[i] = new double[n];
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }
    qsort(values, k, sizeof(struct point), cmp); // сортируем вектор значений функции и запоминаем какой точке принадлежит каждое значение
    for (int p = 0; p < m; ++p)
        for (int t = 0; t < n; ++t)
            ANSWER1[p][t] = ANSWER[values[p].point][t]; // записываем в файл координаты точек, в которых функция принимает минимальное значение

    //fout1.close();
    delete values;
    return ANSWER1;
}




