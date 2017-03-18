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
double QUASI[MAXDIM]; // ���������� ����� �����
double ANSWER[1048576][DIM]; // �������, ���������� ��������� ����� �����

/*
������ POLY ������ ���������������� ������� �������� � �������� ������,
��������: 45 = 100101; ����� ������� �� 5, 2, � 0 ������ (����������� ������ ������)
� ������� ������������: X**5 + X**2 + X**0
��� �������� ������������ � ������� ������ � ���� �������������� ���������� � ������ 16 (1977),
236-242. ����� ����������� ������� ���� ���� ������� � ���������, ������� ����� ����������
����������� �-������ �������� (� ������), �������� �������� ���� ����, 40, ������ 1976.

������������� ������� VINIT �� ������������ ������. ��� ���������� ������� M, M ��������
�������� ����������: ��� �������� ���������� �����. ����������� �������� ������������� � "INSOBL".
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
������� ����������� �������� ������������� �������� "DIMEN" - ����������� ��������������
������� ������ ����� 0 � 30.
����� ����������� "ATMOST", ������� �������� ����� ������� "GOSOBL", ��������� �������������.
��� ������ ���� ������������ � ������ 2**30.
(�� ������������, ��� �������� �� ���������� � ����������� ������������ 31 �����)
����� �������� ������� V, ������� ���������������� ��� MAXCOL = ����� �������� � ATMOST.
� "GOSOBL" �� ���������, ��� �� ����� �� �����������.
������� �������� ������� ���� V ����������������, ��������� "VINIT" �� "BDSOBL".
������ ������� ������������� �������� ����������(�� "BDSOBL").
���� ��������� ������� M, �������� ����� ������� M ��������������.
����� �  V �� ����� ���� �������� ��������� ��������.
"RECIPD" �������� 1/(����� ����������� �� ����).
"INSOBL" ������ ��������� ��� ������ ������� �������, �� �� ���������� �� �� ����� ���������.
��������� ������� ������� �� "GOSOBL".
"LASTQ" ������ ��������� ���������� ���������������� �������.
"TAUS" ��� ����������� "�������������" ��������.
��� � ���������� � BRATLEY/FOX, ��� ����� ��� N = 2**K ,��� K �� ����������� (TAUS+S-1)
��� ����������� K � ������� TAUS ��� ���������� �����������.
������� ������ :
�� ��������� ������������ : DIMEN, ATMOST
��� ����� DATA "BDSOBL" : POLY, VINIT
�������� ������ :
��� ��������� ������������ : FLAG, TAUS
��� "GOSOBL" VIA /SOBOL/ : V, S, MAXCOL, COUNT, LASTQ, RECIPD
*/
int INSOBL(int DIMEN, int ATMOST)
{
    int const MAXBIT = 30;
    int L, M, NEWV;
    int TAU[12] = {0, 1, 3, 5, 8, 11, 15, 19, 23, 27, 31, 35};
    char INCLUD[MAXDEG];
    //�������� ����������:
    if (!((DIMEN >= 1) && (DIMEN <= MAXDIM) && (ATMOST > 0) && (ATMOST < pow(2,MAXBIT))))
    {
        cout << "error"<<endl;
        return 1;
    }
    (DIMEN < MAXDEG)? TAUS = TAU[DIMEN-1] : TAUS = -1;

    int I = ATMOST;
    MAXCOL = int(round(log2(I))); //����� ����� �������� � ATMOST
   /* while (I)
    {
        MAXCOL++;
        I >>= 1;
    }
*/
    //������������� ������ ������ ������� V:
    for (int I = 0; I < MAXCOL; ++I)
        V[0][I] = 1;

    //������������� ��������� ����� ������� V:
    for (int I = 1; I <= DIMEN; ++I)
    {
        //���������� �������� �������� I(�������� ����������� "BDSOBL")
        //����� ������� ���������� I �� ��������� �����������.

        int J = POLY[I-1];
        M = int(round(log2(J)))-1;
        //J >>= 1;
      /*  while (J)
        {
            J >>= 1;
            ++M;
        }*/
        //�� ������������ ���� ��� �� ����������
        J = POLY[I-1];
        for (int K = M - 1; K >= 0; --K)
        {
            INCLUD[K] = J % 2;
            J >>= 1;
        }

        //������� �������� ���� I ������� �� VINIT
        for (int J = 0; J < M; ++J)
            V[I][J] = VINIT[I][J];

        //��������� ��������� �������� ���� I
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

    //����������� ������ V �� ��������������� ������� ������
    L = 1;
    for (int J = MAXCOL - 2; J >= 0; --J)
    {
        L <<= 1;
        for( I = 0; I < DIMEN; ++I)
            V[I][J] = V[I][J]*L;
    }
    RECIPD = 1.0/(L<<1); //RECIPD IS 1/(����� ����������� ��������� � V)
    COUNT = 1;
    //������� ������ ������ � ��� �������� ��� "GOSOBL"
    for(I = 0; I < DIMEN; ++I)
        LASTQ[I] = 0;

    return 0;
}


/*
��� ������������ "GOSOBL" ������� �������������� ������� � ������ �������
������������ ������ �������� "INSOBL" �� ������ "GOSOBL".  ����� ������ "INSOBL",
������� ������:
V        ������� ������������ �����
S        �����������
MAXCOL   �������� ������� V, ������� ����� ������������
COUNT    ���������� ����� ������
LASTQ    ��������� ��� ���������� ���������������� �������
RECIPD   (1/�����������) ��� ���� ���������� ������� ������� ������ �� ���� � COUNT
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

    //������ ����� ����������� QUASI
    for(I = 0; I < DIMEN; I++)
        QUASI[I] = (LASTQ[I]^=V[I][L])*RECIPD;
    return 0;
}





int cmp(const void *a, const void *b){
    struct point *x = (struct point *) a;
    struct point *y = (struct point *) b;
    return 2*((x->value - y->value)>0)-1;
}


// ������� ���������:
// n - ����������� ������������
// a,b - ������ � ����� ����� n ������� ����
// m - ���������� �����, ������� ����� �������� �� ������
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
    int k = pow(r[n-2],n); // ���������� ����� � �����

    f();

   // ofstream fout("out_all.txt"); // ����, � ������� �������� k ����� �����
   // ofstream fout1("out_min.txt"); // ����, � ������� �������� ��� m �����, ��� ������� ����������
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
        GOSOBL(j,n); // ��������� ���� ����� �����, �� ���������� ������������ � ������ QUASI
        for (int i = 0; i < n; ++i)
        {
            ANSWER[j+1][i] = QUASI[i] = a + (b - a) * QUASI[i]; // ���������� � ������� ���������� �����
            //fout << QUASI[i] << " "; // ���������� � ���� ���������� �����
        }
       // fout<<endl;
        values[j+1].value = (*func)(QUASI); //���������� � ������ �������� ������� � ������ �����
        values[j+1].point = j+1; //���������� � ������ ����� ������ �����
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
    qsort(values, k, sizeof(struct point), cmp); // ��������� ������ �������� ������� � ���������� ����� ����� ����������� ������ ��������
    for (int p = 0; p < m; ++p)
        for (int t = 0; t < n; ++t)
            ANSWER1[p][t] = ANSWER[values[p].point][t]; // ���������� � ���� ���������� �����, � ������� ������� ��������� ����������� ��������

    //fout1.close();
    delete values;
    return ANSWER1;
}




