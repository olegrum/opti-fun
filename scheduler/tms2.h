
#include <algorithm>
#include <iomanip>
using namespace std;

// COMM
const int maxdim = 20;
const int maxfig = 20;
const int nbits = 31;

int c[maxdim][maxfig][maxfig/*!*/]; // /*!*/ - ���������� ��� � ��, �.�. C(MAXDIM, MAXFIG, 0:MAXFIG-1)
int count[maxfig/*!*/], d[maxdim][maxfig];
int nextq[maxdim], qpow[maxfig];
int dimen, nfigs;
double recip;

// COMM2
int cj[maxdim][nbits/*!*/];
int count2;
int dimen2;
int nextq2[maxdim];

// FIELD
const int maxq = 50;
const int maxdeg = 50, deg = 0;
int p, q, add[maxq/*!*/][maxq/*!*/], mul[maxq/*!*/][maxq/*!*/], sub[maxq/*!*/][maxq/*!*/];


vector < vector <double> > vResult;
vector <double> vector_;
int ResultCounter = 0;

struct point
{
    int point;
    double value;
};

void SETFLD(int qin)
{
    cout << "SETFLD" << endl;
    int i,j;
    if (qin <= 1 || qin > maxq) cout << "SETFLD: Bad value of Q";
    q = qin;
    p = q;
    if (p == 0) cout << "SETFLD: There is no field of order" << q;

    for(i = 0; i < q; ++i)
        for(j = 0; j < q; ++j)
        {
            add[i][j] = (i + j) % p;
            mul[i][j] = (i * j) % p;
        }

    for (i = 0; i < q; ++i)
        for (j = 0; j < q; ++j) sub[add[i][j]][i] = j;
}

void PLYMUL(int *pa, int *pb, int *pc)
{
    int i,j,dega,degb,degc,term;
    int pt[maxdeg+2];
    dega = pa[deg];
    degb = pb[deg];
    if (dega == -1 || degb == -1) degc = -1;
    else degc = dega + degb;
    if (degc > maxdeg)
        cout << "PLYMUL: Degree of product exceeds MAXDEG" << endl;
    for (i = 0; i <= degc; ++i)
    {
        term = 0;
        for(j = max(0,i-dega); j <= min(degb,i); ++j)
            term = add[term][mul[pa[i-j+1]][pb[j+1]]];
        pt[i+1] = term;
    }
    pc[deg] = degc;
    for(i = 1; i <= degc + 1; ++i) pc[i] = pt[i];
    for(i = degc + 2; i < maxdeg + 2; ++i) pc[i] = 0;
}

void CALCV(int *px,int *b,int *v,int maxv)
{
    int h[maxdeg+2];
    int bigm = 0, m = 0, kj, term;
    int arbit = 1,nonzer = 1;

    for(int i = 0; i < b[deg]+2; ++i)
        h[i] = b[i];
    bigm = h[deg];
    PLYMUL(px,b,b);

    m = b[deg];

    kj = bigm;
    for(int r = 0; r < kj; ++r)
        v[r] = 0;
    v[kj] = 1;

    if(kj < bigm)
    {
        // weak branch
        term = sub[0][h[kj+1]];
        for(int r = kj - 1; r < bigm - 1; ++r)
        {
            v[r] = arbit;
            term = sub[term][ mul[h[r+1]][v[r]] ];
        }
        v[bigm] = add[nonzer][term];
        for(int r = bigm + 1; r < m - 1; ++r)
            v[r] = arbit;
    }
    else
        for(int r = kj + 1; r <= m - 1; ++r)
            v[r] = arbit;
    for(int r = 0; r <= maxv - m; ++r)
    {
        term = 0;
        for(int i = 0; i <= m - 1; ++i)
            term = sub[term][ mul[b[i+1]][v[r+i]] ];
        v[r+m] = term;
    }
}

void CALCC2()
{
    int maxe = 5, maxv = nbits + maxe;
    int px[maxdeg+2], b[maxdeg+2];
    int v[maxv+1], ci[nbits][nbits];
    int e,i,j,r,u,term;
    int irred[maxdim][maxe+2];
    {
        irred[0][0] = 1;
        irred[0][1] = 0;
        irred[0][2] = 1;

        irred[1][0] = 1;
        irred[1][1] = 1;
        irred[1][2] = 1;

        irred[2][0] = 2;
        irred[2][1] = 1;
        irred[2][2] = 1;
        irred[2][3] = 1;

        irred[3][0] = 3;
        irred[3][1] = 1;
        irred[3][2] = 1;
        irred[3][3] = 0;
        irred[3][4] = 1;

        irred[4][0] = 3;
        irred[4][1] = 1;
        irred[4][2] = 0;
        irred[4][3] = 1;
        irred[4][4] = 1;

        irred[5][0] = 4;
        irred[5][1] = 1;
        irred[5][2] = 1;
        irred[5][3] = 0;
        irred[5][4] = 0;
        irred[5][5] = 1;

        irred[6][0] = 4;
        irred[6][1] = 1;
        irred[6][2] = 0;
        irred[6][3] = 0;
        irred[6][4] = 1;
        irred[6][5] = 1;

        irred[7][0] = 4;
        irred[7][1] = 1;
        irred[7][2] = 1;
        irred[7][3] = 1;
        irred[7][4] = 1;
        irred[7][5] = 1;

        irred[8][0] = 5;
        irred[8][1] = 1;
        irred[8][2] = 0;
        irred[8][3] = 1;
        irred[8][4] = 0;
        irred[8][5] = 0;
        irred[8][6] = 1;

        irred[9][0] = 5;
        irred[9][1] = 1;
        irred[9][2] = 0;
        irred[9][3] = 0;
        irred[9][4] = 1;
        irred[9][5] = 0;
        irred[9][6] = 1;

        irred[10][0] = 5;
        irred[10][1] = 1;
        irred[10][2] = 1;
        irred[10][3] = 1;
        irred[10][4] = 1;
        irred[10][5] = 0;
        irred[10][6] = 1;

        irred[11][0] = 5;
        irred[11][1] = 1;
        irred[11][2] = 1;
        irred[11][3] = 1;
        irred[11][4] = 0;
        irred[11][5] = 1;
        irred[11][6] = 1;

        irred[12][0] = 5;
        irred[12][1] = 1;
        irred[12][2] = 1;
        irred[12][3] = 0;
        irred[12][4] = 1;
        irred[12][5] = 1;
        irred[12][6] = 1;

        irred[13][0] = 5;
        irred[13][1] = 1;
        irred[13][2] = 0;
        irred[13][3] = 1;
        irred[13][4] = 1;
        irred[13][5] = 1;
        irred[13][6] = 1;

        irred[14][0] = 6;
        irred[14][1] = 1;
        irred[14][2] = 1;
        irred[14][3] = 0;
        irred[14][4] = 0;
        irred[14][5] = 0;
        irred[14][6] = 0;
        irred[14][7] = 1;

        irred[15][0] = 6;
        irred[15][1] = 1;
        irred[15][2] = 0;
        irred[15][3] = 0;
        irred[15][4] = 1;
        irred[15][5] = 0;
        irred[15][6] = 0;
        irred[15][7] = 1;

        irred[16][0] = 6;
        irred[16][1] = 1;
        irred[16][2] = 1;
        irred[16][3] = 1;
        irred[16][4] = 0;
        irred[16][5] = 1;
        irred[16][6] = 0;
        irred[16][7] = 1;

        irred[17][0] = 6;
        irred[17][1] = 1;
        irred[17][2] = 1;
        irred[17][3] = 0;
        irred[17][4] = 1;
        irred[17][5] = 1;
        irred[17][6] = 0;
        irred[17][7] = 1;

        irred[18][0] = 6;
        irred[18][1] = 1;
        irred[18][2] = 0;
        irred[18][3] = 0;
        irred[18][4] = 0;
        irred[18][5] = 0;
        irred[18][6] = 1;
        irred[18][7] = 1;

        irred[19][0] = 6;
        irred[19][1] = 1;
        irred[19][2] = 1;
        irred[19][3] = 1;
        irred[19][4] = 0;
        irred[19][5] = 0;
        irred[19][6] = 1;
        irred[19][7] = 1;
    }
    SETFLD(2);

    for(i = 0; i < dimen; ++i)
    {
        e = irred[i][deg];
        b[deg] = 0;
        b[1] = 1;
        u = 0;
        for(j = 0; j < e + 2; ++j)
            px[j] = irred[i][j];

        for(j = 0; j <= nbits - 1; ++j)
        {
            if (u == 0) CALCV(px, b, v, maxv);
            for(r = 0; r <= nbits - 1; ++r)
                ci[j][r] = v[r+u];
            ++u;
            if(u == e) u = 0;
        }
        for(r = 0; r <= nbits - 1; ++r)
        {
            term = 0;
            for(j = 1; j <= nbits; ++j)
                term = (term << 1) + ci[j-1][r];
            cj[i][r] = term;
        }
    }
}

void INLO2(int dim, int skip)
{
    int r, gray;
    dimen2 = dim;
    if(dimen2 <= 0 || dimen2 > maxdim)
    {
        cout<<"INLO2 : Bad dimension";
        return;
    }

    CALCC2();
    gray = skip xor skip >> 1;
    for(int i = 1; i <= dimen2; ++i) nextq2[i-1] = 0;
    r = 0;
    while(gray != 0)
    {
        if(gray % 2 != 0)
            for(int i = 1; i <= dimen2; ++i)
                nextq2[i-1] = nextq2[i-1] xor cj[i-1][r];
        gray >>= 1; //����� �� 2
        r++;
    }
    count2 = skip;
}

void GOLO2(double *quasi, vector< vector<double> >& bounds)
{
    int r;
    recip = pow(2, -nbits);
    for(int i = 0; i < dimen2; ++i)
    {
        quasi[i] = nextq2[i] * recip;
        vResult[ResultCounter][i] = bounds[i][0] + (bounds[i][1]-bounds[i][0])*quasi[i];
    }
    ++ResultCounter;

    r = 0;
    int i = count2;
    while(i % 2 != 0)
    {
        r++;
        i >>= 1;
    }
    if(r >= nbits)
    {
        cout<<"GOLO2 : Too many calls";
        return;
    }
    for(int I = 0; I < dimen2; ++I)
        nextq2[I] = nextq2[I] xor cj [I][r];
    ++count2;
}

double func(vector<double>& p)
{
    return p[0]*p[0]+p[1]*p[1];
}

int cmp(const void *a, const void *b)
{
    struct point *x = (struct point *) a;
    struct point *y = (struct point *) b;
    return 2*((x->value - y->value)>0)-1;
}

void GENIN2(int dimen_, int seqlen_, int m_, double (*fun)(vector<double>&), vector< vector<double> >& bounds) // PROGRAM
{
    dimen=dimen_;
    int seqlen=seqlen_;
    int m=m_;

    // �������� ������ ��� �����. ���������� �����: seqlen; ��������� �� �����: dimen
    vResult.resize(seqlen);
    vector_.resize(seqlen);
    for(int i=  0; i < seqlen; ++i)
        vResult[i].resize(dimen);

    int skip = 0;

    INLO2(dimen, skip);
    cout << "GENIN2 :  Initialization complete" << endl;

    double quasi[maxdim];

    for(int i = 1; i <= seqlen; ++i)
        GOLO2(quasi, bounds);

    cout << "GENIN2:  iteration ends" << endl;


    //int k; //���������� ����� � �����
    //ofstream fout("out.txt");
    struct point * values = (struct point*) calloc(seqlen,sizeof(struct point));

    for (int i = 0; i < seqlen; ++i)
    {
        for (int j = 0; j < dimen; ++j)
        {
            vector_[j]=vResult[i][j];
           // fout<< fixed << setprecision(9)<<vResult[i][j]<<" ";
        }
       // fout<<endl;
        values[i].value = fun(vector_); //���������� � ������ �������� ������� � ������ �����
        values[i].point = i; //���������� � ������ ����� ������ �����
    }
    //fout.close();

    qsort(values, seqlen, sizeof(struct point), cmp); // ��������� ������ �������� ������� � ���������� ����� ����� ����������� ������ ��������

   // ofstream fout1("out_min.txt"); // ����, � ������� �������� ��� m �����, ��� ������� ����������

 /*   for (int p = 0; p < m; ++p)
    {
        for (int t = 0; t < dimen; ++t)
            //fout1<< vResult[values[p].point][t]<< " "; // ���������� � ���� ���������� �����, � ������� ������� ��������� ����������� ��������
       // fout1 <<" "<< values[p].value<< endl;
    }
*/
   // fout1.close();
    //cout << "Answer in out_min.txt" << endl;
    free (values);
}

/*
int main()
{
    int m,d,s;
    vector< vector<double> > bounds;
    do
    {
        cout << "Enter dimension: ";
        cin >> d;
        if( d > maxdim ) cout << "Dimension may not exceed " << maxdim << endl;
    }
    while(d > maxdim);

    bounds.resize(d);
    cout << "Enter bound for every dimension: ";
    for(int i=0; i < d; ++i)
    {
        bounds[i].resize(2);
        cout<<"Enter 2 bounds for "<<i+1<<" dimension: ";
        cin>>bounds[i][0]>>bounds[i][1];
    }

    cout << "Choose sequence length from list below: " << endl;
    int pow_;
    pow_ = 1 << 10;
    cout << "2^10 = " << pow_ << endl;
    pow_ = 1 << 15;
    cout << "2^15 = " << pow_ << endl;
    pow_ = 1 << 20;
    cout << "2^20 = " << pow_ << endl;
    do
    {
        cout << "Sequence length: ";
        cin >> s;
        if( s < 0 ) cout << "Length must be strictly positive" << endl;
    }while(s < 0);
    cout << "How many dots you want: ";
    cin >> m;

    double (*f_pointer)(vector<double>&) = &func;

    GENIN2(d, s, m, f_pointer, bounds);

    return 0;
}*/
