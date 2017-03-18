#include <pthread.h>
#include <iostream>
#include <stack>
#include <fstream>
#include <ctime>
#include <time.h>
#include <limits>
#include <new>
#include "tms.h"
#include "powell.h"

#define ERROR_CREATE_THREAD   -11
#define ERROR_JOIN_THREAD     -12
#define SUCCESS                 0
#define ERROR                  -1
#define SLEEP                   10
#define WAIT_TIME               1

#if _WIN64 || _WIN32                            // CHECKING NUMBER OF CORES
#include <windows.h>
int NUM(){
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
}
#define sleep(x) Sleep(x)
#else
#include <unistd.h>
int NUM(){
    return sysconf(_SC_NPROCESSORS_ONLN);
}
#endif
#define NUM NUM()                               // END_CHECKING

using namespace std;


class globe{
    public:
    double** matr;
    int dim;                                  // DIMENSION OF SPACE
    pthread_mutex_t mut;                     // MUTEX
    stack <int> mystack1;                           // STACK FOR FREE CORES
    pthread_cond_t stack_cond;               // CONDITIONAL VARIABLE , WHICH IS EXPECTING THREADS
    double accur;                       // OPTIONS OF FUNCTION AND ACCURACY
};
static globe glb;

typedef class Mine_str{                        // STRUCTER FOR EACH POINT FROM TMS-NETS
    public:
    double *xn;                                  // COORDINATES OF POINT
    int var;                                    // THREAD NUMBER , WHICH IS EXECUTING METHOD IN THIS POINT
    double inf;                                  // MINIMUM IN THIS POINT
    Mine_str()
    {
        try{
            xn=new double[glb.dim];
        }
        catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << '\n';
        }
        var=0;
        inf=std::numeric_limits<double>::infinity();
    }
    ~Mine_str()
    {
        delete [] xn;
    }
} mine_t;

double f0(double* x) {                           // TARGET FUNCTION
    double sum;
    /*int i,j;

    sum=glb.matr[0][0];

    for(i=1; i<=glb.dim; i++)
        for(j=1; j<=glb.dim; j++)
        {
            sum+=glb.matr[i][j]*x[i-1]*x[j-1];
        }*/
    //sum=x[2]*x[2]+x[0]*x[0]+(x[0]-1)*(x[0]-1)+(1-x[1])*(1-x[1]);
    //sum=2*x[0]*x[0]+x[1]*x[1]+x[0]*x[1];
    sum=100*(x[1]-x[0]*x[0]*x[0])*(x[1]-x[0]*x[0]*x[0])+(1-x[0])*(1-x[0]);
    return sum;
}

void* DFP(void* args) {                                                     // METHOD DFP

     if (pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL) != 0) {        // MAIN CAN CANCEL THIS THREAD
      printf("Thread pthread_setcancelstate failed (DFP)");
      exit(EXIT_FAILURE);
     }
     if (pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL) != 0) {   // IN ANY MOMENT
      printf("Thread pthread_setcanceltype failed (DFP)");
      exit(EXIT_FAILURE);
     }
   /* srand(time(0));
    double rnd=(rand()%SLEEP);
    sleep(rnd*2);*/

    mine_t* arg=(mine_t*) args;
    pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);

    pthread_mutex_lock(&glb.mut);
    glb.mystack1.push(arg->var);                                                // RETURN NUMBER OF CORE TO STACK
    pthread_mutex_unlock(&glb.mut);

    if ( pthread_cond_signal(&glb.stack_cond) != 0){                            // SENDING SIGNAL TO MAIN
           printf("Pthread_cond_signal failed (DFP)");
           exit(EXIT_FAILURE);
    }
    pthread_exit (NULL);                                                    // DESTROING THREAD
    return NULL;
}                                                                           // END_DFP

void* Paul(void* args) {                                                    // METHDO POWELL ( WORKING AS DFP )

    if (pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL) != 0) {
      printf("Thread pthread_setcancelstate failed (POWELL)");
      exit(EXIT_FAILURE);
     }
     if (pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL) != 0) {
      printf("Thread pthread_setcanceltype failed (POWELL)");
      exit(EXIT_FAILURE);
     }

   /* srand(time(0));
    double rnd=(rand()%SLEEP);
    sleep(rnd*2);*/

    mine_t* arg=(mine_t*) args;

    //arg->inf=powell(arg->xn,glb.dim+1,glb.accur,&f0);
    pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);

    pthread_mutex_lock(&glb.mut);
    glb.mystack1.push(arg->var);
    pthread_mutex_unlock(&glb.mut);

    if ( pthread_cond_signal(&glb.stack_cond) != 0){
           printf("Pthread_cond_signal failed (POWELL)");
           exit(EXIT_FAILURE);
    }
    pthread_exit (NULL);
    return NULL;
}                                                                           // END POWELL


void create_thr(int* num_fun, mine_t* fun, int* cnt, pthread_t thr[]){      // CREATING THREAD
    int status;

    switch (*num_fun) {                                                     // CHECKING WHOSE TURN

    case 0:{
        fun->var=glb.mystack1.top(); glb.mystack1.pop();                            // TAKE NUMBER OF CORE FROM STACK
        status = pthread_create(&thr[fun->var], NULL, DFP, (void*)fun);     // CREATING THREAD
        if (status != 0) {
                printf("main error: can't create thread, status = %d\n", status);
                exit(ERROR_CREATE_THREAD);
        }
        *cnt+=1;                                                            // iD = iD+1
        *num_fun=1;                                                         // THE NEXT IS POWELL
    }
    break;
    case 1:{
        fun->var=glb.mystack1.top();
        glb.mystack1.pop();
        status = pthread_create(&thr[fun->var], NULL, Paul, (void*)fun);
        if (status != 0) {
                printf("main error: can't create thread, status = %d\n", status);
                exit(ERROR_CREATE_THREAD);
        }
        *cnt+=1;
        *num_fun=0;
    }
    break;
    default: exit(ERROR);

    }
}                                                                           // END CREATING THREAD

int main()  {

    /***********DATA LOADING********************/
    int i,j,k;                                                              // k - AMOUNT OF POINTS , WHICH TMS IS SEARCHING
    double a,b;                                                              // INTERVAL [a,b]

    ifstream fin("args.txt");
    fin >> glb.dim >> glb.accur >> a >> b>> k;

    try{
        glb.matr = new double * [glb.dim+1];
        for (i=0; i<=glb.dim; i++)
            glb.matr[i] = new double[glb.dim+1];
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }

    for(i=0;i<=glb.dim;i++) {
            for(j=0;j<=glb.dim;j++){
                fin >> glb.matr[i][j] ;
            }
    }
    fin.close();


    /*************************************************************************************************************/


    /***********DATA PROCESSING*****************/
    double** points=tms(glb.dim, a, b, k, f0);                      //TMS

    // TMS IS RIGHT HERE

    int status,iP=0,iD=0,num_fun=0,status_addr;;                           // iP,iD - INDICES OF POINTS FOR DFP AND POWELL
    mine_t point_D[k], point_P[k];                                          // ARRAY OF POINTS FOR EACH METHOD
    pthread_t thr[NUM-1];                                                   // NUM-1 THREADS FOR METHODS
    struct timespec stime;

    for (i=1; i<=NUM-1; i++) glb.mystack1.push(i);
    for (i=0; i<=k-1; i++) {                                                // CREATING 2 ARRAYS [k,dim] FOR DFP AND POWELL METHODS
        for (j=0; j<=glb.dim-1; j++){
            point_P[i].xn[j]=points[i][j];
            point_D[i].xn[j]=points[i][j];
        }
    }

    for (i=0; i<k; i++)
        delete [] points[i];
    delete [] points;

    if ( pthread_mutex_init(&glb.mut, NULL) != 0 ){                             // MUTEX INITIALIZATION
        printf(" pthread_mutex_init failed ");
        exit(EXIT_FAILURE);
    }
    if ( pthread_cond_init(&glb.stack_cond, NULL) != 0 ){                       // CONDITIONAL VARIABLE INITIALIZATION
        printf(" pthread_cond_init failed ");
        exit(EXIT_FAILURE);
    }

    while(iP<k || iD<k || (int)glb.mystack1.size()!=NUM-1){                                            // RUNNING OVER ALL POINTS
            pthread_mutex_lock(&glb.mut);                                       // MUTEX LOCK
            if(glb.mystack1.size()>=1){                                         // CHECKING STACK
                pthread_mutex_unlock(&glb.mut);                                 // UNLOCK MUTEX
                switch (num_fun) {                                         // IF 0 THEN DFP GOES ; ELSE POWELL GOES
                case 0: {
                    create_thr(&num_fun , &point_D[iD] ,  &iD , thr);
                }
                break;
                case 1: {
                    create_thr(&num_fun , &point_P[iP] ,  &iP , thr);
                }
                break;
                default: break;
                }
            }else{
                stime.tv_sec = time(NULL) + WAIT_TIME;                                          // ALLOWED TIME FOR WAITING
                if ( pthread_cond_timedwait (&glb.stack_cond,&glb.mut,&stime) == ETIMEDOUT ) {          // MAIN IS WAITING COMPLETED THREADS
                    cout<<"Hello"<<endl;
                    for ( i=1; i<=NUM-1; i++) {
                        pthread_cancel(thr[i]);                                                 // ELSE CANCEL ALL THREADS
                    }
                    for (i = 1; i <= NUM-1; i++) {
                        status = pthread_join(thr[i], (void**)&status_addr);
                        if (status != SUCCESS) {
                            printf("main error: can't join thread, status = %d\n", status);
                            exit(ERROR_JOIN_THREAD);
                        }
                        glb.mystack1.push(i);
                    }
                }
                pthread_mutex_unlock(&glb.mut);                                                     // MUTEX UNLOCK
            }
            num_fun=(iD==k ? (iP==k ? 2 : 1) : (iP==k ? 0 : num_fun) );
    }


    /***********************************************************************************************/



    /***************SAVING DATA***********************/

    for (i=0; i<=glb.dim; i++)
        delete [] glb.matr[i];
    delete [] glb.matr;

    if ( pthread_cond_destroy(&glb.stack_cond) != 0 ){
        printf(" pthread_cond_destroy failed");
        exit(EXIT_FAILURE);
    }
    if ( pthread_mutex_destroy (&glb.mut) !=0 ){
        printf(" pthread_mutex_destroy failed");
        exit(EXIT_FAILURE);
    }
    /***********************************************************************************************/

    return 0;
    }









