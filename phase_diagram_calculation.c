#include<stdio.h>
#include<math.h>
#include"stack.h"           // verified... it's working properly
#define step_size 0.01      // defines step size on composition axis for Graham Scan Algorithm
#define R 0.008314
#define n 2                 // n is the no. of phases
#define dT 10                // temp jump for phase diagram calculation
 
int g_cnt=0;                // global counter for TX array
double TX[2000][6];         // to store the change in phase to draw Phase Diagram
double TA=800, TB=1200;     // TA and TB are the transformation temperatures of pure A and B respectively
double SA=0.01, SB=0.01;    // SA and SB  are entropy change from α to β for A and B respectively
double W1,W2;               // W1 and W2 are interchange energy parameters for α to β phases respectively
int p=1/step_size -1;       // p is the no. of points taken between X=0 to X=1
int T1,T2;                  // T1-T2 is the temp range in which we are calculating Phase Diagram with interval of dT
int q;                      // q is the no. of temp steps
 
// FUNCTIONS START HERE
double G_alpha(double x,double T);   // calculate G_alpha at any given Xb and T  //NO NEED TO CALL THIS, IT'S BEING USED AS NESTING
double G_beta(double x,double T);    // calculate G_beta at any given Xb and T   //NO NEED TO CALL THIS, IT'S BEING USED AS NESTING
double *G_phases(double x,double T); // returning the G values of all phases in an array at constant X and T //NO NEED TO CALL THIS ALSO
int Min(double G[]);                 // returns index corresponding to min energy phase for example 0 for alpha and 1 for beta
double *g_min(double T);             // returns the 2D array containing phases with minimum energy in complete range of X at constant Temp T
double slope(double g1,double x1,double g2,double x2);       // returns slope of line from two points
void Graham_Scan(double T);          // will store data for convex hull in stack at a Temp T
double *PD_data();                   // returns 3D array containing data for PD construction
double *U_array(double x1,double x2,double T);               // returns a 1D array containing values for Newton-Raphson method
 
 
double G_alpha(double x,double T)    // calculate G_alpha at any given Xb and T  //NO NEED TO CALL THIS, IT'S BEING USED AS NESTING
{
    double G;
    G = (1-x)*x*W1 + R*T*((1-x)*log(1-x) + x*log(x));
    return G;
}
double G_beta(double x,double T)    // calculate G_beta at any given Xb and T   //NO NEED TO CALL THIS, IT'S BEING USED AS NESTING
{
    double G=0;
    G=(1-x)*(TA-T)*SA + x*(TB-T)*SB + (1-x)*x*W2 + R*T*((1-x)*log(1-x)+x*log(x));
    return G;
}
double *G_phases(double x,double T)  // returning the G values of all phases in an array at constant X and T //NO NEED TO CALL THIS ALSO
{
    static double G[n];
    G[0]=G_alpha(x,T); G[1]=G_beta(x,T);
    printf("G1 is %0.02lf and G2 is %0.02lf\n",G[0],G[1]);
    return G;
}
int Min(double G[])     //returns index corresponding to min energy phase for example 0 for alpha and 1 for beta
{
    int i=0;
    int min=0;
    for(i=0;i<n;i++)
    {
        if(G[i]<=G[min])
            min=i;
    }
    return min;
}
double *g_min(double T)   // returns the 2D array containing phases with minimum energy in complete range of X at constant Temp T
{
    int i=0,ind;
    double *g;
    static double G[100][4];    // here G[0] is having composition, G[2] having index of min energy phase at G[0], G[1] has that min energy value and G[3] has Temp
    for(i=0;i<p;i++)
    {
        g=G_phases((i+1)*step_size,T);
        ind=Min(g);
        G[i][0]=(i+1)*step_size;
        G[i][2]=ind;
        G[i][1]=g[ind];
        G[i][3]=T;
        printf("At Temp=%0.02lf and comp. X=%0.2lf min G =%0.02lf and min phase is %0.0lf\n",G[i][3],G[i][0],G[i][1],G[i][2]);
    }
    printf("\n");
    return G;
}
double slope(double g1,double x1,double g2,double x2)       //return slope of line from two points
{
    double m=(g1-g2)/(x1-x2);
    return m;
}
 
void Graham_Scan(double T)              // will store data for convex hull in stack at a Temp T
{
    double *g=g_min(T);
    double G[100][4];
    double prev=-999999999999,curr;
    double *g1,*g2;                     // g1 will receive data from Peek() function
    int i=1;
    for(i=1;i<=p;i++)
    {
        G[i][0]=*g;
        g+=1;
        G[i][1]=*g;
        g+=1;
        G[i][2]=*g;
        g+=1;
        G[i][3]=*g;
        g+=1;
    }
 
    Push(G[1]);
    for(i=2;i<=p;i++)
    {
        g1=Peek();
        curr=slope(g1[1],g1[0],G[i][1],G[i][0]);
        if(curr>prev)
        {
            Push(G[i]);
            prev=curr;
        }
        else if(curr<=prev)
        {
            Pop();
            while(curr<prev)
            {
                g1=Peek();
                curr=slope(g1[1],g1[0],G[i][1],G[i][0]);
                Pop();
                g2=Peek();
                prev=slope(g2[1],g2[0],g1[1],g1[0]);
            }
            g1[2]=-1;
            G[i][2]=-1;
            Push(g1);
            Push(G[i]);
            prev=curr;
        }
    }
}
 
double *PD_data()
{
    int max,min,i=0,j,count,k,cnt=2;
    min=T1;max=T2;
    double *g;
    double G[500][100][4];
    if(T1>T2)
    {
        max=T1;min=T2;
    }
    q =(max-min)/dT;
    for(k=0;k<=q;k++)
    {
        Graham_Scan(min+dT*k);
        cnt=2;
        TX[g_cnt][0]=0;
        TX[g_cnt][1]=min+dT*k;
 
        i=r.top;
        count=i;
        printf("******* count= %d\n",i);
        while(i!=0)
        {
            g=Peek();
            for(j=0;j<4;j++)
                G[k][i][j]=g[j];
            Pop();
 
            if(G[k][i][2]==-1)
            {
                TX[g_cnt][0]++;
                TX[g_cnt][cnt++]=G[k][i][0];
            }
            i--;
        }
        printf("******* count= %d\n",i);
        g_cnt++;
    }
    return G;
}
 
// MAIN CODE START HERE
int main()
{
    printf("Enter the value of W1=");
    scanf("%lf",&W1);
    printf("Enter the value of W2=");
    scanf("%lf",&W2);
    printf("Enter the start temp T1 for PD calculation=");
    scanf("%d",&T1);
    printf("Enter the final temp T2 for PD calculation=");
    scanf("%d",&T2);
 
    int i=0,j=0,count=0;      // p contains the no. of points taken on composition axis
    double *g;
    double G[p][4];
    char output[]="PD_data.txt";
    FILE*fp;
    fp=fopen(output,"w");
    printf("\nAlgorithm starts here\n");
    g=PD_data();
    printf("\n");
    for(i=0;i<g_cnt;i++)
    {
        for(j=1;j<TX[i][0]+2;j++)
        {
            printf("%8.3lf  ",TX[i][j]);
            fprintf(fp,"%8.3lf\t",TX[i][j]);
        }
        printf("\n");
        fprintf(fp,"\n");
    }
    fclose(fp);
    return 0;
}
