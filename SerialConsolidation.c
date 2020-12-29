//Consolidation code in serial to run on Aimos

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
//#include "petscmat.h" 

int nnodes;
int nel;
int np = 4;
int *offset;
int *nodeParts;
float *coordinatesX = NULL;
float *coordinatesY = NULL;
float *coordinatesZ = NULL;

float *coordinatesXPow = NULL;
float *coordinatesYPow = NULL;
float *coordinatesZPow = NULL;

int *elementsW = NULL;
int *elementsX = NULL;
int *elementsY = NULL;
int *elementsZ = NULL;

int *elementsWPow = NULL;
int *elementsXPow = NULL;
int *elementsYPow = NULL;
int *elementsZPow = NULL;

float *psi = NULL;
int Tol = 1;
float powder_thick = 30; //microns
float porosity = 0.55;

float N[4];
float dN[12];
float jac;
float ke[16];
float fe[4];
float *a_bar = NULL;

void printLine(int line)
{
  char fileName[100];
  snprintf(fileName,100,"errorAtLineCuda.txt");
  FILE *fp;
  fp = fopen(fileName,"w+");
  fprintf(fp,"Line is %d.\n",line);
  fclose(fp);
}

static inline void gol_swap(char **pA,char **pB)
{
  // You write this function - it should swap the pointers of pA and pB.

  //Create temporary variable "temp" for the use of swapping the pointers pA and pB for the
  //old and new worlds
  char *temp;
    temp = *pA;
    *pA = *pB;
    *pB = temp;
}

void num_ElementsNodes(char BaseName1[80], int myrank)
{
  char BaseName[80];
  strcpy(BaseName,BaseName1);
  
  int nnodesl;
  int ncellsl;
  char fname [100];
  snprintf(fname,100,"%d.vtu",myrank);
  //printf("%s\n",fname);
  strcat(BaseName, fname);
  //printf("%s\n",BaseName);

  FILE *fp;
  fp = fopen(BaseName,"r");
  char tline1[100];
  char tline2[100];
  
  fgets(tline1, 100, fp);
  //printf("%s\n",tline1);
  fgets(tline2, 100, fp);
  char tline[100];
  fgets(tline, 100, fp);
  //printf("%s\n",tline);
  
  sscanf(tline,"<Piece NumberOfPoints=\"%d\"",&nnodesl);
  printf("Number of nodes is %d\n",nnodesl);

  char str1 [200] = "<Piece NumberOfPoints=\"";
  char str2 [30];
  snprintf(str2,30,"%d",nnodesl);
  char str3 [3] = "\"";
  strcat(str1,str2);
  strcat(str1,str3);
  char str4 [100] = " NumberOfCells=\"%d\">";
  strcat(str1,str4);
  //printf("%s\n",str1);
  sscanf(tline,str1, &ncellsl);
  printf("Number of cells is %d\n",ncellsl);

  nnodes = nnodesl;
  nel = ncellsl;

  fclose(fp);
}

void offsetCalc(char BaseName1[80], int numranks, int myrank)
{
  int i;
  nodeParts = calloc(np, sizeof(int));
  offset = calloc(np, sizeof(int));

  printLine(__LINE__);
  for(i=0;i<numranks;i++)
  {
    char BaseName[80];
    strcpy(BaseName,BaseName1);
    
    int nnodesl;
    char fname [100];
    snprintf(fname,100,"%d.vtu",i);
    //printf("%s\n",fname);
    strcat(BaseName, fname);
    //printf("%s\n",BaseName);
    printLine(__LINE__);
    FILE *fp;
    fp = fopen(BaseName,"r");
    char tline1[100];
    char tline2[100];
    printLine(__LINE__);
    fgets(tline1, 100, fp);
    printLine(__LINE__);
    //printf("%s\n",tline1);
    fgets(tline2, 100, fp);
    printLine(__LINE__);
    //printf("%s\n",tline2);
    char tline[100];
    fgets(tline, 100, fp);
    printLine(__LINE__);
    //printf("%s\n",tline);
    printLine(__LINE__);
    sscanf(tline,"<Piece NumberOfPoints=\"%d\"",&nnodesl);
    //printf("Number of nodes is %d\n",nnodesl);
    
    nodeParts[i] = nnodesl;
    
    fclose(fp);
  }
  printLine(__LINE__);
  for(i=1;i<4;i++)
    {
      offset[i] = offset[i-1] + nodeParts[i];
    }
  printf("Offset array is  %d ",offset[0]);
  printf("%d ",offset[1]);
  printf("%d ",offset[2]);
  printf("%d\n",offset[3]);

}

void readCoordinates(char BaseName1[80], int myrank, int nnodes)
{
  coordinatesX = calloc(nnodes,sizeof(float));
  coordinatesY = calloc(nnodes,sizeof(float));
  coordinatesZ = calloc(nnodes,sizeof(float));
  
  printLine(__LINE__);
  char str[100] = "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  //printf("%s\n",str);
  char BaseName[80];
  strcpy(BaseName,BaseName1);
  printLine(__LINE__);

  char fname [100];
  snprintf(fname,100,"%d.vtu",myrank);
  //printf("%s\n",fname);
  strcat(BaseName, fname);
  //printf("%s\n",baseName);
  printLine(__LINE__);

  FILE *fp;
  fp = fopen(BaseName,"r");
  printLine(__LINE__);
  
  int check = 1;
  int *point = NULL;
  char *tline1 = malloc(100);
  char *tline2 = malloc(100);
  printLine(__LINE__);
  fgets(tline1, 100, fp);
  printLine(__LINE__);
  
  int i;
  while (check == 1)
    {
      if (strcmp(str,tline1) == 0)
	{
	  point = &check;
	  *point = 0;
	  //printf("Found it");
	}
      //printf("string comparrison is %d\n",strcmp(str,tline1));
      printLine(__LINE__);
      fgets(tline2, 100, fp);
      printLine(__LINE__);
      printLine(__LINE__);
      gol_swap(&tline1,&tline2);
      printLine(__LINE__);
      //printf("%s",tline1);
    }

  for(i=0;i<nnodes;i++)
    {
      sscanf(tline1,"%f %f %f",&coordinatesX[i],&coordinatesY[i],&coordinatesZ[i]);
      fgets(tline2, 100, fp);
      gol_swap(&tline1,&tline2);
      //printf("%f %f %f\n",coordinatesX[i],coordinatesY[i],coordinatesZ[i]);

    }
   printf("Last row of coordinates are: %f %f %f\n",coordinatesX[nnodes-1],coordinatesY[nnodes-1],coordinatesZ[nnodes-1]);

   free(tline1);
   free(tline2);

   fclose(fp);
}

void readElements(char BaseName1[80], int myrank, int nel)
{
  elementsW = calloc(nel,sizeof(int));
  elementsX = calloc(nel,sizeof(int));
  elementsY = calloc(nel,sizeof(int));
  elementsZ = calloc(nel,sizeof(int));
  
  printLine(__LINE__);
  char str[100] = "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  //printf("%s",str);
  char BaseName[80];
  strcpy(BaseName,BaseName1);
  printLine(__LINE__);

  char fname [100];
  snprintf(fname,100,"%d.vtu",myrank);
  //printf("%s\n",fname);
  strcat(BaseName, fname);
  //printf("%s\n",baseName);
  printLine(__LINE__);

  FILE *fp;
  fp = fopen(BaseName,"r");
  printLine(__LINE__);

  int check = 1;
  int *point = NULL;
  char *tline1 = malloc(100);
  char *tline2 = malloc(100);
  printLine(__LINE__);
  fgets(tline1, 100, fp);
  printLine(__LINE__);
  
  int i;
  while (check == 1)
    {
      if (strcmp(str,tline1) == 0)
	{
	  point = &check;
	  *point = 0;
	  //printf("Found it");
	}
      //printf("string comparrison is %d\n",strcmp(str,tline1));
      //printLine(__LINE__);
      fgets(tline2, 100, fp);
      //printLine(__LINE__);
      gol_swap(&tline1,&tline2);
      //printLine(__LINE__);
      //printf("%s",tline1);
    }
  printLine(__LINE__);
  for(i=0;i<nel;i++)
    {
      //printLine(__LINE__);
      sscanf(tline1,"%d %d %d %d",&elementsW[i],&elementsX[i],&elementsY[i],&elementsZ[i]);
      //printLine(__LINE__);
      fgets(tline2, 100, fp);
      //printLine(__LINE__);
      gol_swap(&tline1,&tline2);
      //printLine(__LINE__);
      //elementsW[i]++;
      //elementsX[i]++;
      //elementsY[i]++;
      //elementsZ[i]++;
      //printf("%f %f %f %f\n",elementsW[i],elementsX[i],elementsY[i],elementsZ[i]);

    }
  printLine(__LINE__);
  printf("Last row of elements are: %d %d %d %d\n",elementsW[nel-1],elementsX[nel-1],elementsY[nel-1],elementsZ[nel-1]);

  free(tline1);
  free(tline2);

  fclose(fp);
}

void readPsi(char BaseName1[80], int myrank, int nel)
{
  psi = calloc(nel,sizeof(float));
  
  printLine(__LINE__);
  char str[100] = "<DataArray type=\"Float64\" Name=\"Psi1_1\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  //printf("%s",str);
  char BaseName[80];
  strcpy(BaseName,BaseName1);
  printLine(__LINE__);

  char fname [100];
  snprintf(fname,100,"%d.vtu",myrank);
  //printf("%s\n",fname);
  strcat(BaseName, fname);
  //printf("%s\n",baseName);
  printLine(__LINE__);

  FILE *fp;
  fp = fopen(BaseName,"r");
  printLine(__LINE__);
  
  int check = 1;
  int *point = NULL;
  char *tline1 = malloc(100);
  char *tline2 = malloc(100);
  printLine(__LINE__);
  fgets(tline1, 100, fp);
  //free(tline1);
  printLine(__LINE__);
  
  int i;
  while (check == 1)
    {
      if (strcmp(str,tline1) == 0)
	{
	  point = &check;
	  *point = 0;
	  //printf("Found it");
	}
      //printf("string comparrison is %d\n",strcmp(str,tline1));
      //printLine(__LINE__);
      fgets(tline2, 100, fp);
      //printLine(__LINE__);
      gol_swap(&tline1,&tline2);
      //printLine(__LINE__);
      //printf("%s",tline1);
    }
  printLine(__LINE__);
  for(i=0;i<nel;i++)
    {
      //printLine(__LINE__);
      sscanf(tline1,"%f",&psi[i]);
      //printLine(__LINE__);
      fgets(tline2, 100, fp);
      //printLine(__LINE__);
      gol_swap(&tline1,&tline2);
      //printLine(__LINE__);
      //printf("%f\n",psi[i]);

    }
  printLine(__LINE__);
  printf("Last psi value is: %f\n",psi[nel-1]);

  free(tline1);
  free(tline2);

  fclose(fp);
}

void shape(float gp[3], float xe[12])
{
  int i;
  //local coordinates
  float r = gp[0];
  float s = gp[1];
  float t = gp[2];
  //shape functions
  N[0] = r; N[1] = s; N[2] = t; N[3] = 1-r-s-t;
  float N_r[4] = {1, 0, 0, -1};
  float N_s[4] = {0, 1, 0, -1};
  float N_t[4] = {0 ,0 , 1, -1};

  float x_r = N_r[0]*xe[0*3+0] + N_r[1]*xe[1*3+0] + N_r[2]*xe[2*3+0] + N_r[3]*xe[3*3+0];
  float x_s = N_s[0]*xe[0*3+0] + N_s[1]*xe[1*3+0] + N_s[2]*xe[2*3+0] + N_s[3]*xe[3*3+0];
  float x_t = N_t[0]*xe[0*3+0] + N_t[1]*xe[1*3+0] + N_t[2]*xe[2*3+0] + N_t[3]*xe[3*3+0];

  float y_r = N_r[0]*xe[0*3+1] + N_r[1]*xe[1*3+1] + N_r[2]*xe[2*3+1] + N_r[3]*xe[3*3+1];
  float y_s = N_s[0]*xe[0*3+1] + N_s[1]*xe[1*3+1] + N_s[2]*xe[2*3+1] + N_s[3]*xe[3*3+1];
  float y_t = N_t[0]*xe[0*3+1] + N_t[1]*xe[1*3+1] + N_t[2]*xe[2*3+1] + N_t[3]*xe[3*3+1];

  float z_r = N_r[0]*xe[0*3+2] + N_r[1]*xe[1*3+2] + N_r[2]*xe[2*3+2] + N_r[3]*xe[3*3+2];
  float z_s = N_s[0]*xe[0*3+2] + N_s[1]*xe[1*3+2] + N_s[2]*xe[2*3+2] + N_s[3]*xe[3*3+2];
  float z_t = N_t[0]*xe[0*3+2] + N_t[1]*xe[1*3+2] + N_t[2]*xe[2*3+2] + N_t[3]*xe[3*3+2];

  //printf("xr %f xs %f xt %f,yr %f ys %f yt %f,zr %f zs %f zt %f\n",x_r,x_s,x_t,y_r,y_s,y_t,z_r,z_s,z_t);
  
  float jacDet = x_r*(y_s*z_t - y_t*z_s) - x_s*(y_r*z_t - y_t*z_r) + x_t*(y_r*z_s - y_s*z_r);
  jac = fabsf(jacDet);

  //Check Jacobian
  if(jac <= 0.0)
    {
      fprintf(stderr, "Negative jacobian, element too distorted!\n");
    }

  float inv_jac[9] = {(y_s*z_t - y_t*z_s)/jacDet, (x_t*z_s - x_s*z_t)/jacDet, (x_s*y_t - x_t*y_s)/jacDet,
		      (y_t*z_r - y_r*z_t)/jacDet, (x_r*z_t - x_t*z_r)/jacDet, (x_t*y_r - x_r*y_t)/jacDet,
		      (y_r*z_s - y_s*z_r)/jacDet, (x_s*z_r - x_r*z_s)/jacDet, (x_r*y_s - x_s*y_r)/jacDet};
  for(i=0;i<4;i++)
    {
      dN[i*3+0] = N_r[i]*inv_jac[0*3+0] + N_s[i]*inv_jac[1*3+0] + N_t[i]*inv_jac[2*3+0];
      dN[i*3+1] = N_r[i]*inv_jac[0*3+1] + N_s[i]*inv_jac[1*3+1] + N_t[i]*inv_jac[2*3+1];
      dN[i*3+2] = N_r[i]*inv_jac[0*3+2] + N_s[i]*inv_jac[1*3+2] + N_t[i]*inv_jac[2*3+2];
    }
  //printf("N is %f %f %f %f\n",N[0],N[1],N[2],N[3]);
  //printf("dN is %f %f %f %f %f %f %f %f %f %f %f %f\n",dN[0],dN[1],dN[2],dN[3],
  //dN[4],dN[5],dN[6],dN[7],dN[8],dN[9],dN[10],dN[11]);
  //printf("jacobian is %f\n",jac);
  
}

void weakform(float xe[12], float Psie, float porosity)
{
  int i,j,k,l;
  // 1 point formula - degree of precision 1
  float gp[3] = {0.25, 0.25, 0.25};
  int w = 1;
  //int ngp = sizeof(w)/sizeof(w[0]);
  int ngp = 1;

  //initialize stiffness matrix
  for(i=0;i<16;i++)
    {
      ke[i] = 0;
    }
  
  //right hand size
  for(i=0;i<4;i++)
    {
      fe[i] = 0;
    }

  //stress strain displacement matrix
  float B[4] = {0,0,0,0};
  //loop over gauss points

  for(i=0;i<ngp;i++)
    {
      float por = porosity;
      shape(gp,xe);
      float z = N[0]*xe[0*3+2] + N[1]*xe[1*3+2] + N[2]*xe[2*3+2] + N[3]*xe[3*3+2];
      if( z > powder_thick)
        {
          por = 0.0;
        }
      for(j=0;j<4;j++)
        {
          B[j] = dN[j*3+2];
        }
      //printf("B is %f %f %f %f\n",B[0],B[1],B[2],B[3]);
      //Transpose of N
      float Ntr[4] = {N[0],N[1],N[2],N[3]};
      //fill ke
      for(k=0;k<4;k++)
	{
	  for(l=0;l<4;l++)
	    {
	      ke[k*4+l] = ke[k*4+l] + Ntr[k] * B[l] * w * jac;
	    }
	}
      //fill fe
      for(k=0;k<4;k++)
	{
	  fe[k] = fe[k] - Ntr[k] * ((por * Psie)/(1 - por * (1 - Psie)))*w*jac;
	}
    }
}

void linearSolve(float M[], float F[],int n)
{
  int i,j,k;
  //int n = 3;
  //float M[9] = {1,1,1,0,2,5,2,5,-1};
  //float F[3] = {6,-4,27};
  //float A[n*6][n*6];
  //printLine(__LINE__);
  /*
  float **A = (float **)malloc(n*4 * sizeof(float *)); 
  for (i=0; i<n; i++) 
    A[i] = (float *)malloc(n*4 * sizeof(float));
  */
  float *A;
  A = calloc(n*n+n,sizeof(float));
  printLine(__LINE__);
  float c,sum=0.0;
  //printf("\nEnter the order of matrix: ");
  //scanf("%d",&n);
  //n = 3;
  //printf("\nEnter the elements of augmented matrix row-wise:\n\n");
  //printf("Elements in the M matrix\n");
  for(i=0; i<n; i++)
    {
      for(j=0; j<=n; j++)
        {
	  //printf("A[%d][%d] : ", i,j);
	  //scanf("%f",&A[i][j]);
	  if(j==n)
	    {
	      A[i*(n+1)+j] = F[i];
	    }
	  else
	    {
	      A[i*(n+1)+j] = M[i*n+j];
	    }
	  //printf("%f\n",A[i*(n+1)+j]);
        }
      //printf("Elements in the M matrix %d\n",i);
    }
  //printLine(__LINE__);
  printf("Read in M and F\n");
  for(j=0; j<n; j++) //loop for the generation of upper triangular matrix
    {
      for(i=0; i<n; i++)
        {
	  if(i>j)
            {
	      c=A[i*(n+1)+j]/A[j*(n+1)+j];
	      for(k=0; k<=n; k++)
                {
		  A[i*(n+1)+k]=A[i*(n+1)+k]-c*A[j*(n+1)+k];
                }
            }
        }
      //printf("Generation of upper triangular matrix %d\n",j);
    }
  printf("Generated upper triangular matrix\n");
  //printLine(__LINE__);
  a_bar[n-1]=A[(n-1)*(n+1)+n]/A[(n-1)*(n+1)+(n-1)];
  /* this loop is for backward substitution*/
  for(i=n-2; i>=0; i--)
    {
      sum=0;
      for(j=i+1; j<n; j++)
        {
	  sum=sum+A[i*(n+1)+j]*a_bar[j];
        }
      a_bar[i]=(A[i*(n+1)+n]-sum)/A[i*(n+1)+i];
      //printf("Backward substitution %d\n",i);
    }
  //printLine(__LINE__);
  printf("\nThe solution is: \n");
  for(i=0; i<n; i++)
    {
      printf("\nx%d=%f\t",i,a_bar[i]); /* x1, x2, x3 are the required solutions*/
    }
  printLine(__LINE__);
}

void luDecomp(float mat[], float F[],int n) 
{
  /*
  printLine(__LINE__);
  float lower[n*n], upper[n*n]; 
  memset(lower, 0, sizeof(lower)); 
  memset(upper, 0, sizeof(upper));
  printLine(__LINE__);
  float y[n], x[n];
  memset(x, 0, sizeof(x)); 
  memset(y, 0, sizeof(y));
  printLine(__LINE__);
  */
  //float b[3] = {12,17,5};

  printLine(__LINE__);
  float *lower = NULL;
  float *upper = NULL;
  lower = calloc(n*n,sizeof(float));
  upper = calloc(n*n,sizeof(float));
  float *y = NULL;
  y = calloc(n*n,sizeof(float));
  printLine(__LINE__);

  /*
  float *mat = NULL;
  mat = calloc(MAX,sizeof(float));
  int i;
  for(i=0;i<n*n;i++) {
    mat[i] = matInput[i];
  }
  */
  
  // Decomposing matrix into Upper and Lower 
  // triangular matrix 
  for (int i = 0; i < n; i++) { 
  
    // Upper Triangular 
    for (int k = i; k < n; k++) { 
      
      // Summation of L(i, j) * U(j, k) 
      float sum = 0; 
      for (int j = 0; j < i; j++) 
	sum += (lower[i*n+j] * upper[j*n+k]); 
      
      // Evaluating U(i, k) 
      upper[i*n+k] = mat[i*n+k] - sum; 
        } 
    
    // Lower Triangular 
    for (int k = i; k < n; k++) { 
      if (i == k) 
	lower[i*n+i] = 1; // Diagonal as 1 
      else { 
	
	// Summation of L(k, j) * U(j, i) 
	float sum = 0; 
	for (int j = 0; j < i; j++) 
	  sum += (lower[k*n+j] * upper[j*n+i]); 
	
	// Evaluating L(k, i) 
	lower[k*n+i] = (mat[k*n+i] - sum) / upper[i*n+i]; 
      } 
    }
    printf("i is %d\n",i);
  } 

  printLine(__LINE__);
  //Display the results for c compatibility
  printf("Lower Triangular\n");
  //Displaying the Results
  for (int i = 0; i < n; i++) { 
    // Lower 
    for (int j = 0; j < n; j++) {
      printf("%f ",lower[i*n+j]);
    }
    printf("\n");
  }
  printf("Upper Triangular\n");
  for (int i = 0; i < n; i++) {
    // Upper 
    for (int j = 0; j < n; j++) {
      printf("%f ",upper[i*n+j]);
    }
    printf("\n");
  }


  //Solve for y intermediate values
  int i, j;
  float sum;
  y[0] = F[0]/lower[0*n+0];
  for(i=1;i<n;i++)
    {
      sum = 0;
      for(j=0;j<n;j++)
	{
	  if(lower[i*n+j] != 0)
	    {
	      sum = sum + lower[i*n+j]*y[j];
	    }
	}
      y[i] = (F[i] - sum)/lower[i*n+i];
    }

  //display y values
  printf("y values for the solutions are: \n");
  for(i=0;i<n;i++)
    {
      printf("%f ",y[i]);
    }
  printf("\n");

  //solve for x values
  a_bar[n-1] = y[n-1]/upper[n*n-1];
  for(i=n-2;i>=0;i--)
    {
      sum = 0;
      for(j=0;j<n;j++)
	{
	  if(upper[i*n+j] != 0)
	    {
	      sum = sum + upper[i*n+j]*a_bar[j];
	    }
	}
      a_bar[i] = (y[i] - sum)/upper[i*n+i];
    }

  //display x values
  printf("x values for the solutions are: \n");
  for(i=0;i<n;i++)
    {
      printf("%f ",a_bar[i]);
    }
  printf("\n");
      
} 

void calcDeformation()
{
  int count = 0;
  int i,j,k,m;

  
  elementsWPow = calloc(nel,sizeof(int));
  elementsXPow = calloc(nel,sizeof(int));
  elementsYPow = calloc(nel,sizeof(int));
  elementsZPow = calloc(nel,sizeof(int));
  for(i=0;i<nel;i++)
    {
      if((coordinatesZ[elementsW[i]] < powder_thick) && (coordinatesZ[elementsX[i]] < powder_thick) &&
	 (coordinatesZ[elementsY[i]] < powder_thick) && (coordinatesZ[elementsZ[i]] < powder_thick))
	{
	  elementsWPow[count] = elementsW[i]; elementsXPow[count] = elementsX[i];
	  elementsYPow[count] = elementsY[i]; elementsZPow[count] = elementsZ[i];
	  count++;
	}
    }
  nel = count;

  /*
  coordinatesXPow = calloc(nnodes,sizeof(float));
  coordinatesYPow = calloc(nnodes,sizeof(float));
  coordinatesZPow = calloc(nnodes,sizeof(float));
  count = 0;
  for(i=0;i<nnodes;i++)
    {
      if(coordinatesZ[i] < powder_thick)
	{
	  coordinatesXPow[count] = coordinatesX[i];
	  coordinatesYPow[count] = coordinatesY[i];
	  coordinatesZPow[count] = coordinatesZ[i];
	  count++;

	  for(j=0;j<nel;j++)
	    {
	      if(elementsWPow[j] == i)
		{
		  elementsWPow[j] = count;
		}
	      if(elementsXPow[j] == i)
		{
		  elementsXPow[j] = count;
		}

	      if(elementsYPow[j] == i)
		{
		  elementsYPow[j] = count;
		}
	      if(elementsZPow[j] == i)
		{
		  elementsZPow[j] = count;
		}
	    }
	}
    }
  nnodes = count;
  */
  for(i=0;i<nel;i++)
    {
      elementsW[i] = elementsWPow[i]; elementsX[i] = elementsXPow[i];
      elementsY[i] = elementsYPow[i]; elementsZ[i] = elementsZPow[i];
    }
  /*
  for(i=0;i<nnodes;i++)
    {
      coordinatesX[i] = coordinatesXPow[i];
      coordinatesY[i] = coordinatesYPow[i];
      coordinatesZ[i] = coordinatesZPow[i];
    }
  */
  printf("New number of nodes is %d, New number of elements is %d\n",nnodes,nel);
  
  
  //count = 0;	    
  float z;
  for(i=0;i<nnodes;i++)
    {
      z = coordinatesZ[i];
      if(fabsf(z - powder_thick) < Tol)
	{
	  count++;
	}
    }
  //printf("count is %d\n",count);
  
  //int *fixnodes = NULL;
  //fixnodes = calloc(count,sizeof(int));
  int fixnodes[count];
  for(i=0;i<count;i++)
    {
      fixnodes[i] = 0;
    }
  count = 0;
  //printf("fixnodes is\n");
  for(i=0;i<nnodes;i++)
    {
      
      z = coordinatesZ[i];
      //printf("%f\n",z);
      if(fabsf(z - powder_thick) < Tol)
	{
	  fixnodes[count] = i;
	  count++;
	}
      //printf("%lu\n",sizeof(fixnodes)/sizeof(fixnodes[0]));
    }
  //===================
  //Assembling ID array
  //===================

  int ID[nnodes];
  for(i=0;i<nnodes;i++)
    {
      ID[i] = 1;
    }
  int ndispl = sizeof(fixnodes)/sizeof(fixnodes[0]);
  //prediscribed displacements
  int nd;
  for(i=0;i<ndispl;i++)
    {
      nd = fixnodes[i];
      ID[nd] = 0;
    }

  //Fill ID array
  count = 0;
  for(j=0;j<nnodes;j++)
    {
      if(ID[j] != 0)
	{
	  ID[j] = count;
	  count++;
	}
      //printf("%d\n",ID[j]);
    }
  printf("assembled ID array\n");

  //=================
  //Generate LM array
  //=================
  int LM[4*nel];
  for(i=0;i<4*nel;i++)
    {
      LM[i] = 0;
    }
  //printf("Length of LM is %lu\n",sizeof(LM)/sizeof(LM[0]));
  for(k=0;k<nel;k++)
    {
      LM[0*nel+k] = ID[elementsW[k]];
      LM[1*nel+k] = ID[elementsX[k]];
      LM[2*nel+k] = ID[elementsY[k]];
      LM[3*nel+k] = ID[elementsZ[k]];
      //printf("LM is %d %d %d %d\n",LM[0*nel+k],LM[1*nel+k],LM[2*nel+k],LM[3*nel+k]);
    }
  /*
  int minLM;
  minLM = 0;
  for(i=0;i<4*nel;i++)
    {
      //printf("LM is %d\n",LM[i]);
      if(LM[i] < minLM)
	{
	  minLM = LM[i];
	}
    }
  printf("Minimum of LM is %d\n",minLM);
  */
  //find max of ID array
  int ndof = ID[0];
  for(i=0;i<nnodes;i++)
    {
      if(ID[i] > ndof)
	{
	  ndof = ID[i];
	}
    }
  printf("max of ID array is %d\n",ndof);

  //displacement vector
  float d[nnodes];
  for(i=0;i<nnodes;i++)
    {
      d[i] = 0;
    }
  printf("generated LM array\n");

  //============================
  //Compute Sparcity
  int nzmax = 0;
  int elem;
  int i_index;
  int j_index;
  for(elem=0;elem<nel;elem++)
    {
      for(k=0;k<4;k++)
	{
	  i_index = LM[k*nel+elem];
	  if(i_index >= 0)
	    {
	      for(m=0;m<4;m++)
		{
		  j_index = LM[m*nel+elem];
		  if(j_index >= 0)
		    {
		      nzmax++;
		    }
		}
	    }
	}
    }
  printf("nzmax is %d\n",nzmax);
  //int i_row[nzmax];
  int *i_row;
  i_row = calloc(nzmax,sizeof(int));
  //int i_col[nzmax];
  int *i_col;
  i_col = calloc(nzmax,sizeof(int));
  printLine(__LINE__);
  for(i=0;i<nzmax;i++)
    {
      i_row[i] = 0;
      i_col[i] = 0;
    }
  //printLine(__LINE__);
  count = 0;
  for(elem=0;elem<nel;elem++)
    {
      for(k=0;k<4;k++)
	{
	  i_index = LM[k*nel+elem];
	  if(i_index >= 0)
	    {
	      for(m=0;m<4;m++)
		{
		  j_index = LM[m*nel+elem];
		  if(j_index >= 0)
		    {
		      i_row[count] = i_index;
		      i_col[count] = j_index;
		      count++;
		    }
		}
	    }
	}
    }
  int irowMin = 0;
  int icolMin = 0;
  int irowMax = 0;
  int icolMax = 0;
  for(i=0;i<nzmax;i++)
    {
      if(i_row[i] < irowMin)
	{
	  irowMin = i_row[i];
	}
      if(i_col[i] < icolMin)
	{
	  icolMin = i_col[i];
	}
      if(i_row[i] > irowMax)
	{
	  irowMax = i_row[i];
	}
      if(i_col[i] > icolMax)
	{
	  icolMax = i_col[i];
	}
      if((i_row[i] < 0) || (i_col[i] < 0))
	{
	  printf("irow is %d and icol is %d\n",i_row[i],i_col[i]);
	}
    }
  printf("irowMin is %d and icolMin is %d and irowMax is %d and icolMax is %d\n",irowMin,icolMin,irowMax,icolMax);
  
  printLine(__LINE__);
  printf("Computed Sparcity\n");

  //============================
  //assembling stiffness matrix
  //============================

  //float K[nzmax];
  float *K;
  K = calloc(nzmax,sizeof(float));
  for(i=0;i<nzmax;i++)
    {
      K[i] = 0;
    }
  //float F[ndof];
  float *F;
  F = calloc((ndof+1),sizeof(float));
  for(i=0;i<(ndof+1);i++)
    {
      F[i] = 0;
    }
  printLine(__LINE__);
  count = 0;
  //loop over elements
  float Psie;
  float xe[12];
  //printf("Last element of LM is %d\n",LM[3*nel+nel+1]);
  for(i=0;i<nel;i++)
    {
      //printLine(__LINE__);
      xe[0] = coordinatesX[elementsW[i]]; xe[1] = coordinatesY[elementsW[i]]; xe[2] = coordinatesZ[elementsW[i]];
      xe[3] = coordinatesX[elementsX[i]]; xe[4] = coordinatesY[elementsX[i]]; xe[5] = coordinatesZ[elementsX[i]];
      xe[6] = coordinatesX[elementsY[i]]; xe[7] = coordinatesY[elementsY[i]]; xe[8] = coordinatesZ[elementsY[i]];
      xe[9] = coordinatesX[elementsZ[i]]; xe[10] = coordinatesY[elementsZ[i]]; xe[11] = coordinatesZ[elementsZ[i]];
      Psie = psi[i];
      //printf("Psie is %f\n",Psie);
      /*
      printf("xe is ");
      for(m=0;m<12;m++)
	{
	  printf("%f ",xe[m]);
	}
      printf("\n");
      */
      weakform(xe,Psie,porosity);
      //printf("fe is %f %f %f %f\n",fe[0],fe[1],fe[2],fe[3]);
      /*
      printf("ke is ");
      for(m=0;m<16;m++)
	{
	  printf("%f ",ke[m]);
	}
      printf("\n");
      */
      //printLine(__LINE__);
      for(j=0;j<4;j++)
	{
	  i_index = LM[j*nel+i];
	  //printLine(__LINE__);
	  //printf("i_index is %d\n",i_index);
	  if(i_index >= 0)
	    {
	      F[i_index] = F[i_index] + fe[j];
	      //printLine(__LINE__);
	      for(k=0;k<4;k++)
		{
		  j_index = LM[k*nel+i];
		  //printLine(__LINE__);
		  if(j_index >= 0)
		    {
		      K[count] = K[count] + ke[j*4+k];
		      //printLine(__LINE__);
		      count++;
		    }
		}
	    }
	}
    }
  //printf("K is %f\n",K[nzmax]);
  /*
  for(i=0;i<nzmax;i++)
    {
      printf("K is %f\n",K[i]);
    }
  */
  /*
  for(i=0;i<ndof;i++)
    {
      printf("F is %f\n",F[i]);
    }
  */
  printLine(__LINE__);
  printf("assembled stiffness matrix\n");
  //ensamble sparse matrix
  //float M[ndof*ndof];
  float *M;
  M = calloc((ndof+1)*(ndof+1),sizeof(float));
  for(i=0;i<=ndof;i++)
    {
      for(j=0;j<=ndof;j++)
	{
	  //printLine(__LINE__);
	  M[i*ndof+j] = 0;
	}
    }
  printf("i_col and i_row value of M is %f\n",M[i_row[nzmax]*ndof+i_col[nzmax]]);
  printLine(__LINE__);
  count = 0;
  //int index;
  //printf("i_row and i_col issues are %d and %d\n",i_row[3397583],i_col[3397583]);
  for(i=0;i<nzmax;i++)
    {
      //count++;
      //index = i_row[i]*(ndof-1)+i_col[i];
      //printf("count is %d index is %d\n",i,index);
      //printLine(__LINE__);
      M[i_row[i]*(ndof+1)+i_col[i]] = M[i_row[i]*(ndof)+i_col[i]] +  K[i];
      //printf("M is %f",M[i]);
    }
  printLine(__LINE__);
  /*
  for(i=0;i<ndof*ndof;i++)
    {
      printf("M is %f\n",M[i]);
    }
  */
  printf("Ensambled sparse matrix\n");
  //a_bar = calloc(ndof, sizeof(float));
  a_bar = calloc(ndof+1, sizeof(float));
  //linearSolve(M,F,ndof);
  luDecomp(M,F,ndof+1);

  free(i_row);
  free(i_col);
  free(K);
  free(F);
  free(M);
}
      

int main()
{
  //char baseName[80] = "/Users/Jacob/Desktop/First_Year_PhD/Research/Powder-Layer-Consolidation/";
  char baseName[80] = "/Users/Jacob/Desktop/First_Year_PhD/Research/Layer_270_03_27_30/0/";

  int myrank = 0;
  int numranks = 4;

  num_ElementsNodes(baseName,myrank);
  offsetCalc(baseName,numranks,myrank);
  readCoordinates(baseName,myrank,nnodes);
  readElements(baseName,myrank,nel);
  readPsi(baseName,myrank,nel);
  calcDeformation();

  free(coordinatesX);
  free(coordinatesY);
  free(coordinatesZ);
  
  free(coordinatesXPow);
  free(coordinatesYPow);
  free(coordinatesZPow);
  
  free(elementsW);
  free(elementsX);
  free(elementsY);
  free(elementsZ);

  free(elementsWPow);
  free(elementsXPow);
  free(elementsYPow);
  free(elementsZPow);
  
  free(a_bar);
  
  return true;
}
