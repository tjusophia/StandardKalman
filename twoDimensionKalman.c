/*kalman  two*/
#include <stdio.h>
#include <math.h>
#define  N 1
#define  M 2


typedef struct {
	float  Q[M][M];         /* 过程噪声方差                   */
	float  R[M][M];         /* 测量噪声方差                   */
	float  A[M][M];				  /* 状态方程系统参数               */
	float  B[M][N];				  /* 状态方程系统参数               */
	float  H[M][M];				  /* 量测方程系统参数               */
	float  x_estimate[M][N];      /*结果输出                        */
	int FirstFlag;                /* 首次进行滤波标志位，1-首次运行 */
} tTwoDimKalmanFcb;
tTwoDimKalmanFcb  gtAccelerationFcbTwo;    /* 二维滤波器控制块 */
float ans[M][M];                 /*求逆方阵输出结果*/
float KalmanFilterInit(tTwoDimKalmanFcb* pfcb)
{
	float t = 0.07;
	pfcb->FirstFlag = 1;
	pfcb -> A[0][0] = 1.0; pfcb->A[0][1] = t;
	pfcb->A[1][0] = 0.0; pfcb->A[1][1] = 1.0;

	pfcb->B[0][0] = 0.5*pow(t,2); pfcb->B[1][0] = t;

	pfcb->H[0][0] = 1.0; pfcb->H[0][1] = 0.0;
	pfcb->H[1][0] = 0.0; pfcb->H[1][1] = 1.0;

	pfcb->Q[0][0] = 0.9; pfcb->Q[0][1] = 0.1;
	pfcb->Q[1][0] = 0.1; pfcb->Q[1][1] = 0.5;

	pfcb->R[0][0] = 0.9; pfcb->R[0][1] = 0.5;
	pfcb->R[1][0] = 0.5; pfcb->R[1][1] = 0.9;
}
/*MXM阶矩阵求逆矩阵函数……begin*/
/*void Matrix_inverse(float arc[M][M], int n, float ans[M][M])//计算矩阵的逆//列*/
void Matrix_inverse(float (*arc)[M])/*arc是一维指针，指向大小为M的数组，因此叫数组指针*/
    {
	int i, j, k;	
    float max, tempA, tempB, P;
    int max_num;	
    float arcs[M][M];
	memset(ans,0,M*M*sizeof(float));
	for(i=0;i<M;i++)
		for (j = 0; j < M; j++)
		{
			arcs[i][j]= arc[i][j];
		}
   /* memcpy(arcs, arc, 288);	*/
      for (i = 0; i < M; i++)	
         {	ans[i][i] = 1;	
	      }	
      for (i = 0; i < M; i++)//第i列	
      {		
		  max = fabs(arcs[i][i]);
          max_num = i;		
          for (j = i + 1; j < M; j++)//选出主元，		
           {			
		     if (fabs(arcs[j][i]) > max)			
             {	
			  max = fabs(arcs[j][i]);		
              max_num = j;			
             }	
            }

  if (max == 0)		
  {			printf("i can't");			
  break;	
  return;
  }
      for (k = 0; k < M; k++)//交换行		
         {	
		  tempA = arcs[i][k];			
          arcs[i][k] = arcs[max_num][k];			
          arcs[max_num][k] = tempA;			
          tempB = ans[i][k];			
          ans[i][k] = ans[max_num][k];			
          ans[max_num][k] = tempB;		
          }		

      for (k = i + 1; k < M; k++)		
         {	
		  P = arcs[k][i] / arcs[i][i];			
          for (j = 0; j < M; j++)			
           {
			  arcs[k][j] = arcs[k][j] - arcs[i][j] * P;				
              ans[k][j] = ans[k][j] - ans[i][j] * P;			
            }		
          }
	
      }	
	  for (i =0; i < M;i++)
	  {
		  for (k = i + 1; k < M; k++)
		  {
			  P = arcs[i][k] / arcs[k][k];
			  for (j = 0; j < M; j++)
			  {
				  arcs[i][j] = arcs[i][j] - arcs[k][j] * P;
				  ans[i][j] = ans[i][j] - ans[k][j] * P;
			  }
		  }
	  }
      for (i = 0; i < M; i++)//行	
      {	
		  P = arcs[i][i];		
		  arcs[i][i] = arcs[i][i] / P;
	      
		  for (j = 0; j < M; j++)
		  {
			  ans[i][j] = ans[i][j] / P;
		  }
         
      }	
	  
	  /*for (i = 0; i < M; i++)
		  for (j = 0; j < M; j++)
		  {
			  printf("%f,",arcs[i][j]);
			 
		  }
	  for (i = 0; i < M; i++)
		  for (j = 0; j < M; j++)
		  {
			 
			  printf("%f,", ans[i][j]);
		  }*/
	  return;
     }
/*MXM阶矩阵求逆矩阵函数……end*/

/*********************************************************************************************************
** Function name:           KalmanFilterTwoDim
** Descriptions:            二维卡尔曼滤波函数-hcs
** input parameters:        pfcb        : 滤波控制块
**                          meas_value  : 输入的量测值
**                          contr_input : 输入的控制量
** output parameters:       none
** Returned value:          滤波值
*********************************************************************************************************/
float KalmanFilterTwoDim(tTwoDimKalmanFcb* pfcb, float (*meas_value)[M], float (*contr_input)[N])
{
	
	float x_predict[M][N],  p_predict[M][M], x_estimate[M][N], p_estimate[M][M], K[M][M],Z[M][N],U[N][N],C[M][M],D[M][M],E[M][N],F[M][N];
	int n, m, l,i,j;
	static float Last_p_estimate[M][M], Last_x_estimate[M][N], Q[M][M], R[M][M], H[M][M], H_rotate[M][M], A_rotate[M][M], A[M][M], B[M][N];
		/*M=状态向量的个数，N=1状态向量的列数一般是1,L=M*/
	memset(x_predict, 0, M * N * sizeof(float));
	memset(p_predict, 0, M * M * sizeof(float));
	memset(x_estimate, 0, M * N * sizeof(float));
	memset(p_estimate, 0, M * M * sizeof(float));
	memset(K, 0, M * M * sizeof(float));
	memset(C, 0, M * M * sizeof(float));
	memset(D, 0, M * M * sizeof(float));
	memset(E, 0, M * N * sizeof(float));
	memset(F, 0, M * N * sizeof(float));

	memcpy(Z, meas_value, M * N * sizeof(float));
	memcpy(U, contr_input, N * N * sizeof(float));
	if (pfcb->FirstFlag == 1)
    {
		memcpy(x_estimate, Z, M * N * sizeof(float));
		memcpy(Last_x_estimate, Z, M * N * sizeof(float));
		memcpy(Q, pfcb->Q, M * M * sizeof(float));
		memcpy(Last_p_estimate, Q, M * M * sizeof(float));
		memcpy(R, pfcb->R, M * M * sizeof(float));
		memcpy(H, pfcb->H, M * M * sizeof(float));
		memcpy(A, pfcb->A, M * M * sizeof(float));
		memcpy(B, pfcb->B, M * N * sizeof(float));
	/*for (i = 0; i < M; i++)
		for(j=0;j<N;j++)
	    {
			x_estimate[i][j]      = Z[i][j];
			Last_x_estimate[i][j] = Z[i][j];
	    }

	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
		{
			Q[i][j]=pfcb->Q[i][j];
			Last_p_estimate[i][j] = Q[i][j];
			p_estimate[i][j] = Q[i][j];
			R[i][j] = pfcb->R[i][j];
			H[i][j] = pfcb->H[i][j];
		}

	/*求H和A转置*/
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
		{
			H_rotate[m][n] = H[n][m];
			A_rotate[m][n] = A[n][m];
		}

		pfcb->FirstFlag = 0;
		memcpy(pfcb->x_estimate, x_estimate, M * N * sizeof(float));
		return;
	}
	
	/*预测模型计算 X(k|k-1)=A*X(k-1|k-1)+B*U(k) */
	for (m = 0; m < M; m++)
		for (n = 0; n < N; n++)
		{
			for (l = 0; l < M; l++)
			{
				x_predict[m][n] += A[m][l] * Last_x_estimate[l][n];
			}
			for(l=0;l<N;l++)
			{
				C[m][n] += B[m][l] * U[l][n];
			}
			x_predict[m][n] = x_predict[m][n] + C[m][n];
		}

	/*预测协方差矩阵计算P(k|k-1)=A*P(k-1|k-1)A'+Q */
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
			for (l = 0; l < M; l++)
			{
				C[m][n]+= A[m][l] * Last_p_estimate[l][n] ;
			}

	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
		{
			for (l = 0; l < M; l++)
			{
				p_predict[m][n] += C[m][l] * A_rotate[l][n];
			}
			p_predict[m][n] = p_predict[m][n] + Q[m][n];
		}

	/*kalman增益矩阵计算K=P(k|k-1)*H'/(H*P(k|k-1)*H'+R) */

	/*求H*P*/
	memset(C, 0, M * M * sizeof(float));
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
		{
			for (l = 0; l < M; l++)
			{
				C[m][n] += H[m][l] * p_predict[l][n];
			}
		}
	/*求H*P*H'+R*/
	memset(D, 0, M* M * sizeof(float));
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
		{
			for (l = 0; l < M; l++)
			{
				D[m][n] += C[m][l] * H_rotate[l][n];
			}
			D[m][n] = D[m][n] + R[m][n];
		}
	/*求H*P*H'+R的逆矩阵:调用求逆函数，逆矩阵保存在ans[][]数组中*/
	Matrix_inverse(D);
	/*求P*H'*/
	memset(D, 0, M* M * sizeof(float));
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
		{
			for (l = 0; l < M; l++)
			{
				D[m][n] += p_predict[m][l]*H_rotate[l][n];
			}
		}
	/*求P*H'*(H*P*H'+R的逆矩阵)/*/
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
		{
			for (l = 0; l < M; l++)
			{
				K[m][n] += D[m][l] * ans[l][n];
			}
		}

	/*更新状态估计值：X(k|k)=X(k|k-1)+K(k)*(Z(k)-H*X(k|k-1)) */
	for (m = 0; m < M; m++)
		for (n = 0; n < N; n++)
		{
			for (l = 0; l < M; l++)
			{
				E[m][n] += H[m][l] * x_predict[l][n];
			}
			E[m][n] = Z[m][n] - E[m][n];
		}
	for (m = 0; m < M; m++)
		for (n = 0; n < N; n++)
		{
			for (l = 0; l < M; l++)
			{
				F[m][n] += K[m][l] * E[l][n];
			}
			x_estimate[m][n] = x_predict[m][n] + F[m][n];
		}
	/*更新协方差矩阵：P(k|k)=(I-K(k)*P(k|k-1))               */
	memset(C, 0, M* M * sizeof(float));
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
		{
			if (m == n)
				C[m][n] = 1 - K[m][n];
			else
				C[m][n] = -K[m][n];
		}
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
			for (l = 0; l < M; l++)
			{
				p_estimate[m][n]+= C[m][l] * p_predict[l][n];
			}
		
	memcpy(pfcb->x_estimate, x_estimate, M* N * sizeof(float));
}



main()
{
	int i, j, m, n, l;

	float arc[M][M] = { 0,-2,1,3 };
	float RelDist[31] = { 24.7000,   24.7000,   23.5000,   23.5000,   22.2000,   22.2000,   21.0000,   21.0000,   21.0000,   21.0000,   19.9000,   19.9000,   18.8000,   17.8000,   17.8000,   16.9000,   16.9000,   16.9000,   16.9000,   16.0000,   16.0000,   15.2000,   15.2000,   14.5000,   13.8000,   13.8000,   13.8000,   13.8000,   13.2000,   13.2000 ,12.5000 };
	float RelVel[31] = {   -9.4000,   - 9.4000,   - 9.2000,   - 9.2000,   - 8.9000,   - 8.9000,   - 8.5000,   - 8.5000,   - 8.5000,   - 8.5000,
                          - 8.2000,   - 8.2000,   - 7.9000,   - 7.5000,   - 7.5000,   - 7.2000,   - 7.2000,   - 7.2000,   - 7.2000,   - 6.9000,   
		                  - 6.9000,   - 6.6000,   - 6.6000,   - 6.3000,   - 6.0000,   - 6.0000,   - 6.0000,   - 6.0000,   - 5.6000,   - 5.6000,   - 5.3000 };
	float meas_value[M][N];
	float contr_input[N][N];
	float t = 0.07;
	float last_Vel = 0;
	float a;
	tTwoDimKalmanFcb  gtAccelerationFcbTwo;    /* 二维滤波器控制块 */
	KalmanFilterInit(&gtAccelerationFcbTwo);
	for (i = 0; i < 31; i++)
	{
		meas_value[0][0] = RelDist[i];
		meas_value[1][0] = RelVel[i];
		a = (RelVel[i] - last_Vel) / t;
		contr_input[0][0] = a;
		last_Vel = RelVel[i];
		KalmanFilterTwoDim(&gtAccelerationFcbTwo, meas_value, contr_input);
		
		for(m=0;m<M;m++)
			for (n = 0; n < N; n++)
			{
				printf("%f ,  ", gtAccelerationFcbTwo.x_estimate[m][n]);
			}
		printf("%f ,  ", RelDist[i]);
		printf("%f ,\n  ", RelVel[i]);
	}

}   
	
