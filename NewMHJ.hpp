template <typename Function>
float NewMHJ(int N, float *x,Function f)
{
	int n_count=1000; // максимальное число итераций
	int i,count,k;
	double eps_proc=1.e-6,eps_opt=1.e-9; // точность вычислений
	double t1=0.618,t2=1-t1; // константы метода золотого сечения
	double a_m=-100.0,b_m=100.0; // интервал поиска
	double s,ss,lam,d,a,b;
	double m,l,fl,fm;
	float* y1;
	float* y;
	float* x1;
	float* dl;
	y=new float[N];
	x1=new float[N];
	dl=new float[N];
	y1=new float[N];

	d=1.f;  //единичный вектор

//Начальное приближение
  for( i=0; i<N; i++)   x[i]=(float)(rand()-RAND_MAX/2)/(float)RAND_MAX;
  for( i=0; i<N; i++)   y[i]=x[i];

	count=0;
 while(1)
 {  
	count++;
	for( i=0; i<N; i++)  x1[i]=x[i];
//=======================================================
//Одномерная оптимизация методом золотого сечения по координатным направлениям

	for(i=0;i<N;i++) y1[i]=y[i];
	for(k=0;k<N;k++)
	{
	 a=a_m; b=b_m;
	 l=a+t2*(b-a);
	 m=a+t1*(b-a);

	 y1[k]=(float)(y[k]+l);//*d;
	 fl=f(y1);
	 y1[k]=(float)(y[k]+m);//*d;	
	 fm=f(y1);
	
	   do
	   {
		if(fm<fl)
		{
			a=l;  l=m;  m=a+t1*(b-a); fl=fm;
			y1[k]=(float)(y[k]+m);//*d;
			fm=f(y1);
		}
		else
		{
			b=m;  m=l;  l=a+t2*(b-a); fm=fl;
			y1[k]=(float)(y[k]+l);//*d;
			fl=f(y1);

		}

	   } while ((b-a)>eps_opt);
//=========================================================
	 lam=(a+b)/2;
     y[k]=(float)(y[k]+lam*d);
	}

	for( i=0; i<N; i++)   x[i]=y[i];

	s=0.; ss=0;
	for(i=0; i<N; i++)
	{
		s+=(x1[i]-x[i])*(x1[i]-x[i]);
		ss+=x[i]*x[i];
	}
	if( (s < ss*eps_proc) || (count > n_count)) break;

	for( i=0; i<N; i++)   dl[i]=x[i]-x1[i];

	a=0; b=1;
	l=a+t2*(b-a);
	m=a+t1*(b-a);
//==========================================================
//Одномерная оптимизация методом золотого сечения по направлению
	for(i=0;i<N;i++)	y1[i]=(float)(x[i]+l*dl[i]);
	fl=f(y1);
	for(i=0;i<N;i++)	y1[i]=(float)(x[i]+m*dl[i]);
	fm=f(y1);

		do
		{
			if(fm<fl)
			{
				a=l;  l=m;  m=a+t1*(b-a); fl=fm;
				for(i=0;i<N;i++) y1[i]=(float)(x[i]+m*dl[i]);
				fm=f(y1);
			}
			else
			{
				b=m;  m=l;  l=a+t2*(b-a); fm=fl;
				for(i=0;i<N;i++) y1[i]=(float)(x[i]+l*dl[i]);
				fl=f(y1);
			}
		} while ((b-a)>eps_opt);

	lam=(a+b)/2;
	for( i=0; i<N; i++)   y[i]=(float)(x[i]+lam*dl[i]);
	
}//end while(1)

	delete[] y;
	delete[] x1;
	delete[] dl;
	delete[] y1;

return f(x);
}
