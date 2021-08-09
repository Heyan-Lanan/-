//�Ŵ��㷨��������Ż������ͨ���㷨
//������Ľ�x��ʾΪ0-1����Ȼ���������ŵĶ����ƴ���ʹ��Ŀ�꺯��ֵf(x)�ﵽ��С
#include <bits/stdc++.h>
using namespace std;
/** �㷨�������� **/
//����0-1���ĳ��ȣ������ڱ�ʾ��ľ���
#define SLEN 200
//�������������
#define MAXGEN 100
//����С����λ��
#define POINT 10
//�������
#define mProb 1.0 / SLEN
//Ⱥ���С����������Ⱥ���С��һ�������Բ�һ����
#define PSIZE 10
//һ����Ķ���
typedef struct
{
	int x[SLEN]; //x:����Ա�����0-1��
	double y;	 //y=f(x),Ҫ�Ż������Ŀ�꺯��ֵ
} Solution;
//����һ���⼯/�����飺��֮ΪȺ��population
Solution pop[PSIZE * 2];		   //�⼯���������Ӵ����洢�����ǰPSIZE���Ǹ���������PSIZE�����Ӵ�
Solution *current = pop;		   //��ǰ����Ҳ���Ǹ���
Solution *offspring = pop + PSIZE; //�Ӵ��⼯
//��������
void printPop(Solution *p, const char *str);
//��0-1��x����Ϊʵ��*xo ,�ٶ�����4bits��SLEN-4bitsΪС�����ֳ���
void decode(int *x, double *xo)
{

	//�������ֻ���10����
	//С�����ֻ���10����
	double unit = 1;
	*xo = 0;
	for (int i = 0; i < POINT; i++)
	{
		double num = unit * x[POINT - i - 1];
		(*xo) += num;
		unit *= 2;
	}
unit = 0.5;
	for (int i = 0; i < SLEN - POINT - 1; i++)
	{
		double num = unit * x[POINT + i];
		(*xo) += num;
		unit /= 2;
	}
if (x[SLEN - 1])
	{
		*xo = *xo * (-1);
	}
}
//����y=f(x) ,  0-1��x�ĳ��� SLEN
//������������4bits��С������SLEN-4bits
double func1(int *x)
{
	double xo;
	decode(x, &xo); //��0-1��x����������Ľ�xo
	//cout << xo << endl;
	return xo * xo * xo * xo - xo * xo * xo + xo * xo - xo; //����Ŀ�꺯��ֵ
}
//����һ��Ⱥ������н��Ŀ�꺯��ֵy ��ʹ�ú���ָ�룬֧�ָ��������Ż�
void evaluate(Solution *P, double (*ptrf)(int *ptr))
{
	//todo
	//y=f(x)
	for (int i = 0; i < PSIZE; i++)
	{
		P[i].y = ptrf(P[i].x);
	}
}
void print(int *p, int n)
{
	while (n)
	{
		printf("%d", *p);
		n--;
		p++;
	}
	cout << endl
		 << endl;
}
//�㷨��ʼ�������������⼯����Ŀռ䣬�������currentPop�еĽ⣬��������yֵ
void initialize()
{
	for (int i = 0; i < PSIZE; ++i)	   //����currentPop.pop��ÿ����
		for (int j = 0; j < SLEN; ++j) //��ÿ�����0-1�����������
			current[i].x[j] = rand() % 2;
	evaluate(current, func1);
}
//�Ӹ�����ѡ�������⣬ͨ���ӽ����������Ӵ�����
//����������ͨ��PKѡ�����������ѡ��
void crossover()
{ //��������
	int k1, k2, l, k = 0;
	while (k < PSIZE)
	{ //�������Ӵ���һ������
		//todo
		//���ѡ��������������
		k1 = rand() % PSIZE;
		k2 = rand() % PSIZE;
		//���ȷ����������Ⱦɫ�彻��λ��
		l = rand() % SLEN;
		//��������λ��ǰ��������Ӵ�
		memcpy(&(offspring[k].x[0]), &(current[k1].x[0]), sizeof(int) * l);
		memcpy(&(offspring[k + 1].x[0]), &(current[k2].x[0]), sizeof(int) * l);
		memcpy(&(offspring[k].x[l]), &(current[k2].x[l]), sizeof(int) * (SLEN - l));
		memcpy(&(offspring[k + 1].x[l]), &(current[k1].x[l]), sizeof(int) * (SLEN - l));
		k = k + 2;
	}
}
//��offspring�еĸ�����б��죺�������ΪmProb
//��ν�������x[j]��ȡֵ 0-1������ 0 <--> 1
void mutate()
{ //��������
	for (int i = 0; i < PSIZE; ++i)
		for (int j = 0; j < SLEN; ++j)
			if ((rand() % 10000) / 10000.0 < mProb)
				(offspring[i].x[j] == 0) ? offspring[i].x[j] = 1 : offspring[i].x[j] = 0;
}
//��currentPop��offspring��ѡ����һ�����壬�ж���ѡ���㷨������ͨ�������Ȱ�����Ⱥ������õı������� ��Ȼ��
//����1��ѡ����õ�PSIZE��Ϊ��һ�����ض�ѡ��
//����2����ÿ������һ��ѡ����ʣ�yֵС���ã��ı�ѡ��ĸ��ʾ͸ߣ�Ȼ�����ݴ˸��ʷֲ��������PSIZE��
//����3������ѡ�����ѡ��k�����໥pk��������õķ�����һ��������ѡ��PSIZE�� ����ɾ����ѡ���˵ģ�
bool cmp(Solution a, Solution b)
{
	return a.y < b.y;
}
void select(int k)
{				  //ѡ������ �����ý���ѡ��
	double besty; //������ѡ�������Ӵ���yֵ
	int best;	  //������ѡ��������Ӵ��±�
	Solution tmp[PSIZE];
	double prob[2 * PSIZE + 1];
	prob[0] = 0;
	//һ��һ���Ӵ�ѡ��
	//��pop[2*psize]���÷���1��2��3ѡ��һ����Ӧ��ֵ��yֵ���õ�
	if (k == 2)
	{
		double sum = 0;
		for (int i = 1; i <= 2 * PSIZE; i++)
		{
			double p = (1 / (current[i].y + 3));
			prob[i] = p;
			sum += p;
		}
		for (int i = 1; i <= 2 * PSIZE; i++)
		{
			prob[i] /= sum;
			prob[i] += prob[i - 1];
		}
		//cout << prob[2 * PSIZE] << endl;
	}
	if (k == 1)
	{
		sort(current, current + 2 * PSIZE, cmp);
	}
	else
	{
		for (int i = 0; i < PSIZE; ++i)
		{
			if (k == 2)
			{
				double p = (rand() % 1000) / 1000.0;
				//cout << p << endl;
				for (int j = 1; j <= 2 * PSIZE; j++)
				{
					if (p >= prob[j - 1] && prob[j] >= p)
					{
						best = j - 1;
						//cout << best << endl;
						break;
					}
				}
			}
			if (k == 3)
			{
				int k1 = rand() % (2 * PSIZE), k2 = rand() % (2 * PSIZE);
				best = ((current[k1].y < current[k2].y) ? k1 : k2);
			}
			memcpy(&(tmp[i]), &(pop[best]), sizeof(Solution)); //ѡ������Ľ⣬���Ƶ���ʱ�⼯��
		}
		memcpy(current, tmp, sizeof(Solution) * PSIZE);
	}
}
//���Ⱥ�����Ϣ
void printPop(Solution *p, const char *str)
{
	printf("%s/Info: \n", str);
	for (int i = 0; i < PSIZE; ++i)
	{
		printf("Individual %3d : y=%10.6lf=f(", i, p[i].y);
		for (int j = 0; j < SLEN; ++j)
			printf("%d", p[i].x[j]);
		printf(")\n");
	}
	cout << endl
		 << endl;
}
int main()
{
	int seed = 991;
	srand(seed); //������������ӣ�ʹ���㷨�����������
	initialize();
	printf("The %dth Population ", 0);
	printPop(current, "The current Population");
for (int gen = 1; gen < MAXGEN; gen++)
	{
		crossover();
		mutate();
		evaluate(offspring, func1);
		select(1);
	}
	printf("The %dth Population ", MAXGEN);
	printPop(current, "The Final Population");
	return 0;
}