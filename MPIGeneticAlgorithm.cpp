/*
 * MPIGeneticAlgorithm.c
 *
 *  Created on: 05.12.2012
 *      Author: djbelyak
 */

//===========================================
/**
  * Блок объявления заголовков
  */
//===========================================

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include "mpi.h"
#include <stddef.h>

using namespace std;
//===========================================
/**
  * Блок объявления фуникций и структур
  */
//===========================================

/**
 * Вещественное случайное число в диапазоне [0,1) (source)
 */
double frand()
{
	return double(rand())/RAND_MAX;
}

/**
 * Структура для представления двоичного генома (source)
 */
struct genom
{
	// Геном
	double* data;
	// Длина генома
	int len;
	genom() { data = 0; }
	~genom() { if( data ) delete[] data; }
};

/**
 * Метод вывода генома на экран
 */
void printGenom(genom In)
{

	for (int i=0; i<In.len; i++)
	{
		fprintf(stdout, "gen[%d]=%8.5f\n", i, In.data[i]);
		fflush(stdout);
	}
}

/**
 * Целевая функция (source)
 */
double eval(genom& x)
{
	double sum = 0;
	for( int i=0; i<x.len; i++ )
		sum += x.data[i]*x.data[i];
	return sum/x.len;
}

/**
 * Создание новой особи со случайным геномом (source)
 */
void create(genom& A, int l)
{
	A.data = new double[l];
	A.len = l;
	for( int i=0; i<l; i++ )
		// случайный выбор гена
		A.data[i] = 200*frand()-100;
}

/**
 * Копирование генома одного индивида (для отбора) (source)
 */
void copy(genom& A, genom& B) // копирование генома одного индивида (для отбора)
{
	for( int i=0; i<A.len; i++ )
		B.data[i] = A.data[i];
}

/**
 * Отбор (турнирная схема) - поединок между двумя индивидами (source)
 */
void select(genom& A, genom& B)
{
	//Вероятность выживания лучшего в паре
	const double pwin = 0.75;
	double fa = eval(A);
	double fb = eval(B);
	double p = frand();
	if( (fa<fb && p<pwin) || (fa>fb && p>pwin) )
		copy(A,B); // победил A
	else
		copy(B,A); // победил B
}

/**
 * Перемешивание популяции (source)
 */
void shuffle(genom* P, int size)
{
	for( int i=0; i<size; i++ )
		swap(P[i].data,P[rand()%size].data);
}

/**
 * Скрещивание двух индивидов (одноточечное) (source)
 */
void crossover(genom& A, genom& B)
{
	int n = A.len;
	//Точка скрещивания
	int k = rand()%(n-1);
	for( int i=k+1; i<n; i++ )
		swap(A.data[i],B.data[i]);
}

/**
 * Мутация особи (source)
 */
void mutate(genom& A)
{
	//Вероятность мутации одного гена
	const double pmutate = 0.1;
	//Степень мутации
    double dx = 0.1;
	for( int i=0; i<A.len; i++ )
		if( frand()<pmutate )
			//Мутация гена
			A.data[i] += dx*(2*frand()-1);
}


/**
 * Главная функция программы (wiki)
 */
int main(int argc, char **argv)
{
	//===========================================
	/**
     * Блок объявления переменных
     */
	//===========================================
	// Служебные переменные MPI
	//===========================================
	//Номер текущего процесса
	int myid;
	//Количество процессов в задаче
	int p;
	//Время запуска программы
	double startwtime = 0.0;
	//Время остановки программы
	double endwtime;
	//Длина имени вычислительного узла
	int namelen;
	//Имя вычислительного узла
	char processor_name[MPI_MAX_PROCESSOR_NAME];

    //===========================================
    // Служебные переменные и константы вычисления
    //===========================================
	//Размерность генома (len)
	const int n = 1000;
	//Размер популяции
	const int size = 4096;
	//Общее количество итераций
	const int tmax = 10000;
	//Количество итерации изолированного развития субпопуляций
	const int dt = 100;
	//Размер миграционной части субпопуляции
	const int fraction = 4;
	//Вероятность проведения поединка отбора
	const double pselect = 0.5;
	//Вероятность скрещивания
	const double pcross = 0.5;
	//Буффер для  особей между субпопуляциями
	double ExchangeBuffer[n*size/fraction];
	//Популяция для обмена
	genom* Exchange = new genom[size/fraction];
	//Лучшие значения целевой функции на каждом шаге алгоритма
	double* best = new double[tmax/dt];
	//Средние по популяции значения целевой функции
	double* average = new double[tmax/dt];



    //Инициализация подсистемы MPI
    MPI_Init(&argc, &argv);
    //Получить размер коммуникатора MPI_COMM_WORLD
    //(общее число процессов в рамках задачи)
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    //Получить номер текущего процесса в рамках
    //коммуникатора MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);

    //Вывод номера потока в общем пуле
    fprintf(stdout, "Process %d of %d is on %s\n", myid, p, processor_name);
    fflush(stdout);


    //===========================================
    /**
      * Блок параллельного выполнения программы
      */
    //===========================================

    //1. (a) Генерация субпопуляции (размер - size/p)
    genom* SubPopulation = new genom[size/p];
    for (int i=0; i<size/p; i++ )
    	create(SubPopulation[i], n);

    for (int extIter=0; extIter<tmax/dt; extIter++) {
    for (int intIter=0; intIter<dt; intIter++)
    {
    	//2. (a) Отбор
    	for (int i=0; i<size/(2*p); i++)
    		if (frand()<pselect)
    			select(SubPopulation[2*i], SubPopulation[2*i+1]);

    	//3. (a) Перемешивание
    	shuffle(SubPopulation, size/p);

    	//4. (a) Скрещивание
    	for(int i=0; i<size/(2*p); i++)
    		if(frand()<pcross)
    			crossover(SubPopulation[2*i],SubPopulation[2*i+1]);

    	//5. (a) Мутация
    	for(int i=0; i<size/p; i++)
    		mutate(SubPopulation[i]);

    	//6. (a) Перемешивание
    	shuffle(SubPopulation, size/p);
    	//7. (a) Goto 2 dt раз
    }

    //8. (a) Отдать часть субпопуляции (размер - size/(p*fraction))
    char* buffer;
    int position=0;
    int bufSize;
    MPI_Pack_size (n*size/(p*fraction), MPI_DOUBLE,
    		MPI_COMM_WORLD, &bufSize);
    buffer = new char[bufSize];
    printf ("bufSize: %d\n", bufSize);
    for (int i=0; i<size/(p*fraction); i++)
    	MPI_Pack(SubPopulation[i].data,SubPopulation[i].len, MPI_DOUBLE,
    		buffer, bufSize, &position, MPI_COMM_WORLD);
    MPI_Gather (buffer, position, MPI_PACKED, ExchangeBuffer, p*position,
    		MPI_PACKED,	0, MPI_COMM_WORLD);

    //9. (m) Сбор буфера
    if (myid == 0)
  	{
    	position = 0;
    	for (int i=0; i<size/(fraction); i++)
    	{
    		create(Exchange[i], n);
    		MPI_Unpack(ExchangeBuffer, sizeof(double)*n*size/fraction,
    	    	&position,Exchange[i].data, n, MPI_DOUBLE, MPI_COMM_WORLD);
    	}
    	//10. (m) Перемешивание
    	shuffle(Exchange, (size/fraction));
    	//11. (m) Статистика
    	average[extIter] = 0;
    	best[extIter] = eval(Exchange[0]);
    	for( int i=0; i<size/fraction; i++ )
    	{
    		double f = eval(Exchange[i]);
    		average[extIter] += f;
    		if( f<best[extIter] )
    			best[extIter] = f;
    	}
    	average[extIter] /= size/fraction;
    }

    }
    //Освобождение подсистемы MPI
    MPI_Finalize();
    return 0;
}
