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
	int numprocs;
	//Время запуска программы
	double startwtime = 0.0;
	//Время остановки программы
	double endwtime;
	//Длина имени вычислительного узла
	int namelen;
	//Имя вычислительного узла
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	//===========================================
	// Служебные переменные приложения
	//===========================================
	// Флаг выхода из приложения (0 - работает, 1 - выход)
    int done = 0;
    //Буфер чтения с клавиатуры
    int n;
    //Итератор
    int i;

    //===========================================
    // Служебные переменные вычисления
    //===========================================


    //===========================================
    /**
     * Блок параллельного выполнения программы
     */
    //===========================================

    // Инициализация подсистемы MPI
    MPI_Init(&argc, &argv);
    // Получить размер коммуникатора MPI_COMM_WORLD
    // (общее число процессов в рамках задачи)
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    // Получить номер текущего процесса в рамках
    // коммуникатора MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);

    // Вывод номера потока в общем пуле
    fprintf(stdout, "Process %d of %d is on %s\n", myid,numprocs,processor_name);
    fflush(stdout);


    // Освобождение подсистемы MPI
    MPI_Finalize();
    return 0;
}
