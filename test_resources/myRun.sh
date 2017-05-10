#!/bin/bash


nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_testmemory -p 1 -t 12 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> logt1.txt &


for ((i=0;i<5;i++));do

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 1 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log01.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 2 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log02.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 4 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log03.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 8 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log04.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 16 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log05.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 24 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log06.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 32 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log07.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.01 -t 48 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log08.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.05 -t 2 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log11.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.05 -t 4 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log12.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.05 -t 8 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log13.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.05 -t 16 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log14.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.05 -t 24 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log15.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.05 -t 32 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log16.txt &

	nohup yhrun -N 1 -n 1 -p work bin/ES_Matrix_ompi_nocom -p 0.05 -t 48 -a 1 -l 50 -1 data/data_for_test.txt -2 data/data_for_test.txt -o data/ES_Matrix >> log17.txt &

done