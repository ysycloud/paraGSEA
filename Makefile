#programs,flags,etc.
INCLUDE	=	-I include
TARGET	=	bin/quick_search_serial bin/quick_search_omp bin/quick_search_mpi bin/ES_Matrix_ompi_nocom bin/ES_Matrix_ompi_p2p bin/ES_Matrix_ompi_cocom bin/Cluster_KMeans_ompi

#ALL Phony Targets
.PHONY:	everything	clean	all

#Default starting position
everything:	$(TARGET)
clean:	
	rm -f $(TARGET)
all:	clean everything

bin/quick_search_serial:	src/quick_search_serial.c	\
			src/GSEA.c include/GSEA.h	\
			src/RandomChange.c include/RandomChange.h	\
			src/IO.c include/IO.h
	gcc $(INCLUDE) -g -o bin/quick_search_serial src/quick_search_serial.c src/GSEA.c src/RandomChange.c src/IO.c
	
bin/quick_search_omp:	src/quick_search_omp.c	\
			src/GSEA.c include/GSEA.h	\
			src/RandomChange.c include/RandomChange.h	\
			src/IO.c include/IO.h
	gcc $(INCLUDE) -g -fopenmp -o bin/quick_search_omp src/quick_search_omp.c src/GSEA.c src/RandomChange.c src/IO.c
	
bin/quick_search_mpi:	src/quick_search_mpi.c	\
			src/GSEA.c include/GSEA.h	\
			src/RandomChange.c include/RandomChange.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -o bin/quick_search_mpi src/quick_search_mpi.c src/GSEA.c src/RandomChange.c src/IO.c
	
bin/ES_Matrix_ompi_nocom:	src/ES_Matrix_ompi_nocom.c	\
			src/GSEA.c include/GSEA.h	\
			src/RandomChange.c include/RandomChange.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/ES_Matrix_ompi_nocom src/ES_Matrix_ompi_nocom.c src/GSEA.c src/RandomChange.c src/IO.c

bin/ES_Matrix_ompi_p2p:	src/ES_Matrix_ompi_p2p.c	\
			src/GSEA.c include/GSEA.h	\
			src/RandomChange.c include/RandomChange.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/ES_Matrix_ompi_p2p src/ES_Matrix_ompi_p2p.c src/GSEA.c src/RandomChange.c src/IO.c

bin/ES_Matrix_ompi_cocom:	src/ES_Matrix_ompi_cocom.c	\
			src/GSEA.c include/GSEA.h	\
			src/RandomChange.c include/RandomChange.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/ES_Matrix_ompi_cocom src/ES_Matrix_ompi_cocom.c src/GSEA.c src/RandomChange.c src/IO.c
	
bin/Cluster_KMeans_ompi:	src/Cluster_KMeans_ompi.c	\
			src/GSEA.c include/GSEA.h	\
			src/RandomChange.c include/RandomChange.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/Cluster_KMeans_ompi src/Cluster_KMeans_ompi.c src/GSEA.c src/RandomChange.c src/IO.c
