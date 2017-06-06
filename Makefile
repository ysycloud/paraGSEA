#programs,flags,etc.
INCLUDE	=	-I include
TARGET	=	bin/getReferences bin/quick_search_serial bin/quick_search_omp bin/quick_search_mpi bin/ES_Matrix_ompi_nocom bin/ES_Matrix_ompi_p2p bin/ES_Matrix_ompi_cocom bin/Cluster_KMediods_ompi bin/Cluster_KMediods++_ompi

#ALL Phony Targets
.PHONY:	everything	clean	all

#Default starting position
everything:	$(TARGET)
clean:	
	rm -f $(TARGET)
all:	clean everything
install: 
	install bin/* /usr/bin/
	
bin/getReferences:	src/getReferences.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	gcc $(INCLUDE) -g -o bin/getReferences src/getReferences.c src/GSEA.c src/Tools.c src/IO.c
	
bin/quick_search_serial:	src/quick_search_serial.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	gcc $(INCLUDE) -g -o bin/quick_search_serial src/quick_search_serial.c src/GSEA.c src/Tools.c src/IO.c
	
bin/quick_search_omp:	src/quick_search_omp.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	gcc $(INCLUDE) -g -fopenmp -o bin/quick_search_omp src/quick_search_omp.c src/GSEA.c src/Tools.c src/IO.c
	
bin/quick_search_mpi:	src/quick_search_mpi.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -o bin/quick_search_mpi src/quick_search_mpi.c src/GSEA.c src/Tools.c src/IO.c
	
bin/ES_Matrix_ompi_nocom:	src/ES_Matrix_ompi_nocom.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/ES_Matrix_ompi_nocom src/ES_Matrix_ompi_nocom.c src/GSEA.c src/Tools.c src/IO.c

bin/ES_Matrix_ompi_p2p:	src/ES_Matrix_ompi_p2p.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/ES_Matrix_ompi_p2p src/ES_Matrix_ompi_p2p.c src/GSEA.c src/Tools.c src/IO.c

bin/ES_Matrix_ompi_cocom:	src/ES_Matrix_ompi_cocom.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/ES_Matrix_ompi_cocom src/ES_Matrix_ompi_cocom.c src/GSEA.c src/Tools.c src/IO.c
	
bin/Cluster_KMediods_ompi:	src/Cluster_KMediods_ompi.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/Cluster_KMediods_ompi src/Cluster_KMediods_ompi.c src/GSEA.c src/Tools.c src/IO.c

bin/Cluster_KMediods++_ompi:	src/Cluster_KMediods++_ompi.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/Cluster_KMediods++_ompi src/Cluster_KMediods++_ompi.c src/GSEA.c src/Tools.c src/IO.c