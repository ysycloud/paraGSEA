#programs,flags,etc.
INCLUDE	=	-I include
TARGET	=	bin/ES_Matrix_ompi_testmemory

#ALL Phony Targets
.PHONY:	everything	clean	all

#Default starting position
everything:	$(TARGET)
clean:	
	rm -f $(TARGET)
all:	clean everything
install: 
	install bin/* /usr/bin/


bin/ES_Matrix_ompi_testmemory:	src/ES_Matrix_ompi_testmemory.c	\
			src/GSEA.c include/GSEA.h	\
			src/Tools.c include/Tools.h	\
			src/IO.c include/IO.h
	mpicc $(INCLUDE) -g -fopenmp -o bin/ES_Matrix_ompi_testmemory src/ES_Matrix_ompi_testmemory.c src/GSEA.c src/Tools.c src/IO.c