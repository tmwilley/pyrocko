
all : libparstack.a

parstack.o : parstack.c parstack.h
	gcc -Wextra -fopenmp -fstrict-aliasing -fPIC -O3 -c parstack.c
	
libparstack.a : parstack.o
	rm -f libparstack.a
	ar -csq libparstack.a parstack.o

test_parstack : test_parstack.c
	gcc -Wextra -o test_parstack test_parstack.c -L. -lparstack -lgomp
			
clean :
	rm -f libparstack.a parstack.o
