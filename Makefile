CC = gcc
CFLAGS = -fopenmp 
LDFLAGS = -lm
EXECS = omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6 jacobi2D-omp gs2D-omp

all: ${EXECS}

omp_solved2: omp_solved2.c
	${CC} ${CFLAGS} omp_solved2.c -o omp_solved2 ${LDFLAGS} 

omp_solved3: omp_solved3.c
	${CC} ${CFLAGS} omp_solved3.c -o omp_solved3 ${LDFLAGS} 

omp_solved4: omp_solved4.c
	${CC} ${CFLAGS} omp_solved4.c -o omp_solved4 ${LDFLAGS} 

omp_solved5: omp_solved5.c
	${CC} ${CFLAGS} omp_solved5.c -o omp_solved5 ${LDFLAGS} 

omp_solved6: omp_solved6.c
	${CC} ${CFLAGS} omp_solved6.c -o omp_solved6 ${LDFLAGS} 

jacobi2D-omp: jacobi2D-omp.c
	${CC} ${CFLAGS} jacobi2D-omp.c -o jacobi2D-omp ${LDFLAGS}

gs2D-omp: gs2D-omp.c
	${CC} ${CFLAGS} gs2D-omp.c -o gs2D-omp ${LDFLAGS} 

clean:
	rm -f ${EXECS}
