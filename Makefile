ANSWER=1
CFLAGS=-Werror -g -Og -std=gnu99

all: rkgrep rkgrep_test bmgrep_test

rkgrep: rkgrep.o bloom.o rkgrep_main.o
	gcc $^ -o $@ -lrt -lm

rkgrep_test: rkgrep_test.o rkgrep.o bloom.o rkgrep_harness.o
	gcc $^ -o $@ -lrt -lm 

bmgrep_test: bmgrep_test.o bmgrep.o bmgrep.h
	gcc $^ -o $@ -lrt -lm

%.o : %.c
	gcc $(CFLAGS) -DANSWER=$(ANSWER) -c ${<}

clean-all: 
	rm -f *.o rkgrep rkgrep_test bmgrep_test
clean :
	rm -f bloom.o rkgrep_test.o rkgrep_main.o rkgrep.o bmgrep.o bmgrep_test.o bmgrep_test rkgrep rkgrep_test 
