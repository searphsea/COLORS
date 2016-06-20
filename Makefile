#CFLAGS=-I/home4/zhanghaicang/contact_prediction/gurobi605/linux64/include

all:
	g++ main_local.cpp lrs.h lrs.cpp matrix.h matrix.cpp sequence.h sequence.cpp type.h -O3 -I/leofs/dbu/zhanghaicang/contact_prediction/colors/lib/include/eigen3 -o colors_local
	
	#g++ main_local.cpp lrs.h lrs.cpp matrix.h matrix.cpp sequence.h sequence.cpp type.h -O3 -I/leofs/dbu/zhanghaicang/contact_prediction/colors/lib/include/eigen3 -o colors_matrix

clean:
	rm -rf *.o
	rm -f colors
