
all:
	g++ main.cpp colors.h colors.cpp matrix.h matrix.cpp sequence.h sequence.cpp type.h -O3 -I/leofs/dbu/zhanghaicang/contact_prediction/colors/lib/include/eigen3 -o colors
	

clean:
	rm -rf *.o
	rm -f colors
