comp = gfortran

pattern = *.f95
source = ${wildcard ${pattern}}
obj = ${patsubst %.f95, %.o, ${source}}


main: ${obj}
	${comp} $^ -o $@

%.mod %.o: %.f95
	touch $@
	${comp} -c $<


code_odin.o: calculation.mod 


run: main
	./main
result: main input
	./main<input>result
new: 
	touch code.f95 fun.f95 sub.f95
clear:
	rm -f *.o *.out *.mod main*
