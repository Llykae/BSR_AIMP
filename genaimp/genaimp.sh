#MISE EN PLACE LA LIBRAIRIE libint-1.1.5
#lib=libint-1.1.5

#lib=../../libint-1.1.5
#if !(test -d ${lib})
#	then
#	 tar -xzvf	${lib}.tar.gz
#	fi
#if !(test -e ${lib}/Obj/lib/libint.a)
#	then
#	 mkdir		${lib}/Obj
#	 cd		${lib}/Obj
#	 ../configure 													\
#	 --with-libint-max-am=3  --with-libderiv-max-am1=2 --with-libderiv-max-am2=1 --with-max-class-size=11		\
#	 --enable-deriv=yes --enable-long-double=no --enable-r12=no --with-cc-optflags=-O3 --with-cxx-optflags=-O3 --enable-debug=no --enable-shared=no
#
#	 make
#	 sudo make install
#	 cd ../..
#	fi


#EXECUTION DU PROGRAMME
prog=genaimp
#prog=derivtst
#pg="-pg"		#Option pour regarder le temps passe dans chaque fonction du programme "-pg"
#Opt="-O3"		#Option d'optimisation "-O3"
#Warn="-Wall"		#Option pour message d'alerte a la compilation "-Wall"
Warn="-Wno-unused-result"		#Option pour message d'alerte a la compilation  sans message de non-utilisation


gcc ${Opt} ${Warn} ${pg} ${prog}.c			\
-I/usr/local/libint/1.1.5-stable/include/		\
/usr/local/libint/1.1.5-stable/lib/libderiv.a		\
/usr/local/libint/1.1.5-stable/lib/libint.a		\
-lstdc++						\
-o ${prog}.exe -lm

#./${prog}.exe

#if (test "${pg}" = "-pg")
#	then
#	 gprof ./${prog}.exe
#	 rm gmon.out
#	fi

#Nettoyage
#rm ${prog}.exe
#rm basis.dat 2> /dev/null
