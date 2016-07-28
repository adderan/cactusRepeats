include ${PWD}/sonLib/include.mk


all: ${PWD}/randomGraphCollapse


${PWD}/randomGraphCollapse: ${PWD}/randomGraphCollapse.c ${PWD}/sonLib/lib/sonLib.a
	${cxx} ${cflags} -Werror -I ${PWD}/sonLib/lib -o ${PWD}/randomGraphCollapse randomGraphCollapse.c ${PWD}/sonLib/lib/sonLib.a -lm


${PWD}/sonLib/lib/sonLib.a:
	cd ${PWD}/sonLib && make

clean:
	rm ${PWD}/randomGraphCollapse

