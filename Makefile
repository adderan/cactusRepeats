include ../sonLib/include.mk


all: randomGraphCollapse


randomGraphCollapse: randomGraphCollapse.c ../sonLib/lib/sonLib.a
	${cxx} ${cflags} -Werror -I ../sonLib/lib -o randomGraphCollapse randomGraphCollapse.c ../sonLib/lib/sonLib.a -lm

clean:
	rm randomGraphCollapse

