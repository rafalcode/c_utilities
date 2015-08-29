CC=gcc
CFLAGS=-O3
DBGCFLAGS=-g -Wall
EXES=cmd2ext extcou cleangrpo matread genread genread_d txtread txtread_d dreadi dreadi_d

# Command to extension ... allows apllication of a command onto a file with a certain extension. Only useful for rare edge cases.
cmd2ext: cmd2ext.c
	${CC} ${DBGCFLAGS} -o $@ $^

# ext cou ... gives a directory listing but with counts for each file extension
extcou: extcou.c
	${CC} ${DBGCFLAGS} -o $@ $^

# test version of cleangrpo, just looks at files in current directory 
cleangrploop: cleangrploop.c
	${CC} ${DBGCFLAGS} -o $@ $^

cleangrpo: cleangrpo.c
	${CC} ${CFLAGS} -o $@ $^

cleangrpo_d: cleangrpo.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

matread: matread.c
	${CC} ${DBGCFLAGS} -o $@ $^

genread: genread.c
	${CC} ${CFLAGS} -o $@ $^
# debug version ... no optimisation, very small buffers.
genread_d: genread.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# a redoing of genread but for for general text
# instead of ignoring whitespace, it must be recorded.
txtread: txtread.c
	${CC} ${DBGCFLAGS} -o $@ $^
txtread_d: txtread.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# dreadi: "Data READ Integers" for text files whose elements are elements are newline-separated and all integers
# actually hived off from 
dreadi: dreadi.c
	${CC} ${DBGCFLAGS} -o $@ $^
dreadi_d: dreadi.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^
.PHONY: clean

clean:
	rm -f ${EXES}
