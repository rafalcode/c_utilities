CC=gcc
CFLAGS=-g -Wall
EXES=cmd2ext extcou

# Command to extension ... allows apllication of a command onto a file with a certain extension. Only useful for rare edge cases.
cmd2ext: cmd2ext.c
	${CC} ${CFLAGS} -o $@ $^

# ext cou ... gives a directory listing but with counts for each file extension
extcou: extcou.c
	${CC} ${CFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
