CC=gcc
CFLAGS=-O3
DBGCFLAGS=-g -Wall
EXES=cmd2ext extcou cleangrpo matread genread genread_d txtread txtread_d dreadn dreadn_d vcolfrcr volfrcr_d

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


# dreadn="Data READ Numbers" we want both ints and floats to be read in. OK: we can use an array of unions for that.
# However, there is the problem that sometimes floats are printed as ints.
dreadn: dreadn.c
	${CC} ${DBGCFLAGS} -o $@ $^
dreadn_d: dreadn.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# How to deal with table type text files, which are not easy on the eye? Visual COLumn FoRCeR
# TODO on this:
# Ensure same number of words for each line: so we are dealing with a table.
# First row and first column ... special status?
# For each column find min and max lengths, calculate a STRICSZ for each column.
# Give a max horizontal extent size and accommodate.
# replace common starting letters with z to allow differences to be seen.
vcolfrcr: vcolfrcr.c
	${CC} ${DBGCFLAGS} -o $@ $^
vcolfrcr_d: vcolfrcr.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
