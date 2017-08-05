CC=gcc
CFLAGS=-O3
DBGCFLAGS=-g -Wall
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!
EXES=cmd2ext extcou cleangrpo matread genread genread_d txtread txtread_d dreadn dreadn_d vcolfrcr volfrcr_d txtread_t bgread bgread2 bgread0 bgread_ bgread0a bgread3 bgmergmc bgmergmc_d bgmergmcstealth bgreadx bgreadx0

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

# if there was ever a program made in hell it was this.
# it's called matread, but in fact it's an array reader. A 1-D matrix..
matread: matread.c
	${CC} ${DBGCFLAGS} -o $@ $^


# for readng bedgraph files
# this one (not the 0,2 qand 3) seems to work
bgread: bgread.c
	${CC} ${DBGCFLAGS} -o $@ $^

# OK try to make it read irrespective of column nums
bgreadx: bgreadx.c
	${CC} ${DBGCFLAGS} -o $@ $^

# OK. ther's a problem ... I was doign to this for intersect Bedgraph .. but it's coutnign repeated entries of course because htere's no control
# so this effort is to control ... only to include entries that are not repeated (judgin by the ranges)
bgreadx0: bgreadx0.c
	${CC} ${DBGCFLAGS} -o $@ $^

# poor start on parsing out different chromosome.
# bgread2 worked inside the process function
bgread2: bgread2.c
	${CC} ${DBGCFLAGS} -o $@ $^

# This is the space inefficient version. BUT it's the one that works! Yes, this is the primary one.
# Note that if given a bed file with column 4 as a string or with more fields, they will be ignored.
bgread0: bgread0.c
	${CC} ${DBGCFLAGS} -o $@ $^
# bgread3 is where I try to merge contiguous sections.
# sorry can't do it, too difficult.
# # needs full immersion. Cna't do it in an hour.
bgread3: bgread3.c
	${CC} ${DBGCFLAGS} -o $@ $^
# Actually bgread0 itself is pretty good. Can you believe 
# I didn't realize I had done soem decent stuff on it.
# Getting much closer to a bgmerge now:
# behold bgmergmc
bgmergmc: bgmergmc.c
	${CC} ${DBGCFLAGS} -o $@ $^
bgmergmc_d: bgmergmc.c
	${CC} ${TDBGCFLAGS} -o $@ $^
bgmergmcstealth: bgmergmcstealth.c
	${CC} ${DBGCFLAGS} -o $@ $^

bgread0a: bgread0a.c
	${CC} ${DBGCFLAGS} -o $@ $^
# Ooops again failed.
# bgread0_.c tries to render a container for our bgr_t's
bgread0_: bgread0_.c
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
# Jan 2017, a way of pulling out mm:ss from text, must start with an integer, include a : and not end with a : i.e. time numbers
# not sophis, but fast.
txtread_t: txtread.c
	${CC} -DTNM -o $@ $^


# dreadn="Data READ Numbers" we want both ints and floats to be read in. OK: we can use an array of unions for that.
# However, there is the problem that sometimes floats are printed as ints.
# by git ver  e37b9163 these two were no longer working.
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

mptl0: mptl0.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^ -lpython2.7

szfread: szfread.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
