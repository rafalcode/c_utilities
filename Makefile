CC=gcc
# CFLAGS=-O3
CFLAGS=-g -Wall
DBGCFLAGS=-g -Wall
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!

LIBS=-lgsl -lgslcblas -lm
EXES=cmd2ext extcou cleangrpo matread genread genread_d txtread txtread_d dreadn dreadn_d vcolfrcr volfrcr_d txtread_t bgread bgread2 bgread0 bgread_ bgread0a bgread3 bgmergmc bgmergmc_d bgmergmcstealth bgreadx bgreadx0 bgfiltf contabrd rg0 macsigf bedsumzr bgmergbl bgmergbl2 vttmrg pwmatr tma2pwma pedread

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

# at reading a nxn pairwise matrix.
pwmatr: pwmatr.c
	${CC} ${DBGCFLAGS} -o $@ $^


# FOllowing facination with contingency tables ... this simple ad transparent prog for them
contabrd: contabrd.c
	${CC} ${DBGCFLAGS} -o $@ $^ -lm

rg0: rg0.c
	${CC} ${DBGCFLAGS} -o $@ $^ $(LIBS)

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

# The mac signal variations ... first called bgfilt, but then forked this to create macsigf more specialised.
bgfiltf: bgfiltf.c
	${CC} ${DBGCFLAGS} -o $@ $^

macsigf: macsigf.c
	${CC} ${CFLAGS} -o $@ $^

# filter and merge a bed file into regions defined by another bed file.
bgmergbl: bgmergbl.c
	${CC} ${DBGCFLAGS} -o $@ $^
# same as one above .. allowing depth file now (from samtools depth)
bgmergbl2: bgmergbl2.c
	${CC} ${DBGCFLAGS} -o $@ $^

# wanted to make somethng robust to just three comlumns bedfiles.
bedsumzr: bedsumzr.c
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

genread0: genread0.c
	${CC} ${CFLAGS} -o $@ $^
# debug version ... no optimisation, very small buffers.
genread0_d: genread0.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# using genread0 to do a ped comparator
pedcmp: pedcmp.c
	${CC} ${CFLAGS} -o $@ $^
pedcmp_d: pedcmp.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

mapedstats: mapedstats.c
	${CC} ${CFLAGS} -o $@ $^
mapedstats_d: mapedstats.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# from genread
mprd: mprd.c
	${CC} ${CFLAGS} -o $@ $^
mprd_d: mprd.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# Nah, fork from matread,: after all map is v close to a Nx4 matrix.
mprd2: mprd2.c
	${CC} ${CFLAGS} -o $@ $^
mprd2_d: mprd2.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# Wow, yet nohter one this time forked from szfind
mprd3: mprd3.c
	${CC} ${CFLAGS} -o $@ $^
mprd3_d: mprd3.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# ped files, these are tricky ... have very long lines, let's modify genread to sanuty check them
pedread: pedread.c
	${CC} ${CFLAGS} -o $@ $^
# debug version ... no optimisation, very small buffers.
pedread_d: pedread.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

vttmrg: vttmrg.c
	${CC} ${CFLAGS} -o $@ $^

# debug version ... no optimisation, very small buffers.
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

# based on txtread.c: rough ..reads in and converts to pairwise
# that doesn't uite mean what it says: actually matches n-1 members against n-1 members
# this is quite rigid: it expects the first row to be n-1 strings
tma2pwma: tma2pwma.c
	${CC} -g -Wall -o $@ $^
tma2pwma_d: tma2pwma.c
	${CC} -g -Wall -DDBG -o $@ $^
tma2pwma_d2: tma2pwma.c
	${CC} -g -Wall -DDBG2 -o $@ $^


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

revperm: revperm.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
