CC=gcc
# CFLAGS=-O3
CFLAGS=-g -Wall
DBGCFLAGS=-g -Wall -DDBG
DBG2CFLAGS=-g -Wall -DDBG2
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!

LIBS=-lgsl -lgslcblas -lm
EXES=cmd2ext extcou cleangrpo matread genread genread_d txtread txtread_d parard ahread ahread_d dreadn dreadn_d vcolfrcr volfrcr_d txtread_t bgread bgread2 bgread0 bgread_ bgread0a bgread3 bgmergmc bgmergmc_d bgmergmcstealth bgreadx bgreadx0 bgfiltf contabrd rg0 macsigf bedsumzr bgmergbl bgmergbl2 vttmrg pwmatr tma2pwma pedread dcou dcou2 dcou3 mapedstats mapedstats_d pedcmp pedcmp_d mprd3 mprd3_d mpdmu mpdmu_d dcou4 pedsta pedsta0 mpdmu2 mpdmu2_d mpdmu3 mpdmu3_d mpdmu4 mpdmu4_d tpedsta_d bglsta bglsta_d bglvset vttrd genrd ssrd ssrd_d blard blard_d blard2 blard2_d rurd fard fard_d fard2 paredown f3 f4 pare0 pare2 pare3 onel morel vrd0 chktasty tabcmatch tpedsta2 namsets repinfam csvrd fintimrd csvrdh srtfix sentin sentinu u0 tmsee pread

# Command to extension ... allows apllication of a command onto a file with a certain extension. Only useful for rare edge cases.
cmd2ext: cmd2ext.c
	${CC} ${DBGCFLAGS} -o $@ $^

# ext cou ... gives a directory listing but with counts for each file extension
extcou: extcou.c
	${CC} ${DBGCFLAGS} -o $@ $^
u0: u0.c
	${CC} -o $@ $^
tmsee: tmsee.c
	${CC} ${CFLAGS} -o $@ $^

# test version of cleangrpo, just looks at files in current directory 
cleangrploop: cleangrploop.c
	${CC} ${DBGCFLAGS} -o $@ $^

# increase or decrease timings in an srt file.
srtfix: srtfix.c
	${CC} ${CFLAGS} -o $@ $^

# treat as sentences: an idea to use full stops to parse .. quickly abandoned
# because it won't work .. too often, full stop appear in names and other things ... hardly robust.
# I went straight into no ascii characters and this was me trying to get the prog to play well with UTF8: unfinished.
sentinu: sentinu.c
	${CC} ${CFLAGS} -o $@ $^
sentin: sentin.c
	${CC} ${CFLAGS} -o $@ $^

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
# I had started on mapedstats.c but it turned out only to do ped
# and also to be entirely unscalable, basically reading everything into a pointer-heavy data structure.
pedsta: pedsta.c
	${CC} ${CFLAGS} -o $@ $^

# slurp whole ped (dangerous) so that per SNP comparisons can be made.
pedsta0: pedsta0.c
	${CC} ${CFLAGS} -o $@ $^

# how easy is it to convert the above to tped
tpedsta: tpedsta.c
	${CC} ${CFLAGS} -o $@ $^
tpedsta_d: tpedsta.c
	${CC} ${DBGCFLAGS} -o $@ $^
# old beagle format
bglsta: bglsta.c
	${CC} ${CFLAGS} -o $@ $^
bglsta_d: bglsta.c
	${CC} ${DBGCFLAGS} -o $@ $^
bglvset: bglvset.c
	${CC} ${CFLAGS} -o $@ $^
bglvset_d: bglvset.c
	${CC} ${DBGCFLAGS} -o $@ $^

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
mpdmu4: mpdmu4.c
	${CC} ${CFLAGS} -o $@ $^
mpdmu4_4: mpdmu4.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

# MAP-PED Matchup: print genotypes for the duplicates.
mpdmu: mpdmu.c
	${CC} ${CFLAGS} -o $@ $^
mpdmu_d: mpdmu.c
	${CC} ${DBGCFLAGS} -o $@ $^
mpdmu2: mpdmu2.c
	${CC} ${CFLAGS} -o $@ $^
mpdmu2_d: mpdmu2.c
	${CC} ${DBGCFLAGS} -o $@ $^
mpdmu2_d2: mpdmu2.c
	${CC} ${DBG2CFLAGS} -o $@ $^

# OK ran into trouble (as sual) with what appeard to be a refinement step
# printing out the snpname of the discarded dup, while prting the genotype of its mater dup
# aha, not so easy.
# best thing it to filter out the snapnames before hand.
mpdmu3: mpdmu3.c
	${CC} ${CFLAGS} -o $@ $^
mpdmu3_d: mpdmu3.c
	${CC} ${DBGCFLAGS} -o $@ $^
mpdmu3_d2: mpdmu3.c
	${CC} ${DBG2CFLAGS} -o $@ $^


# for duplicate coutning unit test
dcou: dcou.c
	${CC} ${CFLAGS} -o $@ $^
dcou2: dcou2.c
	${CC} ${DBGCFLAGS} -o $@ $^
# OK, trying to industrialise the successfully working dcou2
dcou3: dcou3.c
	${CC} ${DBGCFLAGS} -o $@ $^

# OK this is custom made for selecting the best tech rep for a duplicated SNP.
# various duplicates are generated.
# ok various Techrep SNPs are simulated, i.e' each SNP has various TR's which give a number 
# of GTs. Majority wins. Precursor to mpdmu.c
dcou4: dcou4.c
	${CC} ${DBGCFLAGS} -o $@ $^

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
# 12.2018 ... restarting this, trying to improve it
txtrd2: txtrd2.c
	${CC} ${DBGCFLAGS} -o $@ $^
txtrd2_d: txtrd2.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^
# Jan 2017, a way of pulling out mm:ss from text, must start with an integer, include a : and not end with a : i.e. time numbers
# not sophis, but fast.
txtread_t: txtread.c
	${CC} -DTNM -o $@ $^
# How about splitting lines on arrow heads?
# # introducing ahread
ahread: ahread.c
	${CC} ${DBGCFLAGS} -o $@ $^
ahread_d: ahread.c
	${CC} ${DBGCFLAGS} -DDBG -o $@ $^

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


#from txtread
vttrd: vttrd.c
	${CC} ${CFLAGS} -o $@ $^
	#use whcar's to read russi
rurd: rurd.c
	${CC} ${CFLAGS} -o $@ $^

genrd: genrd.c
	${CC} ${CFLAGS} -o $@ $^

# when you use ssconvert on a xls file, convert to fasta
ssrd: ssrd.c
	${CC} ${CFLAGS} -o $@ $^
ssrd_d: ssrd.c
	${CC} ${DBGCFLAGS} -o $@ $^

# super fasta blastoutput checker
# multifasta query, chromosome fasta, bloutput reader (fear memory probs)
blard: blard.c
	${CC} ${CFLAGS} -o $@ $^
blard_d: blard.c
	${CC} ${DBGCFLAGS} -o $@ $^

# read in a fasta based on this method
fard: fard.c
	${CC} ${CFLAGS} -o $@ $^
fard_d: fard.c
	${CC} ${DBGCFLAGS} -o $@ $^
fard2: fard2.c
	${CC} ${CFLAGS} -o $@ $^
f3: f3.c
	${CC} ${CFLAGS} -o $@ $^
#f3 abandoned, still can't unstand the transitions
f4: f4.c
	${CC} ${CFLAGS} -o $@ $^

# the creation of a fasta reader demanded joinging lines, this was tricky 
# so I had to go back to first principles and this is the paredown reader
# # with no complications.
paredown: paredown.c
	${CC} ${CFLAGS} -o $@ $^

# the name of the following is CHecK TAble STYle
# not check tasty or chicken tasty or anythign like that.
chktasty: chktasty.c
	${CC} ${CFLAGS} -o $@ $^
# a general form for comparing table file
tabcmatch: tabcmatch.c
	${CC} ${CFLAGS} -o $@ $^
# read in a dbsnp vcf
vrd0: vrd0.c
	${CC} ${CFLAGS} -o $@ $^

# for experimentation:
pare0: pare0.c
	${CC} ${CFLAGS} -o $@ $^
pare2: pare2.c
	${CC} ${CFLAGS} -o $@ $^
# pare4 came from pare2
pare4: pare4.c
	${CC} ${CFLAGS} -o $@ $^
onel: onel.c
	${CC} ${CFLAGS} -o $@ $^
namsets: namsets.c
	${CC} ${CFLAGS} -o $@ $^
pare3: pare3.c
	${CC} ${CFLAGS} -o $@ $^
# and the winner is...
morel: morel.c
	${CC} ${CFLAGS} -o $@ $^

# checking for repat IIDs in fam files
repinfam: repinfam.c
	${CC} ${CFLAGS} -o $@ $^

blard2: blard2.c
	${CC} ${CFLAGS} -o $@ $^
blard2_d: blard2.c
	${CC} ${DBGCFLAGS} -o $@ $^

# read the final time word ... for discogs timing files.
fintimrd: fintimrd.c
	${CC} ${CFLAGS} -o $@ $^

csvrd: csvrd.c
	${CC} ${CFLAGS} -o $@ $^
# not sure what the h is for
csvrdh: csvrdh.c
	${CC} ${CFLAGS} -o $@ $^
# clarifying csvrd, and taking care of empty cells.
csvrde: csvrde.c
	${CC} ${CFLAGS} -o $@ $^

# analyse paragraphs as well
parard: parard.c
	${CC} ${CFLAGS} -o $@ $^
pread: pread.c
	${CC} ${CFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
