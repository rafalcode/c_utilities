# c\_utilities
Various input/output and system utilities
(note: used to be called "sysaduts" system administration
utilites, but realised they were becomeign more general than that)

# Programs
* matread: reads a text which consists of a raw matrix of arbitrary size. floats assumed all round.
* genread: parses a text file into lines and words.
* dreadi: this is a pretty robust records-set-out-as-paragraphs datafile reader. It assumes that each record starts with
a name (string) and then a series of integers, which can be any number, but of course must be the same for each record.

# TODO
I'd like to deal with filelistings files.
There would have only one word per line.


# mprd3.c
* a plink-alike stats reader
* ran into issue of hashing on chr/pos strings, which is fine for
* plink must be 

# srtfix and vttgo3
And similar sound file are for YT downloaded autosubs (not manual subs)
srtfix (a delay program) actually interprets the timings .. not so sophicated mind you 
ebause they are all the same format in YT VTT so not strtok enecessary.
vttgo also very simplistic, again the lines are highly regular so you 
can rey on indices.
Note no special requirements for Cyrillic, this must be UTF-8 at work behind the scenes.

# csvrdm.c
A program to handle the output of the DNA Methylation pipeline. The CSV is similar
to the EPIC Annotation CSV, though there are important changes including renaming of one
or two columns reordering and other things. Note:
- The pipeline outputs CSVs with quotes, this should be changed, they are not useful for csvrdm.
- The column CpgGrp was GpgGrp for a while causing untold frustration.
- need dupcrd.R to de-duplicate Cpg annotation 

# kmlrd.c
following from gpxrd.c these are rough and ready parsers of xml. Note I do not say "xml parsers"
becuase they are most certainly not that. 
gpxrd is particularly unrobust, relying on garmin line numbers which are hardcoded.
kmlrd is an improvement, first appearance of a certain characters is looked for.
kmlrd gets stuck on distance from lat and long from a plane. Penzance was 16.5 km away on the surface
but at 11km up in the air, it's 26.8? No it should be about 20.
so the mydist() function is not great at all.
