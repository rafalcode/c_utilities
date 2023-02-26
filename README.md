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

