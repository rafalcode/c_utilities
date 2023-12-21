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

I find that subtitles are not delays by a constant factor, that it's more complicated than that.
This program helps witht he beginning. Then they just start delaying wildly.

# csvrdm.c
A program to handle the output of the DNA Methylation pipeline. The CSV is similar
to the EPIC Annotation CSV, though there are important changes including renaming of one
or two columns reordering and other things.
The main idea is to take Illumina's multiply annotated CpG's (tables separated by semicolons) and
make them have a row each for easier manipulation. For some reason only known to Illumina
(and perhaps not even them) the multiple annotation often are identical. dupcrd2 is able to get rid of these.
However sometimes they are not and in fact a CpG is assigned to more than one gene. 
Very unfortunately how to resolve them is hard to find out about, nobody can be drawn on the topic.
So after csvrdm has separated them out into row, it is dupcrd2.R's job to resolve.
In dec 23 I decided to follow my line of reason of "probability of biological impact score" rfbisc.
by which a row get points for being in promoter, and in island etc.After running csvrdm on a file and running through dupcrd2.R, you shoulc actually get the same number of lines.
the difference being that they are "resolved". Wrongly or rightly.
the output of csvrd is often "_gened" postfixed. this is a convention, because output goes to STDOUT.

Note:
- The pipeline outputs CSVs with quotes, this should be changed, they are not useful for csvrdm.
- The column CpgGrp was GpgGrp for a while causing untold frustration.
- need dupcrd.R to de-duplicate Cpg annotation 

NOTE: I incorporated rfBISC (RF Biological Impact score) to help decide best Cpg for a certain gene

Big gotcha:
Illumina's annotation has a row for each cpg (naturally) but som eof them are "gene-less" often because islands
don't always occur beside genes and Illumina wnat to annotate the cpgs .. but they're only island-annotations.

So the sequence of programs is: csvrdm then csvfoc (to get 1 line per gene-annotated cpg) and then it goes back into the pipeline,
NOTE a single CPG per gene has to be seelcted.

# kmlrd.c
following from gpxrd.c these are rough and ready parsers of xml. Note I do not say "xml parsers"
becuase they are most certainly not that. 
gpxrd is particularly unrobust, relying on garmin line numbers which are hardcoded.
kmlrd is an improvement, first appearance of a certain characters is looked for.
kmlrd gets stuck on distance from lat and long from a plane. Penzance was 16.5 km away on the surface
but at 11km up in the air, it's 26.8? No it should be about 20.
so the mydist() function is not great at all.
