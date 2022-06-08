# Unique RNA sequence generator
This app/library contains a blazingly-fast DNA search algorithm and a console app to generate unique RNA sequences not found in a given transcriptome. Transcriptomes can be given in text form or in binary form. To convert a transcriptome to binary (highly reccomended) run the command `Alloc_RNA MakeTranscriptFile Path/To/File.txt`. The app will load the transcriptome and output a binary representation at `Path/To/File.txt.dat` (each transcript must be on a different line, you can convert FASTA files to this format using `import_fasta.py`).

When you run the app, you will be prompted to enter the length of sequence you want to generate (must be greater than 10) and the path to the transcriptome to load. Enter the path to the text version of the transcriptome, if a binary version exists the app will find it. This can all be done in a single command using the command `Alloc_RNA SeqLength Path/To/Transcriptome [NumExtraSequences ...]` with each extra sequence to check against given as its own command line argument. 

## Including in your own projects
To include the RNA search algorithm, just add `FindDNA.cpp` to your project and include `FindDNA.h`. To include the generator and search algorithm, also compile `Alloc_RNA.cpp`. 

In order to convert DNA to the form used by the search algorithm, call `StringToDNAVector` from `FindDNA.h` on the string you're searching through. If the string you pass to this function could contain anything other than A, C, G, or T (capital or lower-case), you **must** define `NON_DNA_INPUT` project-wide. This is slightly slower, though, so only define it if needed. 
