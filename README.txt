[DoLiER 0.1 TOOLS] [Document Listing by Enhanced suffix aRrays TOOLS]

This projects contains a set of tools used for the analysis of kmers in a set of sequences such as enumeration, frequencies and number of documents in which each k-mer appears.
The library is based on enhanced suffix arrays.


K-mers are searched in an approximated way. The analysis is executed for a pair (k,m) where k is the kmers' length and m is the maximum number of allowed mismatches (at most).

The term frequencies is used with two distinct accessions:
(1) kmers' frequencies as number of sequences where k-mers appear (document listing problem).
(2) kmers' frequencies as total number of occurrences of a kmek within the entire collection of sequences.

Special behaviours apply when dial with reverse complement. See supplementary materials.

FASTA sequences are over the alphabet {A,C,G,T,N}.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------
[BUILD]
Type make -B all to build all the executables, they will be putted in the folder bin.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------
[TOOLS]
All tools' sources are locate in src.
Below a short description of each tool. See "usage" for more detailed descriptions.

- dolier-kmer-list:		enumerate kmers present in a set of FASTA sequences. Kmers are enumerated in lexicographic order. 
- dolier-kmer-list-split:		split a file containing a list of kmers in n different files. Usually, it's used on a dolier-kmer-list output file.
- dolier-freqs:			output kmers frequencies, namely number of sequences where they appear (document listing problem).
- dolier-kfreqs:			output kmers frequencies, namely total number of occurrences within the entire collection of sequences.
- dolier-nof:			given an input set of FASTA sequences it can be used to get (i) number of sequences inside the file, (ii) number of distinct present kmers (for a fixed length k), (iii) number of positions (given a sequence of length N, it's usually N-k+1 but sequences can contain Ns). 
- dolier-sort-columns:		sort a input n columns file according to values of a specific column.
- dolier-get-column		extract a specific file from a multi-columns file.
- dolier-select:			filter kmers using thresholds on p-value and number of sequences where the kmer appears.
- dolier-split-csv-by-k:		split a multi-column file, related to kmers of different lengths, according to kmers length. Each output file is a multi-column file containing rows related to the same kmer length.
- dolier-remove-subincl:		given a 1 column file containing kmers of different lengths, output only those kmers that are not a substring of a longer kmer. Duplicated kmers are considered self sub-included and will not be output.
- nodup.sh			put lines in IFILE into OFILE checking for duplicates. Used to remove duplicated kmers. Note that input file must have one kmer per line.
- dolier-put-rc:			given a file containing a list of kmers, output those kmers and their reverse complements.
- dolier-pwm-score:		given a set of kmers (they may have different length) and a PWM (position weight matrix), output the best score obtained, the average score and the standard deviation.
- dolier-clu-seq1:			sequential kmers clustering algorithm. It works with kmers spectra. dolier-fvector is not required.
- dolier-fvector:			output kmers spectra.



--------------------------------------------------------------------------------------------------------------------------------------------------------------------
[EXAMPLE]

Please, note that the tools dolier-korder and dolier-binom-q are not more inclued in the library.

- f.seqs.fasta is a set of 100 foreground sequences in FASTA format
- b.seqs.fasta is a set of 1000 background sequences in FASTA format

We want: 
1) extract frequencies both for foreground and background set, for k in [6,20] and m in [0,4]
2) get p-values by the binomial test
3) select only those kmers with a pvalue <=0.05 and that appear in at least 20% of foreground sequences
4) take only top 100 best kmers (according to their p-value) for each k
5) cluster such kmers to obtain a restricted set of kmers (centroids)
6) remove those kmers that are sub-string of longer selected kmers 
7) analyse the score over the golden PWM


1) extract frequencies both for foreground and background set, for  k in [6,20] and m in [0,4]
1.1) get foreground sequences and save then in fseqs.freqs
		dolier-freqs f.seqs.fasta f.freqs.6.0 6 0 -r -p 10     (output only frequencies, not kmers, in f.freqs.6.0 for k=6, m=0 using 10 pthreads. See usage for -r meaning)
		dolier-freqs f.seqs.fasta f.freqs.6.1 6 1 -r -p 10
		...
		dolier-freqs fseqs.fasta f.freqs.20.4 20 4 -r -p 10
1.2) since we want to get the background frequencies of only those kmers present in the foreground sequences, we can extract the kmers list from the foreground and search for them in the background.
		dolier-kmer-list f.kmers.6 6  --no-rc enum f.seqs.fasta  (save in f.kmers.6 the list of present 6-mer in f.seqs.fasta)  
		...
		dolier-kmer-list f.kmers.20 20  --no-rc enum f.seqs.fasta
1.3) now get background frequencies of such kmers
		dolier-freqs b.seqs.fasta --kmer-list f.kmers.6 b.freqs.6.0 6 0 -r -p 10 
		dolier-freqs b.seqs.fasta --kmer-list f.kmers.6 b.freqs.6.1 6 1 -r -p 10
		...
		dolier-freqs b.seqs.fasta --kmer-list f.kmers.20 b.freqs.20.4 20 4 -r -p 10
		

Now, we have f.kmers.K (list of K-mers), f.freqs.K.M (foreground frequencies for k=K and m=M), b.freqs.K.M (background frequencies for k=K and m=M)

2) get p-values by the binomial test
		dolier-binomial-q f.kmers.6 f.freqs.6.0 100 f.kmers.6 b.freqs.6.0 1000 binom.6.0  
		dolier-binomial-q f.kmers.6 f.freqs.6.1 100 f.kmers.6 b.freqs.6.1 1000 binom.6.1  
		...
		dolier-binomial-q f.kmers.20 f.freqs.20.4 100 f.kmers.6 b.freqs.20.4 1000 binom.20.4  
		
Now, binom.K.M contains p-values by frequencies for k=K and m=M

Note that, wen dealing with dolier-kfreqs frequencies, you need to use dolier-nof to get the real number of trials, example
		dolier-kfreqs f.seqs.fasta f.kfreqs.6.0 6 0 -r -p 10
		dolier-kmer-list f.kmers.6 6  --no-rc enum f.seqs.fasta
		dolier-kfreqs b.seqs.fasta --kmer-list f.kmers.6 b.kfreqs.6.0 6 0 -r -p 10
		NOF_FG=`dolier-nof f.seqs.fasta positions 6`
		NOF_BG=`dolier-nof b.seqs.fasta positions 6`
		dolier-binomial-q f.kmers.6 f.freqs.6.0 $NOF_FG f.kmers.6 b.freqs.6.0 $NOF_BG binom.6.0  

3) select only those kmers with a pvalue <=0.05 and that appear in at least 20% of foreground sequences
3.1) before to use dolier-select we need to merge kmers, p-values and frequences (that in this case correspond to the number of sequences where each kmer appears), in a single file
		dolier-merge columns  table.6.0 f.kmers.6 binom.6.0 f.freqs.6.0  (merge and save in table.6.0)
		dolier-merge columns  table.6.1 f.kmers.6 binom.6.1 f.freqs.6.1
		...
		dolier-merge columns  table.20.4 f.kmers.20 binom.20.4 f.freqs.20.4
		
3.2) since we need to get top 100 best kmers at step 4, we can sort the kmers now because dolier-select does not modify the input order of printed kmers.
		dolier-sort-columns table.6.0 3 1 table.6.0.sorted   (table.6.0 is a 3 columns file, where pvalues are at column 1)
		dolier-sort-columns table.6.1 3 1 table.6.1.sorted
		...
		dolier-sort-columns table.20.4 3 1 table.20.4.sorted

3.3) now we can run the selection on table.K.M
		dolier-select table.6.0.sorted tselect.6.0 0.05 20  (select kmers by pvalue and document listing count and save them in tselect.6.0)
		dolier-select table.6.1.sorted tselect.6.1 0.05 20
		...
		dolier-select table.20.4.sorted tselect.20.4 0.05 20
		
Now, select.K.M contains such selected kmers (only kmers sequences but not pvalues or other columns).

4) take only top 100 best kmers (according to their p-value) for each k
Since we sorted kmers by their pvalue at step 3.2, now it's simply
		head -n 100 tselect.6.0 > select.6.0
		head -n 100 tselect.6.1 > select.6.1
		...
		head -n 100 tselect.20.4 > select.20.4

5) cluster such kmers to obtain a restricted set of kmers (centroids)
5.1) At first we can merge all the kmers with the same length but obtained with different values of m.
Since, dolier-remove-subincl consideres duplicates kmers as self sub-included, we need nodup.sh to create a diplicates-free file.
		nodup.sh cat select.6.0 tclu.6
		nodup.sh cat select.6.1 tclu.6
		...
		nodup.sh cat select.20.4 tclu.20
5.2) we want cluster kmers having the same length
		dolier-clu-seq1 tclu.6 clu.6 0.1 --weights --log-weights  (cluster the kmers in tclu.6 and put centroids in clu.6,  0.1 --weights --log-weights seams to be the best configuration)
		...
		dolier-clu-seq1 tclu.20 clu.20 0.1 --weights --log-weights

6) remove those kmers that are sub-string of longer selected kmers 
6.1) At first we need tu put all the kmers (for all values of k) into the same file
		cat clu.6 >> tnoincl
		...
		cat clu.20 >> tnoincl
6.2) now we can run the exclusion
		dolier-remove-subincl tnoincl noincl		
		
7) analyse the score over the golden PWM
We want to put such scores in a single file called scores.
7.1) analyse pwm score independently from K
		dolier-pwm-scores noincl scores pwm.file 0.01 -r
7.2) analyse pwm score by k
At first we need to split kmers in noincl by their length
		dolier-split-csv-by-k  noincl 1 0 split
Now, we have kmers of the same length X in split.X
		if(split.6 exists) dolier-pwm-scores split.6 scores pwm.file 0.01 -r
		...
		if(split.20 exists) dolier-pwm-scores split.20 scores pwm.file 0.01 -r
	
	
--------------------------------------------------------------------------------------------------------------------------------------------------------------------
[PARALLELIZATION]
The most expensive step of the analysis is the frequencies retrieving. 
For this reason, dolier-freqs and dolier-kfreqs can run in parallel using posix threads (on a single multicore machine).


Otherwise, one can take advance of the option --kmer-list to split the computation by splitting the kmers list, it's useful if one want to run such job on a cluster of computers.
For example, we used such option to get the frequencies of all kmers in the human genome inside human promoters:
	dolier-kmer-list hg19.kmers.6 6 --no-rc enum chr1.fasta chr2.fasta ... chr22.fasts
	dolier-kmer-list-split hg19.kmers.6 hg19.kmer.6 30  (to obtain 30 files, from hg19.kmers.0 to hg19.kmers.6.29, having more or less the same amount of  consecutive 6-mers)
then run in parallel
	dolier-freqs hg19.promoters.fasta --kmer-list hg19.kmers.6.0 b-freqs.6.0.0 6 0 -r
	dolier-freqs hg19.promoters.fasta --kmer-list hg19.kmers.6.1 b-freqs.6.0.1 6 0 -r
	...
	dolier-freqs hg19.promoters.fasta --kmer-list hg19.kmers.6.29 b-freqs.6.0.29 6 0 -r
Since hg19.kmers.6 is lexicographic ordered as well as the global order of  hg19.kmers.6.0 ... hg19.kmers.6.29,
concatenating all the files hg19.kmers.6.0.X, following the order on X,  will produce a lexicographic ordered file.


--------------------------------------------------------------------------------------------------------------------------------------------------------------------
[REVERSE COMPLEMENT]
Every tool that deals with reverse complements (option -r or --no-rc), expecially dolier-freqs and dolier-kfreqs, makes use of these theorems:
Let
- S a set of nucleotide sequences
- w a k-mer
- rc(w) the reverse complement of w
- rc(S) the reverse complements of S
- dl(w,S) the document listing of w in S, namely the set of sequences where w appears
- |dl(w,S)| the size of such set
- D_k(S) the set of k-mers,of length k, (without duplicates) that appear in S
We want to list every kmers that appears in S or rc(S).
We want to calculate 
- frequency(w,S) = |dl(w,S)  U  dl(w,rc(S))|
- frequency(rc(w),S) = |dl(rc(w),S)  U  dl(rc(w),rc(S))|
We can proof that
- frequency(w,S) = |dl(w,S)  U  dl(rc(w),S)|		so we don't need to store and index rc(S)
- frequency(w,S) = frequency(rc(w),S)				so we can calculate only one of them

Similarly, let
- occ(w,S) the set of occurrences of w in S
We want to calculate
- frequency(w,S) = |occ(w,S)  U  occ(w,rc(S))|
- frequency(rc(w),S) = |occ(rc(w),S)  U  occ(rc(w),rc(S))|
 We can proof that
- frequency(w,S) = |occ(w,S)  U  occ(rc(w),S)|		so we don't need to store and index rc(S)
- frequency(w,S) = frequency(rc(w),S)				so we can calculate only one of them


This last consideration is the reason because kmer lists do not contain reverse complements.



--------------------------------------------------------------------------------------------------------------------------------------------------------------------
[dolier- LIBRARY]
The project's core library is located in dolier-lib.

Often, input sequences are internal converted from the textual 8bit ASCII format to a 3bit format. 
The alphabet of such representation is {A,C,G,T,N,$} where $ indicates the sequence's end (it's used for sequences concatenation).

The base index structure is a suffix array (SA) that can be used to search for a kmer in k*log(n) time.
K-mers can be searched with an arbitrary number of mismatches. The algorithm is a binary search over the SA, bounded by a reverse BWA heuristic (BWA works with FM-index which is essentially and index for the reverse of the input sequence, so kmers are searched for the tail to the head).
SA is combined with the longest common prefix (LCP) array to allow for k-mers enumeration. For each position of the suffix array, it stores the length of le longest common prefix with the previous position.
and additional array to deal with Ns (we don't want enumerate kmers with N inside). For each position of the SA, it stores the distance with the nearest next N.
Finally, a supplementary array, IDS, and other structures are used to solve the document listing problem (see Cs5DLIndex). 
The global index is obtained concatenating all the input sequences and building a SA on such concatenation.
So, for each position of SA, IDS stores the ID of the original sequence.

Sequences are concatenated putting a $ character beetwen them.
So, searching for AAA with 1 mismatches returns any occurrence of AAA with 1 mismatch over the alphabet {A,C,G,T,N}   ( ANA is a valid occurrence at distance 1)
but, for example, A$A is not considered a valid occurrence since it resides between two sequences.
The same applies when a set of k-mers is indexed.



Below, a summary description of the most important core files:

dolier-lib/data_ts			base data types definition
dolier-lib/FlatColumnFileIterator	a generic class to iterate over a 1 column file
dolier-lib/FlatFreqsIterator		a class to iterate over a 2 columns file in the format  "kmer frequency"
dolier-lib/FlatKmersIterator		a class to iterate over a 1 column file containing a set of lexicographic ordered kmers 
dolier-lib/FlatKmersValueIterator	a generic class to iterate over a 2 columns file in the format  "kmer value"
dolier-lib/main_common			a non-OOP (c style) collection of most common operations used in executables sources
dolier-lib/pars_t			an abstract class for command line arguments parsing
dolier-lib/test_suite			a set of functions for binomial tests, or other distribution
dolier-lib/trim				functions implementing trim operations for strings

dolier-lib/clustering/fvector		a collection of functions to work with kmers spectra
dolier-lib/clustering/seq1/seq1		implementation of the clustering algorithm described above

dolier-lib/cstyle/CstyleIndex		a collection of common procedures and high-level objects to work with nucleotide sequences, 3bit format (also called dna5_t) and not only.
dolier-lib/cstyle/dimers 		a collection of functions to work with mono-mers and di-mers (no more used)
dolier-lib/cstyle/pwm_supp		a collection of procedures to deal with PWM

dolier-lib/index/dolier-Index		it's the first version of the index to deal with the document listing problem. It's quite obsolete but faster.
dolier-lib/index/Cs5DLIndex		the class to solve the document listing problem. It's suffix array based.
dolier-lib/index/DNA4Tree		an implementation of a prefix tree to store kmers over the 2bit alphabet {A,C,G,T}.
dolier-lib/index/NSAIterator		a suffix array based class to enumerate kmers and their occurrences
dolier-lib/index/sais			suffix arrays construction, Copyright (c) 2008-2010 Yuta Mori
dolier-lib/index/SASearcher		a class to perform kmer searching over a suffix array without any additional data structure (LCP or document listing support)

dolier-lib/io/FASTAReader		to read FASTA files

dolier-lib/sequence/DNA4WordGenerator	can be used to generate all the words at distance m (number of mismatches) of an input kmer.
dolier-lib/sequence/DNA5Alphabet	implementation of the 3bit alphabet {A,C,G,T,N,$}. It's exentially a bijective map from {A,C,G,T,N,$} to [0,6]
dolier-lib/sequence/WordGenerator	the same of DNA4WordGenerator but can input/output words in std::string type


