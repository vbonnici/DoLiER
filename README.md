# DoLiER TOOLS
### Document Listing by Enhanced suffix aRrays TOOLS
An enhanced suffix array (SA + LCP) for indexing collections of documents in order to search and enumerate k-mers in collections of DNA sequences.

<hr />

This projects contains a set of tools used for the analysis of kmers in a set of sequences such as enumeration, frequencies and number of documents in which each k-mer appears.
The library is based on enhanced suffix arrays.


K-mers are searched in an approximated way. The analysis is executed for a pair (k,m) where k is the kmers' length and m is the maximum number of allowed mismatches (at most).

The term frequencies is used with two distinct accessions:<br/>
(1) kmers' frequencies as number of sequences where k-mers appear (document listing problem).<br/>
(2) kmers' frequencies as total number of occurrences of a kmek within the entire collection of sequences.<br/>

Special behaviours apply when dial with reverse complement. See supplementary materials.

FASTA sequences are over the alphabet {A,C,G,T,N}.

<hr />

Read README.txt for a detialed documentation.

<hr />

## Build
Type make -B all to build all the executables, they will be putted in the folder bin.

<hr />

## Tools
All tools' sources are locate in src.
Below a short description of each tool. See "usage" for more detailed descriptions.

Name | Description
----------------------- | -------------
dolier-kmer-list |	enumerate kmers present in a set of FASTA sequences. Kmers are enumerated in lexicographic order. 
dolier-kmer-list-split	|	split a file containing a list of kmers in n different files. Usually, it's used on a dolier-kmer-list output file.
dolier-freqs	|		output kmers frequencies, namely number of sequences where they appear (document listing problem).
 dolier-kfreqs	|		output kmers frequencies, namely total number of occurrences within the entire collection of sequences.
dolier-nof	|		given an input set of FASTA sequences it can be used to get (i) number of sequences inside the file, (ii) number of distinct present kmers (for a fixed length k), (iii) number of positions (given a sequence of length N, it's usually N-k+1 but sequences can contain Ns). 
dolier-sort-columns	|	sort a input n columns file according to values of a specific column.
dolier-get-column	|	extract a specific file from a multi-columns file.
dolier-select	|		filter kmers using thresholds on p-value and number of sequences where the kmer appears.
dolier-split-csv-by-k	|	split a multi-column file, related to kmers of different lengths, according to kmers length. Each output file is a multi-column file containing rows related to the same kmer length.
dolier-remove-subincl	|	given a 1 column file containing kmers of different lengths, output only those kmers that are not a substring of a longer kmer. Duplicated kmers are considered self sub-included and will not be output.
nodup.sh	|		put lines in IFILE into OFILE checking for duplicates. Used to remove duplicated kmers. Note that input file must have one kmer per line.
dolier-put-rc	|		given a file containing a list of kmers, output those kmers and their reverse complements.
dolier-pwm-score	|	given a set of kmers (they may have different length) and a PWM (position weight matrix), output the best score obtained, the average score and the standard deviation.
dolier-clu-seq1 |			sequential kmers clustering algorithm. It works with kmers spectra. dolier-fvector is not required.
dolier-fvector |			output kmers spectra.
