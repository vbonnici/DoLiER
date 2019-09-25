INCLUDES= -I./dolierlib -I./dolierlib/io -I./dolierlib/sequence -I./dolierlib/index -I./dolierlib/cstyle  -I./dolierlib/clustering  -I./dolierlib/clustering/SWM11  -I./dolierlib/clustering/seq1
#-L/usr/lib/libgsl.a
CC=g++
#CFLAGS=-w -O3
CFLAGS=-O3 -pedantic -Wall
#CFLAGS=-O3
BINDIR=bin
SRCDIR=src

all : dolier-freqs dolier-kfreqs\
	dolier-kmer-list dolier-kmer-list-split\
	dolier-nof dolier-merge-columns\
	dolier-sort-columns dolier-get-column\
	dolier-select dolier-remove-subincl \
	dolier-put-rc dolier-pwm-scores \
	dolier-split-csv-by-k dolier-fvector \
	dolier-clu-seq1 dolier-clu-seq2 \
	dolier-times


dolier-freqs: $(SRCDIR)/dolier-freqs.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-freqs $(SRCDIR)/dolier-freqs.o	-pthread

dolier-kfreqs: $(SRCDIR)/dolier-kfreqs.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-kfreqs $(SRCDIR)/dolier-kfreqs.o	-pthread	

dolier-kmer-list: $(SRCDIR)/dolier-kmer-list.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-kmer-list $(SRCDIR)/dolier-kmer-list.o
	
dolier-kmer-list-split: $(SRCDIR)/dolier-kmer-list-split.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-kmer-list-split $(SRCDIR)/dolier-kmer-list-split.o
	
dolier-nof: $(SRCDIR)/dolier-nof.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-nof $(SRCDIR)/dolier-nof.o	-pthread
	
dolier-merge-columns: $(SRCDIR)/dolier-merge-columns.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-merge-columns $(SRCDIR)/dolier-merge-columns.o
	
dolier-sort-columns: $(SRCDIR)/dolier-sort-columns.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-sort-columns $(SRCDIR)/dolier-sort-columns.o	
	
dolier-get-column: $(SRCDIR)/dolier-get-column.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-get-column $(SRCDIR)/dolier-get-column.o
	
dolier-select: $(SRCDIR)/dolier-select.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-select $(SRCDIR)/dolier-select.o
	
dolier-remove-subincl: $(SRCDIR)/dolier-remove-subincl.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-remove-subincl $(SRCDIR)/dolier-remove-subincl.o
	
dolier-put-rc: $(SRCDIR)/dolier-put-rc.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-put-rc $(SRCDIR)/dolier-put-rc.o
	
dolier-pwm-scores: $(SRCDIR)/dolier-pwm-scores.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-pwm-scores $(SRCDIR)/dolier-pwm-scores.o
	
dolier-split-csv-by-k: $(SRCDIR)/dolier-split-csv-by-k.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-split-csv-by-k $(SRCDIR)/dolier-split-csv-by-k.o
	
dolier-fvector: $(SRCDIR)/dolier-fvector.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-fvector $(SRCDIR)/dolier-fvector.o -lgsl  -lgslcblas

dolier-SWM11: $(SRCDIR)/dolier-SWM11.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-SWM11 $(SRCDIR)/dolier-SWM11.o -lgsl  -lgslcblas
	
dolier-clu-seq1: $(SRCDIR)/dolier-clu-seq1.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-clu-seq1 $(SRCDIR)/dolier-clu-seq1.o -lgsl  -lgslcblas
	
dolier-clu-seq2: $(SRCDIR)/dolier-clu-seq2.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-clu-seq2 $(SRCDIR)/dolier-clu-seq2.o -lgsl  -lgslcblas
	
dolier-clu-seq1-dists: $(SRCDIR)/dolier-clu-seq1-dists.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-clu-seq1-dists $(SRCDIR)/dolier-clu-seq1-dists.o -lgsl  -lgslcblas
	
dolier-times: $(SRCDIR)/dolier-times.o
	$(CC) $(CFLAGS) -o $(BINDIR)/dolier-times $(SRCDIR)/dolier-times.o	
	
.cpp.o: 
	$(CC) $(CFLAGS) -c $< $(INCLUDES) -o $@ 
	
	
clean:
	rm $(SRCDIR)/*.o
	rm $(BINDIR)/dolier-freqs
	rm $(BINDIR)/dolier-kfreqs
	rm $(BINDIR)/dolier-kmer-list
	rm $(BINDIR)/dolier-kmer-list-split
	rm $(BINDIR)/dolier-nof
	rm $(BINDIR)/dolier-merge-columns
	rm $(BINDIR)/dolier-sort-columns
	rm $(BINDIR)/dolier-column-column
	rm $(BINDIR)/dolier-select
	rm $(BINDIR)/dolier-put-rc
	rm $(BINDIR)/dolier-pwm-scores
	rm $(BINDIR)/dolier-split-csv-by-k
	rm $(BINDIR)/dolier-fvector
	rm $(BINDIR)/dolier-SWM11
	rm $(BINDIR)/dolier-clu-seq1
	rm $(BINDIR)/dolier-clu-seq2
	rm $(BINDIR)/dolier-clu-seq1-dists
	rm $(BINDIR)/dolier-times
