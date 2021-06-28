project:=bandersnatch
all_projects:=$(project)

all: pdf clean

pdf: bandersnatch.tex
	pdflatex bandersnatch.tex
	bibtex bandersnatch
	pdflatex bandersnatch.tex
	pdflatex bandersnatch.tex

clean: 
	rm -vf bandersnatch.aux
	rm -vf bandersnatch.bbl
	rm -vf bandersnatch.blg
	rm -vf bandersnatch.log
	rm -vf bandersnatch.out

