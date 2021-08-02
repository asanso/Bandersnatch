project:=bandersnatch
all_projects:=$(project)

all: test bench getparams pdf

curvesearch: 
	sage python-ref-impl/small-disc-curves.py

getparams:
	sage python-ref-impl/get_params.py

test:
	sage python-ref-impl/test.py

bench:
	sage python-ref-impl/bench.py

clean_params:
	rm -vf python-ref-impl/params-W.py
	rm -vf python-ref-impl/params-M.py
	rm -vf python-ref-impl/params-TE.py

pdf:
	git submodule init
	git submodule update
	cd paper && pdflatex bandersnatch.tex
	cd paper && bibtex bandersnatch.aux
	cd paper && pdflatex bandersnatch.tex
	cd paper && pdflatex bandersnatch.tex

clean_paper: 
	rm -vf paper/bandersnatch.aux
	rm -vf paper/bandersnatch.bbl
	rm -vf paper/bandersnatch.blg
	rm -vf paper/bandersnatch.log
	rm -vf paper/bandersnatch.out

clean: clean_params clean_paper


