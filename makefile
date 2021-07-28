project:=bandersnatch
all_projects:=$(project)

all: test bench getparams pdf

curvesearch: 
	sage code/small-disc-curves.py

getparams:
	sage code/get_params.py

test:
	sage code/test.py

bench:
	sage code/bench.py

clean_params:
	rm -vf code/params-W.py
	rm -vf code/params-M.py
	rm -vf code/params-TE.py

pdf:
	git submodule init
	git submodule update
	cd paper && pdflatex bandersnatch.tex
	cd paper && pbibtex bandersnatch.aux
	cd paper && pdflatex bandersnatch.tex
	cd paper && pdflatex bandersnatch.tex

clean_paper: 
	rm -vf paper/bandersnatch.aux
	rm -vf paper/bandersnatch.bbl
	rm -vf paper/bandersnatch.blg
	rm -vf paper/bandersnatch.log
	rm -vf paper/bandersnatch.out

clean: clean_params clean_paper


