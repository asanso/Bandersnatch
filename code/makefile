project:=bandersnatch
all_projects:=$(project)

all: test bench getparams

curvesearch: small-disc-curves.py small-disc-curves.txt
	sage small-disc-curves.py

getparams: get_params.py
	sage get_params.py

test:
	sage test.py

bench:
	sage bench.py

clean:
	rm -vf params-W.py
	rm -vf params-M.py
	rm -vf params-TE.py
