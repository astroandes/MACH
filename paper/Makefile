LATEX = pdflatex
MAIN  = paper

all: $(MAIN).pdf

$(MAIN).pdf: 	$(MAIN).tex  
	$(LATEX) $(MAIN).tex 
	bibtex $(MAIN)
	$(LATEX) $(MAIN).tex 
	$(LATEX) $(MAIN).tex 

diff:	$(MAIN)_diff.pdf

$(MAIN)_diff.pdf: $(MAIN).tex
	latexdiff $(MAIN)_old.tex $(MAIN).tex > $(MAIN)_diff.tex
	$(LATEX) $(MAIN)_diff.tex 
	bibtex $(MAIN)_diff
	$(LATEX) $(MAIN)_diff.tex 
	$(LATEX) $(MAIN)_diff.tex 

clean:
	rm -f $(MAIN).pdf $(MAIN).out $(MAIN).bbl $(MAIN).aux

view: 
	open $(MAIN).pdf &

