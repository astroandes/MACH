LATEX = pdflatex
MAIN  = article

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
	rm -f $(MAIN).{aux,bbl,ps,pdf,dvi,blg,log,out}
	rm -f $(MAIN)_diff.{aux,bbl,ps,pdf,dvi,blg,log,out}

view: 
	open $(MAIN).pdf &

