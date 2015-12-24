LATEX = pdflatex
MAIN  = paper

all: $(MAIN).pdf

$(MAIN).pdf: 	$(MAIN).tex  
	$(LATEX) $(MAIN).tex 
#	bibtex $(MAIN)
#	$(LATEX) $(MAIN).tex 
	$(LATEX) $(MAIN).tex 

diff:	$(MAIN)_diff.pdf

$(MAIN)_diff.pdf: $(MAIN).tex
	latexdiff $(MAIN)_old.tex $(MAIN).tex > $(MAIN)_diff.tex
	$(LATEX) $(MAIN)_diff.tex 
	bibtex $(MAIN)_diff
	$(LATEX) $(MAIN)_diff.tex 
	$(LATEX) $(MAIN)_diff.tex 

tar : 
	tar -czvf $(MAIN).tar.gz $(MAIN).tex mn2e.cls  mn2e.bst conc_dens_vel.pdf concentration.pdf conc_mass_dens.pdf conc_mass_vel.pdf error.pdf

clean:
	rm -f $(MAIN).pdf $(MAIN).out $(MAIN).bbl $(MAIN).aux

view: 
	open $(MAIN).pdf &

