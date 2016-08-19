PDFDIR=../_static

BASENAME=muse_psf_fitting
$(PDFDIR)/$(BASENAME).pdf: diagrams $(BASENAME).tex
	pdflatex $(BASENAME).tex
	pdflatex $(BASENAME).tex
	mv $(BASENAME).pdf $(PDFDIR)/

diagrams: blackman_ift.pdf psf_fft.pdf

blackman_ift.pdf: draw_blackman.py
	./draw_blackman.py

psf_fft.pdf: draw_psf_fft.py
	./draw_psf_fft.py

clean:
	for suffix in aux log out; do /bin/rm -fv $(BASENAME).$$suffix; done
