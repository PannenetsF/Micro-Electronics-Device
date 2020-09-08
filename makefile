all:
	make slide
	make note
	make clean

slide:
	cd slides/ \
	&& pwd \
	&& latexmk -xelatex slides.tex \
	&& cd ..

note:
	cd notes/ \
	&& pwd \
	&& latexmk -xelatex notes.tex \
	&& cd ..

clean:
	rm `find ./ -regex ".*\.log\|.*\.aux\|.*\.xdv\|.*\.in\|.*\.out\|.*\.lua\|.*\.fdb_latexmk\|.*\.fls\|.*\.toc"`





