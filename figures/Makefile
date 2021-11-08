.PHONY: all
all: recovery-fit.pdf pc1-condition-spec.pdf pc1-subject-spec.pdf

recovery-fit.pdf: s*-*-*-fit.pdf
	pdftk s*-*-*-fit.pdf cat output recovery-fit.pdf

pc1-condition-spec.pdf: s*-*-*-spec.pdf
	pdftk s*-*-*-spec.pdf cat output pc1-condition-spec.pdf

pc1-subject-spec.pdf: s???-spec.pdf
	pdftk s???-spec.pdf cat output pc1-subject-spec.pdf

.PHONY: clean
clean:
	rm -f recovery-fit.pdf pc1-spec.pdf

.PHONY: clean-individual-plots
clean-individual-plots:
	rm -f s*-*-*-*.pdf s???-spec.pdf
