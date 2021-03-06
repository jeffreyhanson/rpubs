prj=""

clean:
	 @rm -rf *.html

# to render a document locally use the following command:
# `make doc prj=folder`, e.g. `make doc prj=gremlin-subgraphs`
doc:
	@cd $(prj) && \
	R -e "rmarkdown::render('$(prj).Rmd', clean = TRUE)"

# to purl the R code out of a document
# `make purl prj=folder`, e.g. `make purl prj=gremlin-subgraphs`
purl:
	@cd $(prj) && \
	R -e "knitr::purl('$(prj).Rmd')"

# to publish a document on Rpubs use the following command:
# `make publish prj=folder`, e.g. `make publish prj=gremlin-subgraphs`
publish:
	@cd $(prj) && \
	R -e "title=rmarkdown::yaml_front_matter('$(prj).Rmd')[['title']];rsconnect::rpubsUpload(title, '$(prj).html', '$(prj).Rmd')"
