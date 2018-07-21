prj=""

clean:
	 @rm -rf *.html

# to render a document locally use the following command:
# `make doc prj=folder`, e.g. `make doc prj=gremlin-subgraphs`
doc:
	@cd $(prj) && \
	R -e "rmarkdown::render('$(prj).Rmd', clean = TRUE)"

# to publish a document on Rpubs use the following command:
# `make publish prj=folder`, e.g. `make publish prj=gremlin-subgraphs`
publish: doc
	@cd $(prj) && \
	R -e "title=markdown::yaml_front_matter('$(prj).Rmd')[['title']];rsconnect::rpubsUpload(title, '$(prj).html', '$(prj).Rmd')"
