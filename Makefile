PREFIX=/media/jeksterslab/scripts/r
PKG=$(PREFIX)/jeksterslabRds

.PHONY: all clean rm

all : rm
	Rscript -e 'jeksterslabRpkg::pkg_build("$(PKG)", tests2vignettes = FALSE, git = TRUE, github = TRUE, commit_msg = "Automated build")'

clean : rm
	-git add --all
	-git commit -m "[skip ci] Automated clean."
	-git push

rm :
	-Rscript -e 'jeksterslabRpkg::pkg_clean("$(PKG)")'
