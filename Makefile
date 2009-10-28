.PHONY: doc src

all: src

src:
	$(MAKE) -C $@

dist:
	VERSION=$$(git describe --tags HEAD | sed -e 's/-/./g' | sed -e 's/^v\(.*\)/\1/g'); \
	git archive --format=tar --prefix=bucky-$$VERSION/ HEAD | gzip > bucky-$$VERSION.tgz
