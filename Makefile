.PHONY: doc src

all: src

src:
	$(MAKE) -C $@

dist:
	VERSION=$$(git describe --tags HEAD | sed -e 's/-/./g' | sed -e 's/^v\(.*\)/\1/g'); \
	git archive --prefix=bucky-$$VERSION/ -o bucky-$$VERSION.tgz HEAD
